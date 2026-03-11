use std::fmt::Write as _;
use std::fs::File;
use std::io::{self, BufReader, IsTerminal, Write};
use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};
use log::info;
use nu_ansi_term::Color;
use rgb::RGB8;
use textplots::{Chart, ColorPlot, Shape};

use cif2beads::{Bead, ChargeCalc, ChargeResult, coarse_grain};

/// Convert mmCIF protein structures to coarse-grained xyzq representation.
///
/// Each amino acid becomes a single bead at its center of mass.
/// Titratable residues (ASP, GLU, HIS, CYS, TYR, LYS, ARG) get an
/// additional bead at the charge center. Charges are computed using
/// Henderson-Hasselbalch (default) or Monte Carlo titration (--mc).
#[derive(Parser)]
#[command(version, about)]
struct Cli {
    #[command(flatten)]
    common: CommonArgs,

    #[command(subcommand)]
    command: Commands,
}

/// Shared input and condition arguments.
#[derive(Args)]
struct CommonArgs {
    /// Input mmCIF file (reads stdin if omitted).
    input: Option<PathBuf>,

    /// Temperature in Kelvin.
    #[arg(long, default_value = "298.15")]
    temperature: f64,

    /// Ionic strength in mol/L.
    #[arg(long, default_value = "0.1")]
    ionic_strength: f64,

    /// Monte Carlo sweeps per titratable site (0 = Henderson-Hasselbalch only).
    #[arg(long, default_value = "10000")]
    mc: usize,
}

#[derive(Subcommand)]
enum Commands {
    /// Convert mmCIF to coarse-grained xyzq representation.
    Convert {
        /// Output file (writes to stdout if omitted).
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// pH for charge calculation.
        #[arg(long, default_value = "7.0")]
        ph: f64,
    },
    /// Scan average charge, dipole vs pH and plot titration curve.
    Scan {
        /// Start pH.
        #[arg(long, default_value = "3.0")]
        ph_start: f64,

        /// End pH.
        #[arg(long, default_value = "11.0")]
        ph_end: f64,

        /// pH step size.
        #[arg(long, default_value = "0.5")]
        ph_step: f64,

        /// Save titration curve data to file.
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
}

fn read_beads(input: &Option<PathBuf>) -> io::Result<Vec<Bead>> {
    if let Some(path) = input {
        let file = File::open(path)?;
        Ok(coarse_grain(BufReader::new(file)))
    } else if io::stdin().is_terminal() {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "no input file and stdin is a terminal; provide a file or pipe data",
        ))
    } else {
        let stdin = io::stdin();
        Ok(coarse_grain(stdin.lock()))
    }
}

/// Format beads as PQR file.
fn format_pqr(beads: &[Bead], calc: &ChargeCalc) -> String {
    use cif2beads::BeadType;
    let mut out = String::new();
    writeln!(
        out,
        "REMARK cif2beads pH={:.2} T={:.2} I={:.3}",
        calc.ph, calc.temperature, calc.ionic_strength,
    )
    .unwrap();
    for (i, b) in beads.iter().enumerate() {
        let atom_name = match b.bead_type {
            BeadType::Backbone => "CA",
            BeadType::Sidechain => "CB",
            BeadType::Ntr => "N",
            BeadType::Ctr => "O",
            BeadType::Ion => &b.res_name,
        };
        let radius = 2.0;
        writeln!(
            out,
            "{:6}{:5} {:^4.4} {:>3.3} {:1}{:4}    {:8.3}{:8.3}{:8.3}{:8.4}{:7.2}",
            "ATOM",
            i + 1,
            atom_name,
            b.res_name,
            b.chain_id,
            b.res_seq,
            b.x,
            b.y,
            b.z,
            b.charge,
            radius,
        )
        .unwrap();
    }
    writeln!(out, "END").unwrap();
    out
}

/// Format beads as plain XYZ file (no charges).
fn format_xyz(beads: &[Bead], calc: &ChargeCalc) -> String {
    let mut out = String::new();
    writeln!(out, "{}", beads.len()).unwrap();
    writeln!(
        out,
        "cif2beads pH={:.2} T={:.2} I={:.3}",
        calc.ph, calc.temperature, calc.ionic_strength,
    )
    .unwrap();
    for b in beads {
        writeln!(out, "{:<5} {:>10.4} {:>10.4} {:>10.4}", b.res_name, b.x, b.y, b.z).unwrap();
    }
    out
}

/// Generate pH steps over [ph_start, ph_end].
fn ph_steps(ph_start: f64, ph_end: f64, ph_step: f64) -> Vec<f64> {
    let n = ((ph_end - ph_start) / ph_step).round() as usize + 1;
    (0..n)
        .map(|i| (i as f64).mul_add(ph_step, ph_start))
        .collect()
}

fn main() -> io::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();
    let common = &cli.common;

    match cli.command {
        Commands::Convert { output, ph } => {
            let calc = ChargeCalc::new()
                .ph(ph)
                .temperature(common.temperature)
                .ionic_strength(common.ionic_strength)
                .mc(common.mc);

            let beads = read_beads(&common.input)?;
            calc.log_conditions();
            let result = calc.run(&beads);
            let charged = result.apply(&beads);

            info!(
                "⟨Z⟩ = {:.2}, ⟨μ⟩ = {:.1} e·Å",
                result.multipole.charge, result.multipole.dipole
            );

            let is_xyz = output
                .as_ref()
                .and_then(|p| p.extension())
                .is_some_and(|ext| ext == "xyz");
            let text = if is_xyz {
                format_xyz(&charged, &calc)
            } else {
                format_pqr(&charged, &calc)
            };

            if let Some(path) = &output {
                let mut file = File::create(path)?;
                file.write_all(text.as_bytes())?;
            } else {
                io::stdout().write_all(text.as_bytes())?;
            }
        }
        Commands::Scan {
            ph_start,
            ph_end,
            ph_step,
            output,
        } => {
            let beads = read_beads(&common.input)?;
            ChargeCalc::new()
                .temperature(common.temperature)
                .ionic_strength(common.ionic_strength)
                .mc(common.mc)
                .log_conditions();
            let ph_values = ph_steps(ph_start, ph_end, ph_step);

            // HH scan
            let hh_data: Vec<(f64, ChargeResult)> = ph_values
                .iter()
                .map(|&ph| {
                    let result = ChargeCalc::new()
                        .ph(ph)
                        .temperature(common.temperature)
                        .ionic_strength(common.ionic_strength)
                        .run(&beads);
                    (ph, result)
                })
                .collect();

            let hh_f32: Vec<(f32, f32)> = hh_data
                .iter()
                .map(|(ph, r)| (*ph as f32, r.multipole.charge as f32))
                .collect();

            // MC scan (parallelized)
            let mc_data = if common.mc > 0 {
                use rayon::prelude::*;
                let pb = indicatif::ProgressBar::new(ph_values.len() as u64);
                pb.tick(); // force initial draw at 0%
                let data = ph_values
                    .par_iter()
                    .map(|&ph| {
                        let result = ChargeCalc::new()
                            .ph(ph)
                            .temperature(common.temperature)
                            .ionic_strength(common.ionic_strength)
                            .mc(common.mc)
                            .run(&beads);
                        pb.inc(1);
                        (ph, result)
                    })
                    .collect::<Vec<_>>();
                pb.finish_and_clear();
                Some(data)
            } else {
                None
            };

            if let Some(path) = &output {
                let mut file = File::create(path)?;
                if mc_data.is_some() {
                    writeln!(
                        file,
                        "# pH Z(HH) Z2(HH) mu(HH) mu2(HH) Z(MC) Z2(MC) mu(MC) mu2(MC)"
                    )?;
                } else {
                    writeln!(file, "# pH Z(HH) Z2(HH) mu(HH) mu2(HH)")?;
                }
                for (i, (ph, hh)) in hh_data.iter().enumerate() {
                    let m = &hh.multipole;
                    if let Some(mc) = &mc_data {
                        let mc = &mc[i].1.multipole;
                        writeln!(
                            file,
                            "{:.2} {:.4} {:.4} {:.4} {:.4} {:.4} {:.4} {:.4} {:.4}",
                            ph,
                            m.charge,
                            m.charge_sq,
                            m.dipole,
                            m.dipole_sq,
                            mc.charge,
                            mc.charge_sq,
                            mc.dipole,
                            mc.dipole_sq,
                        )?;
                    } else {
                        writeln!(
                            file,
                            "{:.2} {:.4} {:.4} {:.4} {:.4}",
                            ph, m.charge, m.charge_sq, m.dipole, m.dipole_sq,
                        )?;
                    }
                }
                info!("Titration curve saved to {}", path.display());
            }

            const YELLOW: RGB8 = RGB8::new(255, 255, 0);
            const RED: RGB8 = RGB8::new(255, 0, 0);

            if let Some(mc_data) = &mc_data {
                let mc_f32: Vec<(f32, f32)> = mc_data
                    .iter()
                    .map(|(ph, r)| (*ph as f32, r.multipole.charge as f32))
                    .collect();

                info!(
                    "⟨Z⟩ vs pH: {} Henderson-Hasselbalch, {} Monte Carlo ({} sweeps)",
                    Color::Yellow.bold().paint("━━"),
                    Color::Red.bold().paint("━━"),
                    common.mc,
                );

                Chart::new(120, 40, ph_start as f32, ph_end as f32)
                    .linecolorplot(&Shape::Lines(&hh_f32), YELLOW)
                    .linecolorplot(&Shape::Lines(&mc_f32), RED)
                    .nice();
            } else {
                info!(
                    "⟨Z⟩ vs pH: {} Henderson-Hasselbalch",
                    Color::Yellow.bold().paint("━━"),
                );

                Chart::new(120, 40, ph_start as f32, ph_end as f32)
                    .linecolorplot(&Shape::Lines(&hh_f32), YELLOW)
                    .nice();
            }
        }
    }

    Ok(())
}
