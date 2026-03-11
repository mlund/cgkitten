use std::fmt::Write as _;
use std::fs::File;
use std::io::{self, BufReader, IsTerminal, Write};
use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};
use log::info;
use nu_ansi_term::Color;
use rgb::RGB8;
use textplots::{Chart, ColorPlot, Shape};

use cif2top::{Bead, ChargeCalc, ChargeResult, coarse_grain};

/// Convert mmCIF protein structures to coarse-grained representation.
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
    /// Convert mmCIF to coarse-grained representation.
    Convert {
        /// Output file (.pqr or .xyz).
        #[arg(short, long)]
        output: PathBuf,

        /// pH for charge calculation.
        #[arg(long, default_value = "7.0")]
        ph: f64,

        /// Save topology as YAML file.
        #[arg(long, default_value = "topology.yaml")]
        top: PathBuf,

        /// Force field model for topology parameters.
        #[arg(long, default_value = "calvados3")]
        model: String,
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
    use cif2top::BeadType;
    let mut out = String::new();
    writeln!(
        out,
        "REMARK cif2top pH={:.2} T={:.2} I={:.3}",
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
        let radius = 2.0; // placeholder; PQR format requires a radius column
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
        "cif2top pH={:.2} T={:.2} I={:.3}",
        calc.ph, calc.temperature, calc.ionic_strength,
    )
    .unwrap();
    for b in beads {
        writeln!(
            out,
            "{:<5} {:>10.4} {:>10.4} {:>10.4}",
            b.res_name, b.x, b.y, b.z
        )
        .unwrap();
    }
    out
}

/// Charge tolerance for merging titratable site types (1%).
const CHARGE_MERGE_TOL: f64 = 0.01;

/// A registered atom type for topology output.
struct AtomType {
    name: String,
    charge: f64,
    mass: f64,
    res_name: String,
    bead_type: cif2top::BeadType,
}

/// Format topology as YAML with unique atom types.
///
/// Backbone beads use their residue name (e.g., ALA, MET).
/// Titratable sites with the same element, mass, and charge (within 1%)
/// are merged into a single type; otherwise they get unique numbered
/// names (e.g., O1, O2, N1).
fn format_topology(beads: &[Bead], ff: Option<&dyn cif2top::forcefield::ForceField>) -> String {
    use cif2top::BeadType;

    let mut types: Vec<AtomType> = Vec::new();
    let mut counters: std::collections::HashMap<String, usize> = std::collections::HashMap::new();

    // Assign each bead to an existing or new atom type
    for b in beads {
        match b.bead_type {
            BeadType::Backbone | BeadType::Ion => {
                // All beads of the same residue type share identical FF parameters
                if !types.iter().any(|t| t.name == b.res_name) {
                    types.push(AtomType {
                        name: b.res_name.clone(),
                        charge: b.charge,
                        mass: b.mass,
                        res_name: b.res_name.clone(),
                        bead_type: b.bead_type,
                    });
                }
            }
            BeadType::Sidechain | BeadType::Ntr | BeadType::Ctr => {
                // Merge sites with same element, mass, and charge (within tolerance)
                // to reduce topology size without losing physical accuracy
                let merged = types.iter().any(|t| {
                    t.res_name == b.res_name
                        && t.mass == b.mass
                        && (t.charge - b.charge).abs() < CHARGE_MERGE_TOL
                        && t.bead_type.is_titratable()
                });
                if !merged {
                    let counter = counters.entry(b.res_name.clone()).or_insert(0);
                    *counter += 1;
                    types.push(AtomType {
                        name: format!("{}{}", b.res_name, counter),
                        charge: b.charge,
                        mass: b.mass,
                        res_name: b.res_name.clone(),
                        bead_type: b.bead_type,
                    });
                }
            }
        }
    }

    let n_before = beads.iter().filter(|b| b.bead_type.is_titratable()).count();
    let n_after = types.iter().filter(|t| t.bead_type.is_titratable()).count();
    if n_after < n_before {
        info!(
            "Merged {n_before} titratable sites into {n_after} unique types (tolerance {:.0}%)",
            CHARGE_MERGE_TOL * 100.0
        );
    }

    let mut out = String::from("atoms:\n");
    for t in &types {
        if let Some(params) = ff.and_then(|f| f.params(&t.res_name, t.bead_type)) {
            writeln!(
                out,
                "  - {{charge: {:.4}, mass: {:.2}, name: {}, σ: {}, ε: {}, hydrophobicity: !Lambda {}}}",
                t.charge, t.mass, t.name, params.sigma, params.epsilon, params.lambda,
            )
            .unwrap();
        } else {
            writeln!(
                out,
                "  - {{charge: {:.4}, mass: {:.2}, name: {}}}",
                t.charge, t.mass, t.name,
            )
            .unwrap();
        }
    }

    if let Some(f) = ff {
        out.push_str(f.system_yaml());
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
        Commands::Convert {
            output,
            ph,
            top,
            model,
        } => {
            let calc = ChargeCalc::new()
                .ph(ph)
                .temperature(common.temperature)
                .ionic_strength(common.ionic_strength)
                .mc(common.mc);

            let ff = cif2top::forcefield::from_name(&model);
            // "none" is an explicit opt-out; unknown names are errors
            if ff.is_none() && model != "none" {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("unknown force field model: {model}"),
                ));
            }

            let beads = read_beads(&common.input)?;
            calc.log_conditions();
            let result = calc.run(&beads);
            let charged = result.apply(&beads);

            info!(
                "⟨Z⟩ = {:.2}, ⟨μ⟩ = {:.1} e·Å",
                result.multipole.charge, result.multipole.dipole
            );

            let is_xyz = output
                .extension()
                .is_some_and(|ext| ext == "xyz");
            let text = if is_xyz {
                format_xyz(&charged, &calc)
            } else {
                format_pqr(&charged, &calc)
            };

            let mut file = File::create(&output)?;
            file.write_all(text.as_bytes())?;

            {
                let yaml = format_topology(&charged, ff.as_deref());
                let mut file = File::create(&top)?;
                file.write_all(yaml.as_bytes())?;
                info!("Topology saved to {}", top.display());
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
                pb.tick(); // force initial draw before rayon dispatches work
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
