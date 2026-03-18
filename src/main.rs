use std::fmt::Write as _;
use std::fs::File;
use std::io::{self, BufReader, IsTerminal, Write};
use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};
use log::info;
use nu_ansi_term::Color;
use rgb::RGB8;
use textplots::{Chart, ColorPlot, Shape};

use cgkitten::{
    Bead, ChargeCalc, ChargeResult, MultiBead, SingleBead, coarse_grain_pdb_with, coarse_grain_with,
};

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

    /// Coarse-graining policy: "multi" (backbone + sidechain beads) or "single" (one bead per residue).
    #[arg(long, default_value = "multi")]
    cg: CgPolicy,

    /// Include only these chain IDs (default: all chains). Repeat to include multiple: --chain A --chain B
    #[arg(long)]
    chain: Vec<String>,
}

#[derive(Clone, clap::ValueEnum)]
enum CgPolicy {
    Multi,
    Single,
}

impl std::fmt::Display for CgPolicy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CgPolicy::Multi => f.write_str("multi"),
            CgPolicy::Single => f.write_str("single"),
        }
    }
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

        /// Scale hydrophobic pair interactions: lambda:<factor> or epsilon:<factor>.
        #[arg(long)]
        scale_hydrophobic: Option<String>,
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

fn read_beads(
    input: &Option<PathBuf>,
    policy: &dyn cgkitten::CoarseGrain,
) -> io::Result<Vec<Bead>> {
    if let Some(path) = input {
        let file = File::open(path)?;
        // Dispatch on file extension; stdin always assumes mmCIF.
        let is_pdb = path
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("pdb"));
        if is_pdb {
            Ok(coarse_grain_pdb_with(BufReader::new(file), policy))
        } else {
            Ok(coarse_grain_with(BufReader::new(file), policy))
        }
    } else if io::stdin().is_terminal() {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "no input file and stdin is a terminal; provide a file or pipe data",
        ))
    } else {
        let stdin = io::stdin();
        Ok(coarse_grain_with(stdin.lock(), policy))
    }
}

/// Filter beads to the requested chains; if `chains` is empty all beads are kept.
fn filter_chains(beads: Vec<Bead>, chains: &[String]) -> Vec<Bead> {
    if chains.is_empty() {
        return beads;
    }
    let kept: Vec<_> = beads
        .into_iter()
        .filter(|b| chains.iter().any(|c| c == &b.chain_id))
        .collect();
    info!("Chain filter {:?}: {} beads retained", chains, kept.len());
    kept
}

fn cg_policy(p: &CgPolicy) -> &'static dyn cgkitten::CoarseGrain {
    match p {
        CgPolicy::Multi => &MultiBead,
        CgPolicy::Single => &SingleBead,
    }
}

/// Format beads as PQR file.
fn format_pqr(beads: &[Bead], names: &[String], calc: &ChargeCalc) -> String {
    use cgkitten::BeadType;
    let mut out = String::new();
    writeln!(
        out,
        "REMARK cif2top pH={:.2} T={:.2} I={:.3}",
        calc.ph, calc.temperature, calc.ionic_strength,
    )
    .expect("writing to String is infallible");
    for (i, b) in beads.iter().enumerate() {
        let atom_name = match b.bead_type {
            BeadType::Backbone | BeadType::Titratable => "CA",
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
            &names[i],
            b.chain_id,
            b.res_seq,
            b.x,
            b.y,
            b.z,
            b.charge,
            radius,
        )
        .expect("writing to String is infallible");
    }
    writeln!(out, "END").expect("writing to String is infallible");
    out
}

/// Format beads as plain XYZ file (no charges).
fn format_xyz(beads: &[Bead], names: &[String], calc: &ChargeCalc) -> String {
    let mut out = String::new();
    writeln!(out, "{}", beads.len()).expect("writing to String is infallible");
    writeln!(
        out,
        "cif2top pH={:.2} T={:.2} I={:.3}",
        calc.ph, calc.temperature, calc.ionic_strength,
    )
    .expect("writing to String is infallible");
    for (i, b) in beads.iter().enumerate() {
        writeln!(
            out,
            "{:<5} {:>10.4} {:>10.4} {:>10.4}",
            &names[i], b.x, b.y, b.z
        )
        .expect("writing to String is infallible");
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
    bead_type: cgkitten::BeadType,
}

/// Assign each bead to a topology type and return per-bead type names.
///
/// Backbone beads use their residue name (e.g., ALA, MET).
/// Titratable sites with the same element, mass, and charge (within 1%)
/// are merged into a single type; otherwise they get unique numbered
/// names (e.g., O1, O2, N1).
fn assign_types(beads: &[Bead]) -> (Vec<AtomType>, Vec<String>) {
    use cgkitten::BeadType;

    let mut types: Vec<AtomType> = Vec::new();
    let mut counters: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    // Per-bead type name, matching the topology entry
    let mut bead_type_names: Vec<String> = Vec::with_capacity(beads.len());

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
                bead_type_names.push(b.res_name.clone());
            }
            BeadType::Sidechain | BeadType::Ntr | BeadType::Ctr | BeadType::Titratable => {
                // Merge sites with same element/residue and charge (within tolerance).
                // Titratable (single-bead whole-residue) ignores mass like Backbone does,
                // since atom-count variations across instances are not physically meaningful.
                // Sidechain/Ntr/Ctr also require mass to match, as they are sub-residue beads.
                let existing = types.iter().find(|t| {
                    t.res_name == b.res_name
                        && (b.bead_type == BeadType::Titratable || t.mass == b.mass)
                        && (t.charge - b.charge).abs() < CHARGE_MERGE_TOL
                        && t.bead_type.is_titratable()
                });
                if let Some(t) = existing {
                    bead_type_names.push(t.name.clone());
                } else {
                    let counter = counters.entry(b.res_name.clone()).or_insert(0);
                    *counter += 1;
                    let name = format!("{}{}", b.res_name, counter);
                    bead_type_names.push(name.clone());
                    types.push(AtomType {
                        name,
                        charge: b.charge,
                        mass: b.mass,
                        res_name: b.res_name.clone(),
                        bead_type: b.bead_type,
                    });
                }
            }
        }
    }

    // Strip trailing "1" from names that have only a single variant (e.g. "LYS1" -> "LYS")
    for (res_name, &count) in &counters {
        if count == 1 {
            let old = format!("{res_name}1");
            let new = res_name.clone();
            if let Some(t) = types.iter_mut().find(|t| t.name == old) {
                t.name = new.clone();
            }
            for n in bead_type_names.iter_mut() {
                if *n == old {
                    *n = new.clone();
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

    (types, bead_type_names)
}

/// Format topology as YAML with unique atom types.
fn format_topology(
    types: &[AtomType],
    ff: Option<&dyn cgkitten::forcefield::ForceField>,
    pairs: &[cgkitten::forcefield::PairInteraction],
    cg: &CgPolicy,
) -> String {
    let mut out = format!("# model: {cg}\natoms:\n");
    for t in types {
        let sc_comment =
            if matches!(cg, CgPolicy::Multi) && t.bead_type == cgkitten::BeadType::Sidechain {
                " # titratable sidechain"
            } else {
                ""
            };
        let ff_fields = ff
            .and_then(|f| f.params(&t.res_name, t.bead_type))
            .map(|p| {
                format!(
                    ", σ: {}, ε: {}, hydrophobicity: !Lambda {}",
                    p.sigma, p.epsilon, p.lambda
                )
            })
            .unwrap_or_default();
        writeln!(
            out,
            "  - {{charge: {:.4}, mass: {:.2}, name: {}{}}}{}",
            t.charge, t.mass, t.name, ff_fields, sc_comment,
        )
        .expect("writing to String is infallible");
    }

    if let Some(f) = ff {
        out.push_str(f.system_yaml());
    }

    // Pair entries are siblings of `default:` under `nonbonded:`,
    // so Faunus uses these instead of the default mixing rule for
    // the specified type pairs.
    let known: std::collections::HashSet<&str> = types.iter().map(|t| t.name.as_str()).collect();
    for pair in pairs {
        for name in [pair.name_a.as_str(), pair.name_b.as_str()] {
            if !known.contains(name) {
                log::warn!(
                    "Pair references unknown atom type '{name}' — Faunus will silently ignore this entry"
                );
            }
        }
        out.push_str(&cgkitten::forcefield::format_pair_yaml(pair));
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

    if io::stderr().is_terminal() {
        eprintln!(
            " /\\_/\\\n{}~ {}\n         {} v{}\n",
            Color::Yellow.bold().paint("(=^·^=)"),
            Color::Cyan.paint("○-○-○-○-○"),
            Color::Green.bold().paint("cgkitten"),
            env!("CARGO_PKG_VERSION"),
        );
    } else {
        eprintln!(" /\\_/\\");
        eprintln!("(=^·^=)~ ○-○-○-○-○");
        eprintln!("         cgkitten v{}\n", env!("CARGO_PKG_VERSION"));
    }

    let cli = Cli::parse();
    let common = &cli.common;
    let policy = cg_policy(&common.cg);

    match cli.command {
        Commands::Convert {
            output,
            ph,
            top,
            model,
            scale_hydrophobic,
        } => {
            let calc = ChargeCalc::new()
                .ph(ph)
                .temperature(common.temperature)
                .ionic_strength(common.ionic_strength)
                .mc(common.mc);

            let ff = cgkitten::forcefield::from_name(&model);
            // "none" is an explicit opt-out; unknown names are errors
            if ff.is_none() && model != "none" {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("unknown force field model: {model}"),
                ));
            }

            let beads = filter_chains(read_beads(&common.input, policy)?, &common.chain);
            match &common.input {
                Some(path) => info!("Input: {}", path.display()),
                None => info!("Input: stdin"),
            }
            calc.log_conditions();
            let result = calc.run(&beads);
            let charged = result.apply(&beads);

            info!(
                "⟨Z⟩ = {:.2}, ⟨μ⟩ = {:.1} e·Å",
                result.multipole.charge, result.multipole.dipole
            );

            let scaling: cgkitten::forcefield::HydrophobicScaling = scale_hydrophobic
                .as_deref()
                .map(|s| {
                    s.parse()
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
                })
                .transpose()?
                .unwrap_or_default();

            if common.mc > 0 {
                info!("Titration: MC ({} steps)", common.mc);
            } else {
                info!("Titration: Henderson-Hasselbalch");
            }
            info!("Hydrophobic scaling: {scaling}");

            // Assign topology type names first so coordinate files use matching names
            let (types, names) = assign_types(&charged);

            // Collect (name, params) for all types that have FF parameters.
            // hydrophobic_pairs() filters this down to hydrophobic residues only,
            // so it's safe to pass everything — non-hydrophobic types are ignored.
            let hp_types: Vec<(String, cgkitten::forcefield::BeadParams)> = types
                .iter()
                .filter_map(|t| {
                    ff.as_deref()
                        .and_then(|f| f.params(&t.res_name, t.bead_type))
                        .map(|p| (t.name.clone(), p))
                })
                .collect();
            let pairs = cgkitten::forcefield::hydrophobic_pairs(
                &hp_types,
                cgkitten::residue::HYDROPHOBIC_RESIDUES,
                &scaling,
            );

            let is_xyz = output.extension().is_some_and(|ext| ext == "xyz");
            let text = if is_xyz {
                format_xyz(&charged, &names, &calc)
            } else {
                format_pqr(&charged, &names, &calc)
            };

            let mut file = File::create(&output)?;
            file.write_all(text.as_bytes())?;

            {
                let mut types = types;
                types.sort_by(|a, b| a.name.cmp(&b.name));
                let yaml = format_topology(&types, ff.as_deref(), &pairs, &common.cg);
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
            let beads = filter_chains(read_beads(&common.input, policy)?, &common.chain);
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
