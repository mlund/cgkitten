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
    Bead, ChargeCalc, ChargeResult, MultiBead, SingleBead, coarse_grain_pdb_with,
    coarse_grain_with, topology::Topology,
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

#[derive(Args)]
struct ConvertArgs {
    /// Output file (.pqr or .xyz). Defaults to input basename with .pqr extension.
    #[arg(short, long)]
    output: Option<PathBuf>,

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

    /// Charge tolerance for merging titratable site types into a shared topology entry.
    #[arg(long, default_value = "0.02")]
    merge_tol: f64,
}

#[derive(Subcommand)]
enum Commands {
    /// Convert mmCIF to coarse-grained representation.
    Convert(ConvertArgs),
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

/// Derive a default output path from the input path by replacing its extension with `ext`.
/// Falls back to `output.<ext>` when reading from stdin.
fn default_output(input: &Option<PathBuf>, ext: &str) -> PathBuf {
    input
        .as_ref()
        .and_then(|p| p.file_stem())
        .map(|stem| PathBuf::from(stem).with_extension(ext))
        .unwrap_or_else(|| PathBuf::from(format!("output.{ext}")))
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
    debug_assert_eq!(beads.len(), names.len(), "beads and names must be 1:1");
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
            BeadType::Residue | BeadType::Titratable => "CA",
            BeadType::Virtual => "CB",
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
fn format_xyz(beads: &[Bead], names: &[String], cmdline: &str) -> String {
    debug_assert_eq!(beads.len(), names.len(), "beads and names must be 1:1");
    let mut out = String::new();
    writeln!(out, "{}", beads.len()).expect("writing to String is infallible");
    writeln!(out, "{cmdline}").expect("writing to String is infallible");
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

/// Format topology as YAML with unique atom types.
fn format_topology(
    topo: &Topology,
    ff: Option<&dyn cgkitten::forcefield::ForceField>,
    pairs: &[cgkitten::forcefield::PairInteraction],
    cg: &CgPolicy,
) -> String {
    let mut out = format!("# model: {cg}\natoms:\n");
    let mut known: std::collections::HashSet<&str> = std::collections::HashSet::new();
    for t in topo.types() {
        known.insert(t.name);
        let sc_comment =
            if matches!(cg, CgPolicy::Multi) && t.bead_type == cgkitten::BeadType::Virtual {
                " # virtual titratable site"
            } else {
                ""
            };
        let ff_params = ff.and_then(|f| f.params(t.res_name, t.bead_type));
        // Use canonical model mass when available (mass > 0 in the force field).
        // mass == 0.0 is the sentinel for virtual/terminal/ion beads whose mass comes
        // from atomic coordinates instead (ions) or is genuinely zero (virtual sites).
        let mass = ff_params
            .and_then(|p| (p.mass > 0.0).then_some(p.mass))
            .unwrap_or(t.mass);
        let ff_fields = ff_params
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
            t.charge, mass, t.name, ff_fields, sc_comment,
        )
        .expect("writing to String is infallible");
    }

    if let Some(f) = ff {
        out.push_str(f.system_yaml());
    }

    // Pair entries are siblings of `default:` under `nonbonded:`,
    // so Faunus uses these instead of the default mixing rule for
    // the specified type pairs.
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
    assert!(
        ph_step > 0.0 && ph_start <= ph_end,
        "ph_step must be positive and ph_start ≤ ph_end"
    );
    let n = ((ph_end - ph_start) / ph_step).round() as usize + 1;
    (0..n)
        .map(|i| (i as f64).mul_add(ph_step, ph_start))
        .collect()
}

fn print_logo() {
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
}

fn log_chains(beads: &[Bead]) {
    let mut chains: Vec<&str> = beads.iter().map(|b| b.chain_id.as_str()).collect();
    chains.sort_unstable();
    chains.dedup();
    info!("Chains: {} ({})", chains.join(", "), chains.len());
}

fn run_convert(
    common: &CommonArgs,
    policy: &dyn cgkitten::CoarseGrain,
    args: ConvertArgs,
) -> io::Result<()> {
    let ConvertArgs {
        output,
        ph,
        top,
        model,
        scale_hydrophobic,
        merge_tol,
    } = args;
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

    match &common.input {
        Some(path) => info!("Input: {}", path.display()),
        None => info!("Input: stdin"),
    }
    let beads = filter_chains(read_beads(&common.input, policy)?, &common.chain);
    log_chains(&beads);
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

    let cmdline = format!(
        "cgkitten v{} | {}",
        env!("CARGO_PKG_VERSION"),
        std::env::args().collect::<Vec<_>>().join(" ")
    );

    // Assign topology type names first so coordinate files use matching names
    let topo = Topology::new(&charged, merge_tol);
    let names = topo.bead_names();

    // Collect (name, params) for all types that have FF parameters.
    // hydrophobic_pairs() filters this down to hydrophobic residues only,
    // so it's safe to pass everything — non-hydrophobic types are ignored.
    let hp_types: Vec<(String, cgkitten::forcefield::BeadParams)> = topo
        .types()
        .filter_map(|t| {
            ff.as_deref()
                .and_then(|f| f.params(t.res_name, t.bead_type))
                .map(|p| (t.name.to_string(), p))
        })
        .collect();
    let pairs = cgkitten::forcefield::hydrophobic_pairs(
        &hp_types,
        cgkitten::residue::HYDROPHOBIC_RESIDUES,
        &scaling,
    );

    if let Some(path) = output {
        // Explicit output: write only the requested format.
        let text = if path.extension().is_some_and(|e| e == "xyz") {
            format_xyz(&charged, names, &cmdline)
        } else {
            format_pqr(&charged, names, &calc)
        };
        let mut file = File::create(&path)?;
        file.write_all(text.as_bytes())?;
        info!("Output saved to {}", path.display());
    } else {
        // No explicit output: save both PQR and XYZ derived from input basename.
        for ext in ["pqr", "xyz"] {
            let path = default_output(&common.input, ext);
            let text = if ext == "xyz" {
                format_xyz(&charged, names, &cmdline)
            } else {
                format_pqr(&charged, names, &calc)
            };
            let mut file = File::create(&path)?;
            file.write_all(text.as_bytes())?;
            info!("Output saved to {}", path.display());
        }
    }

    let yaml =
        format!("# {cmdline}\n") + &format_topology(&topo, ff.as_deref(), &pairs, &common.cg);
    let mut file = File::create(&top)?;
    file.write_all(yaml.as_bytes())?;
    info!("Topology saved to {}", top.display());

    Ok(())
}

fn run_scan(
    common: &CommonArgs,
    policy: &dyn cgkitten::CoarseGrain,
    ph_start: f64,
    ph_end: f64,
    ph_step: f64,
    output: Option<PathBuf>,
) -> io::Result<()> {
    match &common.input {
        Some(path) => info!("Input: {}", path.display()),
        None => info!("Input: stdin"),
    }
    let beads = filter_chains(read_beads(&common.input, policy)?, &common.chain);
    log_chains(&beads);
    let base_calc = ChargeCalc::new()
        .temperature(common.temperature)
        .ionic_strength(common.ionic_strength)
        .mc(common.mc);
    base_calc.log_conditions();
    let ph_values = ph_steps(ph_start, ph_end, ph_step);

    // HH scan (always Henderson-Hasselbalch, ignoring mc setting)
    let hh_data: Vec<(f64, ChargeResult)> = ph_values
        .iter()
        .map(|&ph| {
            let result = base_calc.clone().ph(ph).mc(0).run(&beads);
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
                let result = base_calc.clone().ph(ph).run(&beads);
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

    Ok(())
}

fn main() -> io::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    print_logo();

    let cli = Cli::parse();
    let common = &cli.common;
    let policy = cg_policy(&common.cg);

    match cli.command {
        Commands::Convert(args) => run_convert(common, policy, args),
        Commands::Scan {
            ph_start,
            ph_end,
            ph_step,
            output,
        } => run_scan(common, policy, ph_start, ph_end, ph_step, output),
    }
}
