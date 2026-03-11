# cif2beads

Convert mmCIF protein structures to a coarse-grained representation.

Each amino acid residue becomes a single bead at its geometric center.
Titratable residues (ASP, GLU, HIS, CYS, TYR, LYS, ARG) get an additional bead
at the charge center. N- and C-terminal charges are always included.
Water and non-protein molecules are removed; coordinated metal ions are retained.
Multi-chain structures are fully supported with per-chain terminal beads.

Charges are computed using Metropolis Monte Carlo titration (default, 10 000 sweeps)
with screened Coulomb (Yukawa) electrostatics. Use `--mc 0` for Henderson-Hasselbalch only.

## Install

```bash
cargo install --git https://github.com/mlund/cif2beads
```

## Usage

Shared flags (`input`, `--temperature`, `--ionic-strength`, `--mc`) are placed
before the subcommand.

```bash
# Convert at pH 7 with MC titration (default: 10000 sweeps)
cif2beads structure.cif convert

# Henderson-Hasselbalch only (no MC)
cif2beads structure.cif --mc 0 convert

# Custom pH, sweeps, and output file
cif2beads structure.cif --mc 50000 convert --ph 4.5 -o output.pqr

# Pipe from stdin
cat structure.cif | cif2beads convert -o output.pqr

# Plain XYZ output (no charges)
cif2beads structure.cif convert -o output.xyz

# Custom conditions
cif2beads structure.cif --temperature 310 --ionic-strength 0.15 convert --ph 4.5

# pH scan with terminal plot
cif2beads structure.cif scan

# pH scan with HH only
cif2beads structure.cif --mc 0 scan

# pH scan with custom range and save to file
cif2beads structure.cif scan --ph-start 2 --ph-end 12 --ph-step 0.25 -o curve.dat
```

### Output formats

The output format is determined by the `-o` filename extension:

- **`.pqr`** (default, also used for stdout) —
  [PQR format](https://userguide.mdanalysis.org/stable/formats/reference/pqr.html),
  easily visualized in VMD or PyMOL. Includes charges and radii.
- **`.xyz`** — Plain XYZ format (coordinates only, no charges).

```
REMARK cif2beads pH=7.00 T=298.15 I=0.100
ATOM      1  CA  MET A   1      25.906  25.079   4.445  0.0000   2.00
ATOM      2  CB  ASP A   2      27.340  24.430   2.614 -0.9994   2.00
ATOM      3  CA  ZN  A   3       5.000   5.000   5.000  2.0000   2.00
END
```

## Monte Carlo titration

Charges are computed using Metropolis Monte Carlo titration by default
(`--mc 10000`), which captures many-body electrostatic coupling between
titratable sites. Use `--mc 0` for the faster Henderson-Hasselbalch
independent-site approximation.

### Method

1. **Initialization**: Protonation states are set from Henderson-Hasselbalch
   at the given pH, providing a physically reasonable starting configuration.

2. **Yukawa electrostatics**: Pairwise interactions use a screened Coulomb
   (Yukawa) potential:

   U/kT = λ\_B · q\_i · Σ\_j z\_j · exp(−r\_ij / λ\_D) / r\_ij

   where λ\_B is the Bjerrum length (~7.1 Å at 298 K in water) and λ\_D is
   the Debye screening length (~9.6 Å at 0.1 M ionic strength). The kernel
   exp(−r/λ\_D)/r is precomputed once for all bead pairs.

3. **Metropolis sweeps**: Each sweep proposes N random protonation state
   changes (one per titratable site on average). A trial move toggles a
   site's protonation and is accepted with probability:

   P = min(1, exp(−ΔU/kT))

   where ΔU/kT = Δq · φ\_i ± ln(10) · (pH − pKa), with φ\_i being the
   electrostatic potential at site i from all other charges.

4. **Ensemble averages**: The mean charge ⟨Z⟩, ⟨Z²⟩, dipole moment ⟨μ⟩,
   and ⟨μ²⟩ are accumulated over all sweeps as true ensemble averages.

## pH scan

The `scan` subcommand plots average net charge ⟨Z⟩ as a function of pH in
the terminal. Henderson-Hasselbalch is always shown; when MC is enabled
(the default), the Monte Carlo result is overlaid in a different color.
MC pH points are computed in parallel using rayon.

With `-o`, scan data is saved to a space-separated file with columns:

```
# pH Z(HH) Z2(HH) mu(HH) mu2(HH) [Z(MC) Z2(MC) mu(MC) mu2(MC)]
```

where Z is net charge, Z2 is ⟨Z²⟩, mu is dipole moment ⟨μ⟩ (e·Å), and
mu2 is ⟨μ²⟩. MC columns are included when `--mc` > 0.

## Library usage

```rust
use cif2beads::{ChargeCalc, coarse_grain};

let cif_data = std::fs::read("structure.cif").unwrap();
let beads = coarse_grain(cif_data.as_slice());

// Henderson-Hasselbalch (default)
let result = ChargeCalc::new().ph(7.4).run(&beads);
println!("⟨Z⟩ = {:.2}, ⟨μ⟩ = {:.1} e·Å", result.multipole.charge, result.multipole.dipole);

// Monte Carlo titration
let result = ChargeCalc::new().ph(7.4).mc(10000).run(&beads);
let charged_beads = result.apply(&beads);

// Serde support (with "serde" feature)
// let calc: ChargeCalc = serde_json::from_str(r#"{"ph": 7.4, "mc": 10000}"#)?;
```

Use without the CLI:

```toml
[dependencies]
cif2beads = { git = "https://github.com/mlund/cif2beads", default-features = false }
```

Enable serde support:

```toml
[dependencies]
cif2beads = { git = "https://github.com/mlund/cif2beads", default-features = false, features = ["serde"] }
```

## License

Apache-2.0
