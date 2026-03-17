//! Convert mmCIF protein structures to coarse-grained xyzq representation.
//!
//! Each amino acid residue is represented as a single bead at its center of mass.
//! Titratable side-chains (ASP, GLU, HIS, CYS, TYR, LYS, ARG) generate an
//! additional bead at the charge center. Charges are computed via
//! [`ChargeCalc`], which supports Henderson-Hasselbalch and Monte Carlo titration.
//!
//! ```no_run
//! let cif_data = std::fs::read("structure.cif").unwrap();
//! let beads = cif2top::coarse_grain(cif_data.as_slice());
//! let result = cif2top::ChargeCalc::new()
//!     .ph(7.0)
//!     .mc(10000)
//!     .run(&beads);
//! let charged = result.apply(&beads);
//! ```

#![warn(missing_docs)]

pub(crate) mod charge;
/// Force field model definitions.
pub mod forcefield;
pub(crate) mod mc;
pub(crate) mod mmcif;
/// Amino acid residue definitions and constants.
pub mod residue;

use std::collections::{HashMap, HashSet};
use std::io::BufRead;

use log::{debug, info};

use charge::{ChargeCalculator, Conditions, HendersonHasselbalch};
use mmcif::{AtomRecord, DisulfideBond, ParseOptions};
use residue::{CTERM_GROUP, NTERM_GROUP, TitratableGroup};

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// A coarse-grained bead with position and charge.
#[derive(Clone, Debug)]
pub struct Bead {
    /// X coordinate in Ångströms.
    pub x: f64,
    /// Y coordinate in Ångströms.
    pub y: f64,
    /// Z coordinate in Ångströms.
    pub z: f64,
    /// Partial charge in elementary charge units.
    pub charge: f64,
    /// Mass in Da (sum of constituent atom masses).
    pub mass: f64,
    /// Source residue name (e.g., "ALA", "ZN").
    pub res_name: String,
    /// Chain identifier.
    pub chain_id: String,
    /// Residue sequence number.
    pub res_seq: i32,
    /// Bead classification.
    pub bead_type: BeadType,
}

/// Classification of a coarse-grained bead.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BeadType {
    /// Center-of-mass bead for the residue backbone + side-chain.
    Backbone,
    /// Extra bead for a titratable side-chain group.
    Sidechain,
    /// N-terminal group.
    Ntr,
    /// C-terminal group.
    Ctr,
    /// Metal ion.
    Ion,
    /// Single-bead residue that is also titratable (used by SingleBead policy).
    Titratable,
}

impl BeadType {
    /// True for bead types that carry a titratable charge (sidechain, N-/C-terminal, titratable).
    pub fn is_titratable(self) -> bool {
        matches!(
            self,
            Self::Sidechain | Self::Ntr | Self::Ctr | Self::Titratable
        )
    }
}

impl std::fmt::Display for BeadType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Backbone => write!(f, "backbone"),
            Self::Sidechain => write!(f, "sidechain"),
            Self::Ntr => write!(f, "NTR"),
            Self::Ctr => write!(f, "CTR"),
            Self::Ion => write!(f, "ion"),
            Self::Titratable => write!(f, "titratable"),
        }
    }
}

/// Charge calculation builder.
///
/// Configures pH, temperature, ionic strength, and charge method (Henderson-Hasselbalch
/// or Monte Carlo). Use [`run`](ChargeCalc::run) to compute per-bead charges and multipole moments.
///
/// ```ignore
/// // Builder pattern
/// let result = ChargeCalc::new().ph(7.4).mc(10000).run(&beads);
///
/// // Direct construction
/// let result = ChargeCalc { ph: 7.4, mc: 10000, ..Default::default() }.run(&beads);
///
/// // Serde (with "serde" feature)
/// let calc: ChargeCalc = serde_json::from_str(r#"{"ph": 7.4, "mc": 10000}"#)?;
/// ```
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
#[cfg_attr(feature = "serde", serde(default))]
pub struct ChargeCalc {
    /// pH for charge calculation (default 7.0).
    pub ph: f64,
    /// Temperature in Kelvin (default 298.15).
    pub temperature: f64,
    /// Ionic strength in mol/L (default 0.1).
    pub ionic_strength: f64,
    /// Monte Carlo sweeps per titratable site (0 = Henderson-Hasselbalch only).
    pub mc: usize,
}

impl Default for ChargeCalc {
    fn default() -> Self {
        Self {
            ph: 7.0,
            temperature: 298.15,
            ionic_strength: 0.1,
            mc: 0,
        }
    }
}

impl ChargeCalc {
    /// Create a new charge calculator with default settings (pH 7.0, HH).
    pub fn new() -> Self {
        Self::default()
    }

    /// Set pH.
    pub fn ph(mut self, ph: f64) -> Self {
        self.ph = ph;
        self
    }

    /// Set temperature in Kelvin.
    pub fn temperature(mut self, t: f64) -> Self {
        self.temperature = t;
        self
    }

    /// Set ionic strength in mol/L.
    pub fn ionic_strength(mut self, i: f64) -> Self {
        self.ionic_strength = i;
        self
    }

    /// Set Monte Carlo sweeps (0 = Henderson-Hasselbalch only).
    pub fn mc(mut self, sweeps: usize) -> Self {
        self.mc = sweeps;
        self
    }

    /// Log physical conditions (pH, T, I, Bjerrum/Debye lengths).
    pub fn log_conditions(&self) {
        let lb = mc::bjerrum_length(self.temperature);
        let ld = mc::debye_length(self.temperature, self.ionic_strength);
        info!(
            "pH={:.2} T={:.1} K I={:.3} M λ_B={:.2} Å λ_D={:.2} Å",
            self.ph, self.temperature, self.ionic_strength, lb, ld,
        );
    }

    /// Compute per-bead charges and multipole moments.
    pub fn run(&self, beads: &[Bead]) -> ChargeResult {
        let conditions = Conditions::from(self);

        if self.mc > 0 {
            let mut engine = mc::TitrationMC::new(beads, &conditions);
            let n_sites = engine.num_sites();
            debug!(
                "MC titration: {} sweeps × {} sites = {} trial moves",
                self.mc,
                n_sites,
                self.mc * n_sites,
            );
            let mut rng = rand::thread_rng();
            let mc_result = engine.run(self.mc, &mut rng);
            ChargeResult {
                multipole: Multipole {
                    charge: mc_result.mean_charge,
                    charge_sq: mc_result.mean_charge_sq,
                    dipole: mc_result.mean_dipole,
                    dipole_sq: mc_result.mean_dipole_sq,
                },
                charges: mc_result.charges,
            }
        } else {
            let hh = HendersonHasselbalch;
            let rn = res_name_map(beads);
            // Start from existing charges so non-titratable beads (e.g. metal ions) keep theirs
            let mut charges: Vec<f64> = beads.iter().map(|b| b.charge).collect();
            for (i, b) in beads.iter().enumerate() {
                if let Some(group) = titratable_group_for_bead(b, &rn) {
                    charges[i] = hh.charge(group, &conditions);
                }
            }
            let multipole = compute_multipole(beads, &charges);
            ChargeResult { multipole, charges }
        }
    }
}

/// Multipole moments: net charge and dipole with their squares.
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde", derive(serde::Deserialize, serde::Serialize))]
pub struct Multipole {
    /// Net charge ⟨Z⟩.
    pub charge: f64,
    /// ⟨Z²⟩ (equals Z² for HH; ensemble average for MC).
    pub charge_sq: f64,
    /// Dipole moment magnitude ⟨μ⟩ in e·Å.
    pub dipole: f64,
    /// ⟨μ²⟩ in (e·Å)².
    pub dipole_sq: f64,
}

/// Result of a charge calculation.
#[derive(Clone, Debug)]
pub struct ChargeResult {
    /// Multipole moments.
    pub multipole: Multipole,
    /// Per-bead charges (same indexing as input beads).
    pub charges: Vec<f64>,
}

impl ChargeResult {
    /// Apply computed charges to beads, returning new beads.
    pub fn apply(&self, beads: &[Bead]) -> Vec<Bead> {
        beads
            .iter()
            .enumerate()
            .map(|(i, b)| Bead {
                charge: self.charges[i],
                ..b.clone()
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Coarse-graining trait and implementations
// ---------------------------------------------------------------------------

/// Convert one amino acid residue into beads (backbone/sidechain only, not terminals).
pub trait CoarseGrain {
    /// Convert one amino acid residue into beads.
    fn residue_to_beads(
        &self,
        key: &ResidueKey,
        atoms: &[&AtomRecord],
        is_ss_bonded: bool,
    ) -> Vec<Bead>;
}

/// Look up the titratable group for a residue, returning `None` for SS-bonded CYS.
fn titratable_group_unless_ss(
    key: &ResidueKey,
    res_name: &str,
    is_ss_bonded: bool,
) -> Option<&'static residue::TitratableGroup> {
    if res_name == "CYS" && is_ss_bonded {
        debug!(
            "Skipping titration for SS-bonded CYS {}:{}",
            key.chain_id, key.res_seq
        );
        return None;
    }
    residue::find_titratable_group(res_name)
}

/// Create a residue bead at the geometric center of the given atoms.
fn make_residue_bead(key: &ResidueKey, atoms: &[&AtomRecord], bead_type: BeadType) -> Bead {
    let ((cx, cy, cz), mass) = center_and_mass(atoms).unwrap();
    Bead {
        x: cx,
        y: cy,
        z: cz,
        charge: 0.0,
        mass,
        res_name: atoms[0].res_name.clone(),
        chain_id: key.chain_id.clone(),
        res_seq: key.res_seq,
        bead_type,
    }
}

/// Multi-bead coarse-graining: one backbone bead + optional sidechain bead per titratable residue.
pub struct MultiBead;

impl CoarseGrain for MultiBead {
    fn residue_to_beads(
        &self,
        key: &ResidueKey,
        atoms: &[&AtomRecord],
        is_ss_bonded: bool,
    ) -> Vec<Bead> {
        let mut beads = vec![make_residue_bead(key, atoms, BeadType::Backbone)];

        if let Some(group) = titratable_group_unless_ss(key, &atoms[0].res_name, is_ss_bonded)
            && let Some(bead) = make_titratable_bead(key, atoms, group, BeadType::Sidechain)
        {
            beads.push(bead);
        }

        beads
    }
}

/// Single-bead coarse-graining: one bead per residue. Titratable residues get `BeadType::Titratable`.
pub struct SingleBead;

impl CoarseGrain for SingleBead {
    fn residue_to_beads(
        &self,
        key: &ResidueKey,
        atoms: &[&AtomRecord],
        is_ss_bonded: bool,
    ) -> Vec<Bead> {
        let bead_type =
            if titratable_group_unless_ss(key, &atoms[0].res_name, is_ss_bonded).is_some() {
                BeadType::Titratable
            } else {
                BeadType::Backbone
            };

        vec![make_residue_bead(key, atoms, bead_type)]
    }
}

// ---------------------------------------------------------------------------
// Public functions
// ---------------------------------------------------------------------------

/// Convert mmCIF data to coarse-grained beads using the default multi-bead policy.
///
/// Reads atom records, groups by residue, and creates backbone and titratable
/// side-chain beads. Charges are set to 0.0 for titratable sites; use
/// [`ChargeCalc::run`] to compute charges.
///
/// Metal ions retain their formal charges.
pub fn coarse_grain<R: BufRead>(reader: R) -> Vec<Bead> {
    coarse_grain_with(reader, &MultiBead)
}

/// Convert mmCIF data to coarse-grained beads using the given policy.
pub fn coarse_grain_with<R: BufRead>(reader: R, policy: &dyn CoarseGrain) -> Vec<Bead> {
    let parse_options = ParseOptions::default();
    let parsed = mmcif::parse_mmcif(reader, &parse_options);
    debug!(
        "Parsed {} atom records, {} disulfide bonds from _struct_conn",
        parsed.atoms.len(),
        parsed.disulfide_bonds.len(),
    );

    records_to_beads(&parsed.atoms, &parsed.disulfide_bonds, policy)
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Unique residue key for grouping atoms.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct ResidueKey {
    /// Chain identifier.
    pub chain_id: String,
    /// Residue sequence number.
    pub res_seq: i32,
    /// Insertion code.
    pub i_code: String,
}

/// Convert atom records to coarse-grained beads (geometry only).
fn records_to_beads(
    records: &[AtomRecord],
    disulfide_bonds: &[DisulfideBond],
    policy: &dyn CoarseGrain,
) -> Vec<Bead> {
    // Group atoms by residue
    let mut residue_groups: Vec<(ResidueKey, Vec<&AtomRecord>)> = Vec::new();
    let mut key_to_index: HashMap<ResidueKey, usize> = HashMap::new();

    // Separate metal ions
    let mut metal_beads = Vec::new();

    for record in records {
        if residue::is_metal_ion(record) {
            metal_beads.push(Bead {
                x: record.x,
                y: record.y,
                z: record.z,
                charge: metal_charge(&record.res_name),
                mass: residue::atomic_mass(&record.element),
                res_name: record.res_name.clone(),
                chain_id: record.chain_id.clone(),
                res_seq: record.res_seq,
                bead_type: BeadType::Ion,
            });
            continue;
        }

        if !residue::is_amino_acid(&record.res_name) {
            continue;
        }

        let key = ResidueKey {
            chain_id: record.chain_id.clone(),
            res_seq: record.res_seq,
            i_code: record.i_code.clone(),
        };

        if let Some(&idx) = key_to_index.get(&key) {
            residue_groups[idx].1.push(record);
        } else {
            let idx = residue_groups.len();
            key_to_index.insert(key.clone(), idx);
            residue_groups.push((key, vec![record]));
        }
    }

    debug!(
        "{} residues, {} metal ions",
        residue_groups.len(),
        metal_beads.len()
    );

    let ss_bonded = if disulfide_bonds.is_empty() {
        debug!("No _struct_conn disulfides, using geometric detection");
        find_disulfide_bonds_geometric(&residue_groups)
    } else {
        disulfide_bonds_to_keys(disulfide_bonds)
    };
    debug!("{} disulfide-bonded CYS residues", ss_bonded.len());

    let chain_first_last = find_chain_terminals(&residue_groups);

    let mut beads = Vec::new();

    for (key, atoms) in &residue_groups {
        let is_ss_bonded = atoms[0].res_name == "CYS" && ss_bonded.contains(key);

        // Delegate per-residue bead creation to the policy
        beads.extend(policy.residue_to_beads(key, atoms, is_ss_bonded));

        // Terminal beads (shared across all policies)
        if let Some(&(first_idx, last_idx)) = chain_first_last.get(&key.chain_id) {
            let idx = key_to_index[key];
            if idx == first_idx
                && let Some(bead) = make_titratable_bead(key, atoms, &NTERM_GROUP, BeadType::Ntr)
            {
                beads.push(bead);
            }
            if idx == last_idx
                && let Some(bead) = make_titratable_bead(key, atoms, &CTERM_GROUP, BeadType::Ctr)
            {
                beads.push(bead);
            }
        }
    }

    beads.extend(metal_beads);

    let n_titratable = beads.iter().filter(|b| b.bead_type.is_titratable()).count();

    info!(
        "{} residues, {} titratable sites",
        residue_groups.len(),
        n_titratable,
    );

    beads
}

/// Compute multipole moments from beads and external charges.
fn compute_multipole(beads: &[Bead], charges: &[f64]) -> Multipole {
    let charge: f64 = charges.iter().sum();
    let (mut mx, mut my, mut mz) = (0.0, 0.0, 0.0);
    for (b, &q) in beads.iter().zip(charges) {
        mx += q * b.x;
        my += q * b.y;
        mz += q * b.z;
    }
    let dipole = mz.mul_add(mz, mx.mul_add(mx, my * my)).sqrt();
    Multipole {
        charge,
        charge_sq: charge * charge,
        dipole,
        dipole_sq: dipole * dipole,
    }
}

/// Map (chain_id, res_seq) → residue name from backbone/titratable beads.
pub(crate) fn res_name_map(beads: &[Bead]) -> HashMap<(String, i32), String> {
    beads
        .iter()
        .filter(|b| matches!(b.bead_type, BeadType::Backbone | BeadType::Titratable))
        .map(|b| ((b.chain_id.clone(), b.res_seq), b.res_name.clone()))
        .collect()
}

/// Look up the titratable group for a bead, if any.
pub(crate) fn titratable_group_for_bead(
    bead: &Bead,
    res_names: &HashMap<(String, i32), String>,
) -> Option<&'static TitratableGroup> {
    match bead.bead_type {
        BeadType::Sidechain => res_names
            .get(&(bead.chain_id.clone(), bead.res_seq))
            .and_then(|name| residue::find_titratable_group(name)),
        BeadType::Titratable => residue::find_titratable_group(&bead.res_name),
        BeadType::Ntr => Some(&NTERM_GROUP),
        BeadType::Ctr => Some(&CTERM_GROUP),
        _ => None,
    }
}

/// Compute geometric center (unweighted) and total mass of atoms in a single pass.
/// Uses geometric center rather than mass-weighted center; conventional for CG bead placement.
fn center_and_mass(atoms: &[&AtomRecord]) -> Option<((f64, f64, f64), f64)> {
    if atoms.is_empty() {
        return None;
    }
    let (center, mass) = atoms
        .iter()
        .fold(((0.0, 0.0, 0.0), 0.0), |((sx, sy, sz), m), a| {
            let am = residue::atomic_mass(&a.element);
            ((sx + a.x, sy + a.y, sz + a.z), m + am)
        });
    let n = atoms.len() as f64;
    Some(((center.0 / n, center.1 / n, center.2 / n), mass))
}

/// Create a titratable bead at the group's charge center (charge = 0.0).
fn make_titratable_bead(
    key: &ResidueKey,
    atoms: &[&AtomRecord],
    group: &TitratableGroup,
    bead_type: BeadType,
) -> Option<Bead> {
    let center_atoms: Vec<&AtomRecord> = atoms
        .iter()
        .filter(|a| group.center_atoms.contains(&a.name.as_str()))
        .copied()
        .collect();
    let ((x, y, z), mass) = center_and_mass(&center_atoms)?;
    Some(Bead {
        x,
        y,
        z,
        charge: 0.0,
        mass,
        // Named by element (O, N, S) to distinguish from backbone beads (which use residue name)
        res_name: group.element.to_string(),
        chain_id: key.chain_id.clone(),
        res_seq: key.res_seq,
        bead_type,
    })
}

/// Convert parsed `_struct_conn` disulfide bonds to a set of residue keys.
fn disulfide_bonds_to_keys(bonds: &[DisulfideBond]) -> HashSet<ResidueKey> {
    bonds
        .iter()
        .flat_map(|bond| {
            [
                ResidueKey {
                    chain_id: bond.chain_id_1.clone(),
                    res_seq: bond.res_seq_1,
                    i_code: String::new(),
                },
                ResidueKey {
                    chain_id: bond.chain_id_2.clone(),
                    res_seq: bond.res_seq_2,
                    i_code: String::new(),
                },
            ]
        })
        .collect()
}

/// Maximum SG-SG distance (Å) to consider a disulfide bond.
/// Typical S-S bond is ~2.05 Å; 2.5 Å allows for crystallographic uncertainty.
const DISULFIDE_CUTOFF: f64 = 2.5;

/// Find CYS residues involved in disulfide bonds by SG-SG proximity.
/// Used as fallback when _struct_conn records are absent from the mmCIF file.
fn find_disulfide_bonds_geometric(
    residue_groups: &[(ResidueKey, Vec<&AtomRecord>)],
) -> HashSet<ResidueKey> {
    let cys_sg: Vec<(&ResidueKey, f64, f64, f64)> = residue_groups
        .iter()
        .filter(|(_, atoms)| atoms[0].res_name == "CYS")
        .filter_map(|(key, atoms)| {
            atoms
                .iter()
                .find(|a| a.name == "SG")
                .map(|sg| (key, sg.x, sg.y, sg.z))
        })
        .collect();

    let mut ss_bonded = HashSet::new();
    for i in 0..cys_sg.len() {
        for j in (i + 1)..cys_sg.len() {
            let dx = cys_sg[i].1 - cys_sg[j].1;
            let dy = cys_sg[i].2 - cys_sg[j].2;
            let dz = cys_sg[i].3 - cys_sg[j].3;
            let dist = dz.mul_add(dz, dx.mul_add(dx, dy * dy)).sqrt();
            if dist < DISULFIDE_CUTOFF {
                ss_bonded.insert(cys_sg[i].0.clone());
                ss_bonded.insert(cys_sg[j].0.clone());
            }
        }
    }
    ss_bonded
}

/// Find first and last residue index for each chain.
fn find_chain_terminals(
    residue_groups: &[(ResidueKey, Vec<&AtomRecord>)],
) -> HashMap<String, (usize, usize)> {
    let mut terminals: HashMap<String, (usize, usize)> = HashMap::new();
    for (idx, (key, _)) in residue_groups.iter().enumerate() {
        terminals
            .entry(key.chain_id.clone())
            .and_modify(|(_, last)| *last = idx)
            .or_insert((idx, idx));
    }
    terminals
}

/// Common formal charges for metal ions.
fn metal_charge(res_name: &str) -> f64 {
    match res_name {
        "ZN" | "MG" | "CA" | "MN" | "CU" | "CO" | "NI" | "CD" => 2.0,
        "FE" => 3.0, // Fe3+ more common in proteins
        "NA" | "K" => 1.0,
        _ => 0.0,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const TEST_CIF: &str = r#"
data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.auth_atom_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.type_symbol
_atom_site.label_alt_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N ALA A 1 0.000 0.000 0.000 N . 1
ATOM 2 CA ALA A 1 1.458 0.000 0.000 C . 1
ATOM 3 C ALA A 1 2.400 -1.100 0.000 C . 1
ATOM 4 O ALA A 1 2.200 -2.300 0.000 O . 1
ATOM 5 CB ALA A 1 2.000 1.000 0.000 C . 1
ATOM 6 N GLU A 2 3.600 -0.500 0.000 N . 1
ATOM 7 CA GLU A 2 4.900 -1.100 0.000 C . 1
ATOM 8 C GLU A 2 5.800 -2.200 0.000 C . 1
ATOM 9 O GLU A 2 5.600 -3.400 0.000 O . 1
ATOM 10 OXT GLU A 2 6.800 -1.800 0.500 O . 1
ATOM 11 CB GLU A 2 5.500 -0.500 1.200 C . 1
ATOM 12 CG GLU A 2 6.800 -1.000 1.800 C . 1
ATOM 13 CD GLU A 2 7.200 -0.200 3.000 C . 1
ATOM 14 OE1 GLU A 2 6.500 0.700 3.400 O . 1
ATOM 15 OE2 GLU A 2 8.200 -0.500 3.700 O . 1
HETATM 99 ZN ZN A 100 5.000 5.000 5.000 ZN . 1
"#;

    #[test]
    fn coarse_grain_basic() {
        let beads = coarse_grain(TEST_CIF.as_bytes());

        // ALA backbone + NTR + GLU backbone + CTR + GLU sidechain + ZN ion = 6 beads
        assert_eq!(beads.len(), 6, "beads: {beads:#?}");

        // ALA backbone has zero charge
        let ala = &beads[0];
        assert_eq!(ala.res_name, "ALA");
        assert_eq!(ala.bead_type, BeadType::Backbone);
        assert!(ala.charge.abs() < 1e-10);

        // GLU sidechain bead exists, named by element
        let glu_sc = beads
            .iter()
            .find(|b| b.bead_type == BeadType::Sidechain && b.res_seq == 2)
            .unwrap();
        assert_eq!(glu_sc.res_name, "O");

        // ZN ion has formal charge
        let zn = beads.iter().find(|b| b.res_name == "ZN").unwrap();
        assert_eq!(zn.bead_type, BeadType::Ion);
        assert!((zn.charge - 2.0).abs() < 1e-10);
    }

    #[test]
    fn charge_calc_hh() {
        let beads = coarse_grain(TEST_CIF.as_bytes());
        let result = ChargeCalc::new().run(&beads);
        let charged = result.apply(&beads);

        // GLU sidechain should have negative charge at pH 7
        let glu_sc = charged
            .iter()
            .find(|b| b.bead_type == BeadType::Sidechain && b.res_seq == 2)
            .unwrap();
        assert!(
            glu_sc.charge < -0.9,
            "GLU sidechain charge at pH 7: {}",
            glu_sc.charge
        );

        // NTR bead positive at pH 7
        let ntr = charged
            .iter()
            .find(|b| b.bead_type == BeadType::Ntr)
            .unwrap();
        assert!(ntr.charge > 0.9, "NTR charge at pH 7: {}", ntr.charge);

        // CTR bead negative at pH 7
        let ctr = charged
            .iter()
            .find(|b| b.bead_type == BeadType::Ctr)
            .unwrap();
        assert!(ctr.charge < -0.9, "CTR charge at pH 7: {}", ctr.charge);

        // Backbone beads should have zero charge
        assert!(
            charged
                .iter()
                .filter(|b| b.bead_type == BeadType::Backbone)
                .all(|b| b.charge.abs() < 1e-10),
            "backbone beads should have zero charge"
        );
    }

    #[test]
    fn terminal_beads_created() {
        let beads = coarse_grain(TEST_CIF.as_bytes());

        let ntr = beads.iter().find(|b| b.bead_type == BeadType::Ntr).unwrap();
        assert_eq!(ntr.res_seq, 1);

        let ctr = beads.iter().find(|b| b.bead_type == BeadType::Ctr).unwrap();
        assert_eq!(ctr.res_seq, 2);
    }

    #[test]
    fn sidechain_bead_position() {
        let beads = coarse_grain(TEST_CIF.as_bytes());

        let glu_sc = beads
            .iter()
            .find(|b| b.bead_type == BeadType::Sidechain && b.res_seq == 2)
            .unwrap();
        assert!((glu_sc.x - 7.35).abs() < 0.01, "GLU SC x: {}", glu_sc.x);
        assert!((glu_sc.y - 0.1).abs() < 0.01, "GLU SC y: {}", glu_sc.y);
        assert!((glu_sc.z - 3.55).abs() < 0.01, "GLU SC z: {}", glu_sc.z);
    }

    #[test]
    fn water_and_nonprotein_removed() {
        let cif = r#"
data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.auth_atom_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.type_symbol
_atom_site.label_alt_id
_atom_site.pdbx_PDB_model_num
ATOM 1 CA ALA A 1 1.0 2.0 3.0 C . 1
HETATM 2 O HOH B 1 4.0 5.0 6.0 O . 1
HETATM 3 C1 LIG C 1 7.0 8.0 9.0 C . 1
"#;
        let beads = coarse_grain(cif.as_bytes());
        assert_eq!(beads.len(), 1);
        assert_eq!(beads[0].res_name, "ALA");
    }

    #[test]
    fn struct_conn_disulfide_suppresses_titration() {
        let cif = r#"
data_test
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_seq_id
disulf1 disulf A 5 A 20
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.auth_atom_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.type_symbol
_atom_site.label_alt_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N CYS A 5 0.000 0.000 0.000 N . 1
ATOM 2 CA CYS A 5 1.458 0.000 0.000 C . 1
ATOM 3 SG CYS A 5 3.500 1.500 0.000 S . 1
ATOM 4 N CYS A 20 50.000 50.000 50.000 N . 1
ATOM 5 CA CYS A 20 51.458 50.000 50.000 C . 1
ATOM 6 SG CYS A 20 53.500 51.500 50.000 S . 1
"#;
        let beads = coarse_grain(cif.as_bytes());
        assert_eq!(beads.len(), 3, "beads: {beads:#?}");
        assert!(!beads.iter().any(|b| b.bead_type == BeadType::Sidechain));
    }

    #[test]
    fn disulfide_bonded_cys_not_titrated() {
        let cif = r#"
data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.auth_atom_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.type_symbol
_atom_site.label_alt_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N CYS A 5 0.000 0.000 0.000 N . 1
ATOM 2 CA CYS A 5 1.458 0.000 0.000 C . 1
ATOM 3 CB CYS A 5 2.000 1.200 0.000 C . 1
ATOM 4 SG CYS A 5 3.500 1.500 0.000 S . 1
ATOM 5 N CYS A 20 5.000 0.000 0.000 N . 1
ATOM 6 CA CYS A 20 6.458 0.000 0.000 C . 1
ATOM 7 CB CYS A 20 7.000 1.200 0.000 C . 1
ATOM 8 SG CYS A 20 5.500 1.500 0.000 S . 1
"#;
        let beads = coarse_grain(cif.as_bytes());
        assert_eq!(beads.len(), 3, "beads: {beads:#?}");
        assert!(!beads.iter().any(|b| b.bead_type == BeadType::Sidechain));
    }

    #[test]
    fn free_cys_titrates() {
        let cif = r#"
data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.auth_atom_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.type_symbol
_atom_site.label_alt_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N CYS A 5 0.000 0.000 0.000 N . 1
ATOM 2 CA CYS A 5 1.458 0.000 0.000 C . 1
ATOM 3 CB CYS A 5 2.000 1.200 0.000 C . 1
ATOM 4 SG CYS A 5 3.500 1.500 0.000 S . 1
"#;
        let beads = coarse_grain(cif.as_bytes());
        assert_eq!(beads.len(), 3, "beads: {beads:#?}");
        assert!(beads.iter().any(|b| b.bead_type == BeadType::Sidechain));
    }

    #[test]
    fn terminal_with_titratable_sidechain() {
        let cif = r#"
data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.auth_atom_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.type_symbol
_atom_site.label_alt_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N LYS A 1 0.000 0.000 0.000 N . 1
ATOM 2 CA LYS A 1 1.458 0.000 0.000 C . 1
ATOM 3 C LYS A 1 2.400 -1.100 0.000 C . 1
ATOM 4 O LYS A 1 2.200 -2.300 0.000 O . 1
ATOM 5 CB LYS A 1 2.000 1.200 0.000 C . 1
ATOM 6 CG LYS A 1 3.500 1.500 0.000 C . 1
ATOM 7 CD LYS A 1 5.000 1.200 0.000 C . 1
ATOM 8 CE LYS A 1 6.500 1.500 0.000 C . 1
ATOM 9 NZ LYS A 1 8.000 1.200 0.000 N . 1
ATOM 10 N ALA A 2 3.600 -0.500 0.000 N . 1
ATOM 11 CA ALA A 2 4.900 -1.100 0.000 C . 1
ATOM 12 C ALA A 2 5.800 -2.200 0.000 C . 1
ATOM 13 O ALA A 2 5.600 -3.400 0.000 O . 1
ATOM 14 OXT ALA A 2 6.800 -1.800 0.500 O . 1
"#;
        let beads = coarse_grain(cif.as_bytes());

        // LYS backbone + NTR(N) + sidechain(N at NZ) + ALA backbone + CTR(O) = 5
        assert_eq!(beads.len(), 5, "beads: {beads:#?}");

        let lys_beads: Vec<_> = beads.iter().filter(|b| b.res_seq == 1).collect();
        assert_eq!(
            lys_beads.len(),
            3,
            "LYS should have 3 beads: {lys_beads:#?}"
        );
        assert!(lys_beads.iter().any(|b| b.bead_type == BeadType::Backbone));
        assert!(lys_beads.iter().any(|b| b.bead_type == BeadType::Ntr));
        assert!(lys_beads.iter().any(|b| b.bead_type == BeadType::Sidechain));

        // NTR bead at N atom position
        let ntr = lys_beads
            .iter()
            .find(|b| b.bead_type == BeadType::Ntr)
            .unwrap();
        assert_eq!(ntr.res_name, "N");
        assert!((ntr.x - 0.0).abs() < 0.01);

        // Sidechain bead at NZ position
        let sc = lys_beads
            .iter()
            .find(|b| b.bead_type == BeadType::Sidechain)
            .unwrap();
        assert_eq!(sc.res_name, "N");
        assert!((sc.x - 8.0).abs() < 0.01);
    }

    #[test]
    fn multipole_moments() {
        let beads = coarse_grain(TEST_CIF.as_bytes());
        let result = ChargeCalc::new().run(&beads);

        // HH is deterministic: charge_sq == charge²
        let m = &result.multipole;
        assert!(
            (m.charge_sq - m.charge * m.charge).abs() < 1e-10,
            "HH: ⟨Z²⟩ should equal ⟨Z⟩²"
        );
        assert!(
            (m.dipole_sq - m.dipole * m.dipole).abs() < 1e-10,
            "HH: ⟨μ²⟩ should equal ⟨μ⟩²"
        );
    }
}
