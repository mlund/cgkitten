//! Force field model definitions for coarse-grained beads.
//!
//! Each model provides per-residue parameters (σ, ε, hydrophobicity λ) and
//! a system YAML block for energy definitions.

use crate::BeadType;

/// Per-bead force field parameters.
#[derive(Clone, Copy)]
pub struct BeadParams {
    /// Lennard-Jones diameter σ (Å).
    pub sigma: f64,
    /// Lennard-Jones well depth ε (kJ/mol).
    pub epsilon: f64,
    /// Hydrophobicity parameter λ (Ashbaugh-Hatch).
    pub lambda: f64,
}

/// A force field model provides parameters for coarse-grained beads.
pub trait ForceField {
    /// Look up parameters for a bead by residue name and type.
    fn params(&self, res_name: &str, bead_type: BeadType) -> Option<BeadParams>;

    /// YAML `system:` block appended to topology output.
    fn system_yaml(&self) -> &'static str;
}

/// Calvados 3 coarse-grained amino acid model.
///
/// Parameters from <https://doi.org/10.1093/database/baz024>.
pub struct Calvados3;

// (res_name, lambda, sigma, epsilon)
const CALVADOS3_BACKBONE: &[(&str, f64, f64, f64)] = &[
    ("ALA", 0.3377244362031627, 5.04, 0.8368),
    ("ARG", 0.7407902764839954, 6.56, 0.8368),
    ("ASN", 0.3706962163690402, 5.68, 0.8368),
    ("ASP", 0.092587557536158, 5.58, 0.8368),
    ("CYS", 0.5922529084601322, 5.48, 0.8368),
    ("GLN", 0.3143449791669133, 6.02, 0.8368),
    ("GLU", 0.000249590539426, 5.92, 0.8368),
    ("GLY", 0.7538308115197386, 4.50, 0.8368),
    ("HIS", 0.4087176216525476, 6.08, 0.8368),
    ("ILE", 0.5130398874425708, 6.18, 0.8368),
    ("LEU", 0.5548615312993875, 6.18, 0.8368),
    ("LYS", 0.1380602542039267, 6.36, 0.8368),
    ("MET", 0.5170874160398543, 6.18, 0.8368),
    ("PHE", 0.8906449355499866, 6.36, 0.8368),
    ("PRO", 0.3469777523519372, 5.56, 0.8368),
    ("SER", 0.4473142572693176, 5.18, 0.8368),
    ("THR", 0.2672387936544146, 5.62, 0.8368),
    ("TRP", 1.033450123574512, 6.78, 0.8368),
    ("TYR", 0.950628687301107, 6.46, 0.8368),
    ("VAL", 0.2936174211771383, 5.86, 0.8368),
];

/// Sidechain/terminal beads all share the same σ, ε, λ in Calvados 3.
/// λ = 0 because titratable site beads are point charges with no hydrophobic interaction.
const CALVADOS3_SITE: BeadParams = BeadParams {
    sigma: 2.0,
    epsilon: 0.8368,
    lambda: 0.0,
};

impl ForceField for Calvados3 {
    fn params(&self, res_name: &str, bead_type: BeadType) -> Option<BeadParams> {
        match bead_type {
            BeadType::Backbone => CALVADOS3_BACKBONE
                .iter()
                .find(|(name, _, _, _)| *name == res_name)
                .map(|(_, lambda, sigma, epsilon)| BeadParams {
                    sigma: *sigma,
                    epsilon: *epsilon,
                    lambda: *lambda,
                }),
            BeadType::Ion => None, // ions carry only Coulomb charge, no LJ/AH params
            BeadType::Sidechain | BeadType::Ntr | BeadType::Ctr => Some(CALVADOS3_SITE),
        }
    }

    fn system_yaml(&self) -> &'static str {
        "\nsystem:\n  energy:\n    nonbonded:\n      # A Coulomb term is automatically added\n      default:\n        - !AshbaughHatch {mixing: arithmetic, cutoff: 20.0}\n"
    }
}

/// Create a force field model by name.
pub fn from_name(name: &str) -> Option<Box<dyn ForceField>> {
    match name {
        "calvados3" => Some(Box::new(Calvados3)),
        _ => None,
    }
}
