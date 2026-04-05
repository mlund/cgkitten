//! Force field models for coarse-grained beads.
//!
//! Each model implements [`ForceField`] to provide per-residue parameters and
//! generate the nonbonded YAML section for Faunus topology files.

mod calvados3;
mod kimhummer;
mod pasquier;

use crate::BeadType;
pub use calvados3::Calvados3;
pub use kimhummer::KimHummer;
pub use pasquier::Pasquier;

/// Per-bead force field parameters.
#[derive(Clone, Copy)]
pub struct BeadParams {
    /// Residue mass (Da). 0.0 for site beads whose mass is carried by the backbone bead.
    pub mass: f64,
    /// Lennard-Jones diameter σ (Å).
    pub sigma: f64,
    /// Lennard-Jones well depth ε (kJ/mol).
    pub epsilon: f64,
    /// Hydrophobicity parameter λ (Ashbaugh-Hatch). Zero for non-AH models.
    pub lambda: f64,
}

/// A force field model provides parameters for coarse-grained beads
/// and generates the nonbonded YAML section.
pub trait ForceField {
    /// Look up parameters for a bead by residue name and type.
    fn params(&self, res_name: &str, bead_type: BeadType) -> Option<BeadParams>;

    /// FF-specific per-atom fields for YAML output (appended after charge/mass/name).
    fn format_atom_fields(&self, params: &BeadParams) -> String;

    /// Complete `\nsystem:\n  energy:\n    nonbonded:\n...` YAML block.
    /// `types` provides (name, charge, params) for each atom type in the system.
    fn nonbonded_yaml(&self, types: &[(&str, f64, BeadParams)]) -> String;
}

/// Hydrophobic scaling policy for nonbonded pair overrides.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub enum HydrophobicScaling {
    /// No scaling applied.
    #[default]
    NoScale,
    /// Scale the λ (hydrophobicity) parameter by the given factor.
    ScaleLambda(f64),
    /// Scale the ε (well depth) parameter by the given factor.
    ScaleEpsilon(f64),
}

impl std::str::FromStr for HydrophobicScaling {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some(rest) = s.strip_prefix("lambda:") {
            return rest
                .parse::<f64>()
                .map(Self::ScaleLambda)
                .map_err(|e| format!("bad lambda value: {e}"));
        }
        if let Some(rest) = s.strip_prefix("epsilon:") {
            return rest
                .parse::<f64>()
                .map(Self::ScaleEpsilon)
                .map_err(|e| format!("bad epsilon value: {e}"));
        }
        Err(format!("expected 'lambda:<f>' or 'epsilon:<f>'; got '{s}'"))
    }
}

impl std::fmt::Display for HydrophobicScaling {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HydrophobicScaling::NoScale => f.write_str("none"),
            HydrophobicScaling::ScaleLambda(c) => write!(f, "lambda×{c}"),
            HydrophobicScaling::ScaleEpsilon(c) => write!(f, "epsilon×{c}"),
        }
    }
}

/// Create a force field model by name.
///
/// `scaling` is passed to models that support hydrophobic pair scaling.
/// Returns `Err` if the scaling is incompatible with the model.
/// Returns `Ok(None)` for unknown names (use `"none"` to explicitly opt out).
pub fn from_name(
    name: &str,
    scaling: HydrophobicScaling,
) -> Result<Option<Box<dyn ForceField>>, String> {
    match name {
        "calvados3" => Ok(Some(Box::new(Calvados3::new(scaling)))),
        "kimhummer" | "kh" => {
            if matches!(scaling, HydrophobicScaling::ScaleLambda(_)) {
                return Err(
                    "lambda scaling is not supported for Kim-Hummer (use epsilon:<factor>)"
                        .to_string(),
                );
            }
            Ok(Some(Box::new(KimHummer::new(scaling))))
        }
        "pasquier" => {
            if matches!(scaling, HydrophobicScaling::ScaleLambda(_)) {
                return Err(
                    "lambda scaling is not supported for Pasquier (use epsilon:<factor>)"
                        .to_string(),
                );
            }
            Ok(Some(Box::new(Pasquier::new(scaling))))
        }
        _ => Ok(None),
    }
}
