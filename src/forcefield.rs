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
            BeadType::Backbone | BeadType::Titratable => CALVADOS3_BACKBONE
                .iter()
                .find(|(name, _, _, _)| *name == res_name)
                .map(|(_, lambda, sigma, epsilon)| BeadParams {
                    sigma: *sigma,
                    epsilon: *epsilon,
                    lambda: *lambda,
                }),
            // Ions get a minimal LJ entry so downstream tools (e.g. duello) have σ/ε defined
            BeadType::Ion => Some(CALVADOS3_SITE),
            BeadType::Sidechain | BeadType::Ntr | BeadType::Ctr => Some(CALVADOS3_SITE),
        }
    }

    fn system_yaml(&self) -> &'static str {
        "\nsystem:\n  energy:\n    nonbonded:\n      # A Coulomb term is automatically added\n      default:\n        - !AshbaughHatch {mixing: arithmetic, cutoff: 20.0}\n"
    }
}

/// Hydrophobic scaling policy for nonbonded pair overrides.
#[derive(Clone, Debug, Default)]
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
        Err(format!(
            "expected 'lambda:<f>' or 'epsilon:<f>'; got '{s}'"
        ))
    }
}

/// A pairwise nonbonded interaction override.
pub struct PairInteraction {
    /// First bead type name.
    pub name_a: String,
    /// Second bead type name.
    pub name_b: String,
    /// Lorentz-Berthelot mixed σ (Å).
    pub sigma: f64,
    /// Lorentz-Berthelot mixed ε (kJ/mol), possibly scaled.
    pub epsilon: f64,
    /// Mixed λ (Ashbaugh-Hatch), possibly scaled.
    pub lambda: f64,
}

/// Generate pair interaction overrides for hydrophobic residue pairs.
///
/// Filters `types` to those whose name appears in `hydrophobic`, computes
/// Lorentz-Berthelot mixing rules for all unique pairs (including self-pairs),
/// and applies the scaling policy.
pub fn hydrophobic_pairs(
    types: &[(String, BeadParams)],
    hydrophobic: &[&str],
    scaling: &HydrophobicScaling,
) -> Vec<PairInteraction> {
    let hp: Vec<_> = types
        .iter()
        .filter(|(name, _)| hydrophobic.contains(&name.as_str()))
        .collect();

    // Triangular loop (i <= j) to emit each unique pair once.
    // Self-pairs (i == j) are included because Faunus pair entries
    // override the default for those types entirely.
    let mut pairs = Vec::new();
    for i in 0..hp.len() {
        for j in i..hp.len() {
            let (na, pa) = &hp[i];
            let (nb, pb) = &hp[j];
            // Lorentz-Berthelot combining rules:
            // σ_mix = arithmetic mean, ε_mix = geometric mean, λ_mix = arithmetic mean
            let sigma = (pa.sigma + pb.sigma) / 2.0;
            let mut epsilon = (pa.epsilon * pb.epsilon).sqrt();
            let mut lambda = (pa.lambda + pb.lambda) / 2.0;

            match scaling {
                HydrophobicScaling::NoScale => {}
                HydrophobicScaling::ScaleLambda(c) => lambda *= c,
                HydrophobicScaling::ScaleEpsilon(c) => epsilon *= c,
            }

            pairs.push(PairInteraction {
                name_a: na.clone(),
                name_b: nb.clone(),
                sigma,
                epsilon,
                lambda,
            });
        }
    }
    pairs
}

/// Format a pair interaction as a YAML entry for the Faunus `nonbonded:` section.
/// Indented to sit as a sibling of `default:` inside `system: energy: nonbonded:`.
pub fn format_pair_yaml(pair: &PairInteraction) -> String {
    format!(
        "      [{}, {}]:\n        - !AshbaughHatch {{σ: {:.4}, ε: {:.4}, λ: {:.4}}}\n",
        pair.name_a, pair.name_b, pair.sigma, pair.epsilon, pair.lambda,
    )
}

/// Create a force field model by name.
pub fn from_name(name: &str) -> Option<Box<dyn ForceField>> {
    match name {
        "calvados3" => Some(Box::new(Calvados3)),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lorentz_berthelot_mixing() {
        let types = vec![
            ("ALA".into(), BeadParams { sigma: 5.04, epsilon: 0.8368, lambda: 0.3377 }),
            ("VAL".into(), BeadParams { sigma: 5.86, epsilon: 0.8368, lambda: 0.2936 }),
        ];
        let pairs = hydrophobic_pairs(&types, &["ALA", "VAL"], &HydrophobicScaling::NoScale);
        assert_eq!(pairs.len(), 3); // ALA-ALA, ALA-VAL, VAL-VAL

        // ALA-VAL cross pair
        let cross = &pairs[1];
        assert_eq!(cross.name_a, "ALA");
        assert_eq!(cross.name_b, "VAL");
        assert!((cross.sigma - 5.45).abs() < 1e-10);
        assert!((cross.epsilon - 0.8368).abs() < 1e-4);
        assert!((cross.lambda - 0.31565).abs() < 1e-4);
    }

    #[test]
    fn scale_lambda() {
        let types = vec![
            ("PHE".into(), BeadParams { sigma: 6.36, epsilon: 0.8368, lambda: 0.8906 }),
        ];
        let pairs = hydrophobic_pairs(&types, &["PHE"], &HydrophobicScaling::ScaleLambda(1.2));
        assert_eq!(pairs.len(), 1);
        assert!((pairs[0].lambda - 0.8906 * 1.2).abs() < 1e-10);
        assert!((pairs[0].epsilon - 0.8368).abs() < 1e-10);
    }

    #[test]
    fn scale_epsilon() {
        let types = vec![
            ("TRP".into(), BeadParams { sigma: 6.78, epsilon: 0.8368, lambda: 1.0335 }),
        ];
        let pairs = hydrophobic_pairs(&types, &["TRP"], &HydrophobicScaling::ScaleEpsilon(0.8));
        assert_eq!(pairs.len(), 1);
        assert!((pairs[0].epsilon - 0.8368 * 0.8).abs() < 1e-10);
        assert!((pairs[0].lambda - 1.0335).abs() < 1e-10);
    }

    #[test]
    fn filters_non_hydrophobic() {
        let types = vec![
            ("ALA".into(), BeadParams { sigma: 5.04, epsilon: 0.8368, lambda: 0.3377 }),
            ("GLU".into(), BeadParams { sigma: 5.92, epsilon: 0.8368, lambda: 0.0002 }),
        ];
        let pairs = hydrophobic_pairs(&types, &["ALA"], &HydrophobicScaling::NoScale);
        assert_eq!(pairs.len(), 1); // only ALA-ALA
        assert_eq!(pairs[0].name_a, "ALA");
        assert_eq!(pairs[0].name_b, "ALA");
    }
}
