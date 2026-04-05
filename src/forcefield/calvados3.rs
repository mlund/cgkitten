//! Calvados 3 coarse-grained amino acid model.
//!
//! Parameters from <https://doi.org/10.1093/database/baz024>.

use super::{BeadParams, ForceField, HydrophobicScaling};
use crate::BeadType;
use crate::residue::HYDROPHOBIC_RESIDUES;

/// Calvados 3 force field with optional hydrophobic pair scaling.
pub struct Calvados3 {
    scaling: HydrophobicScaling,
}

impl Calvados3 {
    /// Create with the given hydrophobic pair scaling policy.
    pub fn new(scaling: HydrophobicScaling) -> Self {
        Self { scaling }
    }
}

// Masses are free amino acid mass minus one water (18.02 Da) per peptide bond.
const RESIDUES: &[(&str, f64, f64, f64, f64)] = &[
    ("ALA", 71.07, 0.3377244362031627, 5.04, 0.8368),
    ("ARG", 156.19, 0.7407902764839954, 6.56, 0.8368),
    ("ASN", 114.10, 0.3706962163690402, 5.68, 0.8368),
    ("ASP", 115.09, 0.092587557536158, 5.58, 0.8368),
    ("CYS", 103.14, 0.5922529084601322, 5.48, 0.8368),
    ("GLN", 128.13, 0.3143449791669133, 6.02, 0.8368),
    ("GLU", 129.11, 0.000249590539426, 5.92, 0.8368),
    ("GLY", 57.05, 0.7538308115197386, 4.50, 0.8368),
    ("HIS", 137.14, 0.4087176216525476, 6.08, 0.8368),
    ("ILE", 113.16, 0.5130398874425708, 6.18, 0.8368),
    ("LEU", 113.16, 0.5548615312993875, 6.18, 0.8368),
    ("LYS", 128.17, 0.1380602542039267, 6.36, 0.8368),
    ("MET", 131.20, 0.5170874160398543, 6.18, 0.8368),
    ("PHE", 147.18, 0.8906449355499866, 6.36, 0.8368),
    ("PRO", 97.12, 0.3469777523519372, 5.56, 0.8368),
    ("SER", 87.08, 0.4473142572693176, 5.18, 0.8368),
    ("THR", 101.11, 0.2672387936544146, 5.62, 0.8368),
    ("TRP", 186.22, 1.033450123574512, 6.78, 0.8368),
    ("TYR", 163.18, 0.950628687301107, 6.46, 0.8368),
    ("VAL", 99.13, 0.2936174211771383, 5.86, 0.8368),
];

/// Virtual/terminal beads: point charges with no hydrophobic interaction.
/// λ=0 so AH reduces to pure repulsion; mass=0 because the residue bead carries it.
const SITE: BeadParams = BeadParams {
    mass: 0.0,
    sigma: 2.0,
    epsilon: 0.8368,
    lambda: 0.0,
};

impl ForceField for Calvados3 {
    fn params(&self, res_name: &str, bead_type: BeadType) -> Option<BeadParams> {
        match bead_type {
            BeadType::Residue | BeadType::Titratable => RESIDUES
                .iter()
                .find(|(name, _, _, _, _)| *name == res_name)
                .map(|(_, mass, lambda, sigma, epsilon)| BeadParams {
                    mass: *mass,
                    sigma: *sigma,
                    epsilon: *epsilon,
                    lambda: *lambda,
                }),
            BeadType::Ion | BeadType::Virtual | BeadType::Ntr | BeadType::Ctr => Some(SITE),
        }
    }

    fn format_atom_fields(&self, p: &BeadParams) -> String {
        format!(
            "σ: {}, ε: {}, hydrophobicity: !Lambda {}",
            p.sigma, p.epsilon, p.lambda
        )
    }

    fn nonbonded_yaml(&self, types: &[(&str, f64, BeadParams)]) -> String {
        let mut yaml = "\nsystem:\n  energy:\n    nonbonded:\n      \
            # A Coulomb term is automatically added\n      \
            default:\n        - !AshbaughHatch {mixing: arithmetic, cutoff: 20.0}\n"
            .to_string();

        // Hydrophobic residues are all neutral, so `replace:` is safe —
        // no Coulomb interaction to inherit from default.
        let pairs = hydrophobic_pairs(types, &self.scaling);
        if !pairs.is_empty() {
            yaml.push_str("      replace:\n");
            for (na, nb, sigma, epsilon, lambda) in &pairs {
                yaml.push_str(&format!(
                    "        [{na}, {nb}]:\n          \
                     - !AshbaughHatch {{σ: {sigma:.4}, ε: {epsilon:.4}, λ: {lambda:.4}, cutoff: 20.0}}\n"
                ));
            }
        }
        yaml
    }
}

/// Generate scaled hydrophobic pairs using Lorentz-Berthelot mixing.
fn hydrophobic_pairs<'a>(
    types: &[(&'a str, f64, BeadParams)],
    scaling: &HydrophobicScaling,
) -> Vec<(&'a str, &'a str, f64, f64, f64)> {
    if *scaling == HydrophobicScaling::NoScale {
        return Vec::new();
    }
    let hp: Vec<_> = types
        .iter()
        .filter(|(name, _, _)| HYDROPHOBIC_RESIDUES.contains(name))
        .collect();

    let mut pairs = Vec::new();
    for (i, (na, _, pa)) in hp.iter().enumerate() {
        for (nb, _, pb) in &hp[i..] {
            let sigma = (pa.sigma + pb.sigma) / 2.0;
            let mut epsilon = (pa.epsilon * pb.epsilon).sqrt();
            let mut lambda = (pa.lambda + pb.lambda) / 2.0;

            match scaling {
                HydrophobicScaling::ScaleLambda(c) => lambda *= c,
                HydrophobicScaling::ScaleEpsilon(c) => epsilon *= c,
                HydrophobicScaling::NoScale => unreachable!(),
            }
            pairs.push((*na, *nb, sigma, epsilon, lambda));
        }
    }
    pairs
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_types<'a>(names: &[&'a str]) -> Vec<(&'a str, f64, BeadParams)> {
        names
            .iter()
            .filter_map(|name| {
                RESIDUES.iter().find(|(n, _, _, _, _)| n == name).map(
                    |(_, mass, lambda, sigma, epsilon)| {
                        (
                            *name,
                            0.0,
                            BeadParams {
                                mass: *mass,
                                sigma: *sigma,
                                epsilon: *epsilon,
                                lambda: *lambda,
                            },
                        )
                    },
                )
            })
            .collect()
    }

    #[test]
    fn lorentz_berthelot_mixing() {
        let types = make_types(&["ALA", "VAL"]);
        let pairs = hydrophobic_pairs(&types, &HydrophobicScaling::ScaleLambda(1.0));
        assert_eq!(pairs.len(), 3); // ALA-ALA, ALA-VAL, VAL-VAL

        let (_, _, sigma, epsilon, lambda) = &pairs[1]; // ALA-VAL
        assert!((sigma - 5.45).abs() < 1e-10);
        assert!((epsilon - 0.8368).abs() < 1e-4);
        assert!((lambda - 0.31565).abs() < 1e-4);
    }

    #[test]
    fn scale_lambda() {
        let types = make_types(&["PHE"]);
        let pairs = hydrophobic_pairs(&types, &HydrophobicScaling::ScaleLambda(1.2));
        assert_eq!(pairs.len(), 1);
        assert!((pairs[0].4 - 0.8906449355499866 * 1.2).abs() < 1e-10);
    }

    #[test]
    fn scale_epsilon() {
        let types = make_types(&["TRP"]);
        let pairs = hydrophobic_pairs(&types, &HydrophobicScaling::ScaleEpsilon(0.8));
        assert_eq!(pairs.len(), 1);
        assert!((pairs[0].3 - 0.8368 * 0.8).abs() < 1e-10);
    }

    #[test]
    fn filters_non_hydrophobic() {
        let types = make_types(&["ALA", "GLU"]);
        let pairs = hydrophobic_pairs(&types, &HydrophobicScaling::ScaleLambda(1.0));
        assert_eq!(pairs.len(), 1); // only ALA-ALA
    }
}
