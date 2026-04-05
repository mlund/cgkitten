//! Pasquier LJ model with uniform hydrophobic interaction.
//!
//! - Pasquier et al. (2023), <https://doi.org/10.1016/j.jcis.2022.08.054>
//! - Sigma values from Kim & Hummer, <https://doi.org/10.1016/j.jmb.2007.11.063>
//!
//! LJ with two ε classes (Table 3, last row: open conformation):
//! hydrophobic–hydrophobic ε_hh = 0.485 kT, all other ε = 0.005 kT.
//!
//! Note: σ values are taken from Kim & Hummer rather than the original
//! Pasquier model.

use super::{BeadParams, ForceField, HydrophobicScaling};
use crate::BeadType;
use crate::residue::HYDROPHOBIC_RESIDUES;

const R_KJ_PER_MOL: f64 = 8.314462618e-3;
const T_REF: f64 = 298.15;
const KT_TO_KJ_MOL: f64 = R_KJ_PER_MOL * T_REF;

/// ε for hydrophobic–hydrophobic pairs (kT at 298.15 K).
const EPSILON_HH_KT: f64 = 0.485;
/// ε for all other pairs (kT at 298.15 K).
const EPSILON_OTHER_KT: f64 = 0.005;

/// Pasquier LJ force field with optional hydrophobic epsilon scaling.
pub struct Pasquier {
    scaling: HydrophobicScaling,
}

impl Pasquier {
    /// Create with the given hydrophobic epsilon scaling policy.
    pub fn new(scaling: HydrophobicScaling) -> Self {
        Self { scaling }
    }
}

// Reuse KH sigma and mass values: (name, mass_da, sigma_ii)
const RESIDUES: &[(&str, f64, f64)] = &[
    ("ALA", 71.07, 5.0),
    ("ARG", 156.19, 6.7),
    ("ASN", 114.10, 5.7),
    ("ASP", 115.09, 5.6),
    ("CYS", 103.14, 5.5),
    ("GLN", 128.13, 6.0),
    ("GLU", 129.11, 5.9),
    ("GLY", 57.05, 4.5),
    ("HIS", 137.14, 6.1),
    ("ILE", 113.16, 6.2),
    ("LEU", 113.16, 6.2),
    ("LYS", 128.17, 6.4),
    ("MET", 131.20, 6.2),
    ("PHE", 147.18, 6.4),
    ("PRO", 97.12, 5.6),
    ("SER", 87.08, 5.2),
    ("THR", 101.11, 5.6),
    ("TRP", 186.22, 6.8),
    ("TYR", 163.18, 6.5),
    ("VAL", 99.13, 5.9),
];

/// Virtual/terminal site beads — weak LJ to avoid dominating at short range.
const SITE: BeadParams = BeadParams {
    mass: 0.0,
    sigma: 2.0,
    epsilon: EPSILON_OTHER_KT * KT_TO_KJ_MOL,
    lambda: 0.0,
};

impl ForceField for Pasquier {
    fn params(&self, res_name: &str, bead_type: BeadType) -> Option<BeadParams> {
        match bead_type {
            BeadType::Residue | BeadType::Titratable => RESIDUES
                .iter()
                .find(|(name, _, _)| *name == res_name)
                .map(|(_, mass, sigma)| BeadParams {
                    mass: *mass,
                    sigma: *sigma,
                    epsilon: 0.0,
                    lambda: 0.0,
                }),
            BeadType::Ion | BeadType::Virtual | BeadType::Ntr | BeadType::Ctr => Some(SITE),
        }
    }

    fn format_atom_fields(&self, p: &BeadParams) -> String {
        format!("σ: {}", p.sigma)
    }

    fn nonbonded_yaml(&self, types: &[(&str, f64, BeadParams)]) -> String {
        let mut yaml = "\nsystem:\n  energy:\n    nonbonded:\n".to_string();

        // Neutral pairs → replace (skip redundant Coulomb), charged → append (inherit Coulomb)
        let mut replace_pairs = Vec::new();
        let mut append_pairs = Vec::new();

        for (i, (na, qa, pa)) in types.iter().enumerate() {
            for (nb, qb, pb) in &types[i..] {
                let sigma = (pa.sigma + pb.sigma) / 2.0;
                let is_hh = is_hydrophobic(na) && is_hydrophobic(nb);

                let mut epsilon_kt = if is_hh {
                    EPSILON_HH_KT
                } else {
                    EPSILON_OTHER_KT
                };

                if let HydrophobicScaling::ScaleEpsilon(c) = &self.scaling
                    && is_hh
                {
                    epsilon_kt *= c;
                }

                let epsilon = epsilon_kt * KT_TO_KJ_MOL;
                let entry = format!(
                    "        [{na}, {nb}]:\n          \
                     - !LennardJones {{sigma: {sigma}, epsilon: {epsilon:.4}}}\n"
                );

                if qa.abs() < f64::EPSILON && qb.abs() < f64::EPSILON {
                    replace_pairs.push(entry);
                } else {
                    append_pairs.push(entry);
                }
            }
        }

        if !replace_pairs.is_empty() {
            yaml.push_str("      replace:\n");
            for entry in &replace_pairs {
                yaml.push_str(entry);
            }
        }
        if !append_pairs.is_empty() {
            yaml.push_str("      append:\n");
            for entry in &append_pairs {
                yaml.push_str(entry);
            }
        }
        yaml
    }
}

fn is_hydrophobic(name: &str) -> bool {
    HYDROPHOBIC_RESIDUES.contains(&name)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residue_lookup() {
        let ff = Pasquier::new(HydrophobicScaling::NoScale);
        let p = ff.params("ALA", BeadType::Residue).unwrap();
        assert!((p.sigma - 5.0).abs() < 1e-10);
        assert!((p.mass - 71.07).abs() < 1e-10);
    }

    #[test]
    fn hydrophobic_pair_epsilon() {
        let ff = Pasquier::new(HydrophobicScaling::NoScale);
        let types = vec![
            ("ALA", 0.0, ff.params("ALA", BeadType::Residue).unwrap()),
            ("VAL", 0.0, ff.params("VAL", BeadType::Residue).unwrap()),
        ];
        let yaml = ff.nonbonded_yaml(&types);
        let expected = EPSILON_HH_KT * KT_TO_KJ_MOL;
        assert!(yaml.contains(&format!("{expected:.4}")));
    }

    #[test]
    fn non_hydrophobic_pair_epsilon() {
        let ff = Pasquier::new(HydrophobicScaling::NoScale);
        let types = vec![
            ("ALA", 0.0, ff.params("ALA", BeadType::Residue).unwrap()),
            ("GLU", -1.0, ff.params("GLU", BeadType::Residue).unwrap()),
        ];
        let yaml = ff.nonbonded_yaml(&types);
        let expected = EPSILON_OTHER_KT * KT_TO_KJ_MOL;
        // ALA-GLU is not hydrophobic-hydrophobic
        assert!(yaml.contains(&format!("{expected:.4}")));
    }

    #[test]
    fn epsilon_scaling() {
        let ff = Pasquier::new(HydrophobicScaling::ScaleEpsilon(0.5));
        let types = vec![("ALA", 0.0, ff.params("ALA", BeadType::Residue).unwrap())];
        let yaml = ff.nonbonded_yaml(&types);
        let expected = EPSILON_HH_KT * 0.5 * KT_TO_KJ_MOL;
        assert!(yaml.contains(&format!("{expected:.4}")));
    }

    #[test]
    fn charged_pairs_in_append() {
        let ff = Pasquier::new(HydrophobicScaling::NoScale);
        let types = vec![
            ("ALA", 0.0, ff.params("ALA", BeadType::Residue).unwrap()),
            ("ARG", 1.0, ff.params("ARG", BeadType::Residue).unwrap()),
        ];
        let yaml = ff.nonbonded_yaml(&types);
        assert!(yaml.contains("replace:"));
        assert!(yaml.contains("append:"));
    }
}
