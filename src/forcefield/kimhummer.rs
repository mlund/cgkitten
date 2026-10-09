//! Kim-Hummer coarse-grained amino acid model.
//!
//! - Kim & Hummer (2008), <https://doi.org/10.1016/j.jmb.2007.11.063>
//! - Miyazawa & Jernigan (1996), <https://doi.org/10.1006/jmbi.1996.0114>
//! - Parameters from <https://github.com/bio-phys/complexespp>

use super::{BeadParams, ForceField, HydrophobicScaling};
use crate::BeadType;

const R_KJ_PER_MOL: f64 = 8.314462618e-3; // kJ/(mol·K)
const KT_TO_KJ_MOL: f64 = R_KJ_PER_MOL * 300.0; // T_ref = 300 K

/// Kim-Hummer force field with optional hydrophobic epsilon scaling.
pub struct KimHummer {
    scaling: HydrophobicScaling,
}

impl KimHummer {
    /// Create with the given hydrophobic epsilon scaling policy.
    pub fn new(scaling: HydrophobicScaling) -> Self {
        Self { scaling }
    }

    /// Look up residue index in AA_ORDER, or None for non-residue types.
    fn residue_index(name: &str) -> Option<usize> {
        RESIDUES.iter().position(|(n, _, _)| *n == name)
    }
}

// (name, mass_da, sigma_ii) — self-pair diameter from KH forcefield
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

/// Virtual/terminal site beads — small σ and weak ε to avoid
/// dominating the pair interaction when sites overlap at short range.
const SITE: BeadParams = BeadParams {
    mass: 0.0,
    sigma: 2.0,
    epsilon: 0.1,
    lambda: 0.0,
    residue: None,
};

/// Pre-scaled Miyazawa-Jernigan contact energies: ε_ij = λ(e_ij − e_0)
/// with λ = 0.159, e_0 = −2.27 kT. Stored in kT at T_ref = 300 K;
/// converted to kJ/mol during YAML generation.
/// Rows/columns follow RESIDUES order (ALA, ARG, ..., VAL).
#[rustfmt::skip]
const EPSILON_KT: [[f64; 20]; 20] = [
    [-0.07155, 0.06996, 0.06837, 0.09063, -0.2067, 0.06042, 0.12084, -0.00636, -0.02226, -0.36729, -0.41976, 0.15264, -0.26553, -0.40386, 0.03816, 0.04134, -0.00795, -0.24645, -0.17331, -0.28143],
    [0.06996, 0.11448, 0.10017, -0.00318, -0.0477, 0.07473, 0.0, 0.08745, 0.01749, -0.21624, -0.27984, 0.26712, -0.13515, -0.27189, 0.09063, 0.10335, 0.05883, -0.18126, -0.14151, -0.1272],
    [0.06837, 0.10017, 0.09381, 0.09381, -0.05088, 0.08904, 0.12084, 0.08427, 0.03021, -0.15423, -0.23373, 0.16854, -0.10812, -0.23532, 0.11766, 0.10971, 0.06201, -0.1272, -0.07791, -0.08904],
    [0.09063, -0.00318, 0.09381, 0.16854, -0.02226, 0.12879, 0.19875, 0.10812, -0.00795, -0.1431, -0.17967, 0.09381, -0.0477, -0.19239, 0.14946, 0.10176, 0.07473, -0.09063, -0.07791, -0.03339],
    [-0.2067, -0.0477, -0.05088, -0.02226, -0.50403, -0.09222, 0.0, -0.14151, -0.21147, -0.51357, -0.56604, 0.05088, -0.43248, -0.56127, -0.1272, -0.09381, -0.13356, -0.42612, -0.30051, -0.42771],
    [0.06042, 0.07473, 0.08904, 0.12879, -0.09222, 0.11607, 0.13515, 0.09699, 0.04611, -0.2226, -0.28143, 0.15582, -0.16377, -0.29097, 0.08586, 0.12402, 0.05883, -0.13356, -0.1113, -0.1272],
    [0.12084, 0.0, 0.12084, 0.19875, 0.0, 0.13515, 0.21624, 0.16695, 0.01908, -0.159, -0.20988, 0.07473, -0.09858, -0.20511, 0.16059, 0.12561, 0.08427, -0.11448, -0.08268, -0.0636],
    [-0.00636, 0.08745, 0.08427, 0.10812, -0.14151, 0.09699, 0.16695, 0.00477, 0.01908, -0.24009, -0.30051, 0.17808, -0.17808, -0.29574, 0.0636, 0.07155, 0.03021, -0.18285, -0.11766, -0.17649],
    [-0.02226, 0.01749, 0.03021, -0.00795, -0.21147, 0.04611, 0.01908, 0.01908, -0.12402, -0.29733, -0.36093, 0.14628, -0.27189, -0.3975, 0.00318, 0.02544, -0.02385, -0.27189, -0.19875, -0.20829],
    [-0.36729, -0.21624, -0.15423, -0.1431, -0.51357, -0.2226, -0.159, -0.24009, -0.29733, -0.67893, -0.75843, -0.11766, -0.59625, -0.72663, -0.23691, -0.19875, -0.27984, -0.55809, -0.47382, -0.60102],
    [-0.41976, -0.27984, -0.23373, -0.17967, -0.56604, -0.28143, -0.20988, -0.30051, -0.36093, -0.75843, -0.8109, -0.1749, -0.65826, -0.79659, -0.30687, -0.26235, -0.32913, -0.61533, -0.5406, -0.66939],
    [0.15264, 0.26712, 0.16854, 0.09381, 0.05088, 0.15582, 0.07473, 0.17808, 0.14628, -0.11766, -0.1749, 0.34185, -0.03339, -0.17331, 0.2067, 0.19398, 0.15264, -0.06678, -0.05247, -0.03498],
    [-0.26553, -0.13515, -0.10812, -0.0477, -0.43248, -0.16377, -0.09858, -0.17808, -0.27189, -0.59625, -0.65826, -0.03339, -0.50721, -0.68211, -0.18762, -0.12084, -0.19716, -0.52152, -0.41976, -0.48495],
    [-0.40386, -0.27189, -0.23532, -0.19239, -0.56127, -0.29097, -0.20511, -0.29574, -0.3975, -0.72663, -0.79659, -0.17331, -0.68211, -0.79341, -0.31482, -0.27825, -0.31959, -0.61851, -0.53901, -0.63918],
    [0.03816, 0.09063, 0.11766, 0.14946, -0.1272, 0.08586, 0.16059, 0.0636, 0.00318, -0.23691, -0.30687, 0.2067, -0.18762, -0.31482, 0.08268, 0.1113, 0.05883, -0.23214, -0.14628, -0.16695],
    [0.04134, 0.10335, 0.10971, 0.10176, -0.09381, 0.12402, 0.12561, 0.07155, 0.02544, -0.19875, -0.26235, 0.19398, -0.12084, -0.27825, 0.1113, 0.0954, 0.04929, -0.11448, -0.08109, -0.12402],
    [-0.00795, 0.05883, 0.06201, 0.07473, -0.13356, 0.05883, 0.08427, 0.03021, -0.02385, -0.27984, -0.32913, 0.15264, -0.19716, -0.31959, 0.05883, 0.04929, 0.02385, -0.15105, -0.11766, -0.18921],
    [-0.24645, -0.18126, -0.1272, -0.09063, -0.42612, -0.13356, -0.11448, -0.18285, -0.27189, -0.55809, -0.61533, -0.06678, -0.52152, -0.61851, -0.23214, -0.11448, -0.15105, -0.44361, -0.38001, -0.46269],
    [-0.17331, -0.14151, -0.07791, -0.07791, -0.30051, -0.1113, -0.08268, -0.11766, -0.19875, -0.47382, -0.5406, -0.05247, -0.41976, -0.53901, -0.14628, -0.08109, -0.11766, -0.38001, -0.3021, -0.37365],
    [-0.28143, -0.1272, -0.08904, -0.03339, -0.42771, -0.1272, -0.0636, -0.17649, -0.20829, -0.60102, -0.66939, -0.03498, -0.48495, -0.63918, -0.16695, -0.12402, -0.18921, -0.46269, -0.37365, -0.51675],
];

impl ForceField for KimHummer {
    fn params(&self, res_name: &str, bead_type: BeadType) -> Option<BeadParams> {
        match bead_type {
            BeadType::Residue | BeadType::Titratable => RESIDUES
                .iter()
                .find(|(name, _, _)| *name == res_name)
                .map(|(name, mass, sigma)| BeadParams {
                    mass: *mass,
                    sigma: *sigma,
                    epsilon: 0.0,
                    lambda: 0.0,
                    residue: Some(name),
                }),
            BeadType::Ion | BeadType::Virtual | BeadType::Ntr | BeadType::Ctr => Some(SITE),
        }
    }

    /// Only σ — KH epsilon is pair-specific (from the Miyazawa-Jernigan matrix),
    /// not a per-atom property.
    fn format_atom_fields(&self, p: &BeadParams) -> String {
        format!("σ: {}", p.sigma)
    }

    fn nonbonded_yaml(&self, types: &[(&str, f64, BeadParams)]) -> String {
        let mut yaml = "\nsystem:\n  energy:\n    nonbonded:\n".to_string();

        // Neutral pairs go under `replace:` to skip redundant Coulomb(0×0) evaluation.
        // Charged pairs go under `append:` to inherit Coulomb from the input's default.
        let mut replace_pairs = Vec::new();
        let mut append_pairs = Vec::new();

        for (i, (na, qa, pa)) in types.iter().enumerate() {
            for (nb, qb, pb) in &types[i..] {
                let sigma = (pa.sigma + pb.sigma) / 2.0;
                let mut epsilon = pair_epsilon(pa, pb);

                if let HydrophobicScaling::ScaleEpsilon(c) = &self.scaling
                    && pa.is_hydrophobic()
                    && pb.is_hydrophobic()
                {
                    epsilon *= c;
                }

                let entry = format!(
                    "        [{na}, {nb}]:\n          \
                     - !KimHummer {{sigma: {sigma}, epsilon: {epsilon:.4}}}\n"
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

/// Get pair epsilon in kJ/mol. Falls back to SITE epsilon for non-residue types.
fn pair_epsilon(a: &BeadParams, b: &BeadParams) -> f64 {
    match (
        a.residue.and_then(KimHummer::residue_index),
        b.residue.and_then(KimHummer::residue_index),
    ) {
        (Some(i), Some(j)) => EPSILON_KT[i][j] * KT_TO_KJ_MOL,
        _ => SITE.epsilon,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn residue_lookup() {
        let kh = KimHummer::new(HydrophobicScaling::NoScale);
        let p = kh.params("ALA", BeadType::Residue).unwrap();
        assert!((p.sigma - 5.0).abs() < 1e-10);
        assert!((p.mass - 71.07).abs() < 1e-10);
        assert!(p.epsilon == 0.0);
        assert!(p.lambda == 0.0);
    }

    #[test]
    fn site_lookup() {
        let kh = KimHummer::new(HydrophobicScaling::NoScale);
        let p = kh.params("ASP", BeadType::Virtual).unwrap();
        assert!((p.sigma - 2.0).abs() < 1e-10);
        assert!((p.epsilon - 0.1).abs() < 1e-10);
    }

    #[test]
    fn epsilon_symmetric() {
        for (i, row) in EPSILON_KT.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - EPSILON_KT[j][i]).abs() < 1e-10,
                    "epsilon matrix not symmetric at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn epsilon_ala_ala() {
        let kh = KimHummer::new(HydrophobicScaling::NoScale);
        let ala = kh.params("ALA", BeadType::Residue).unwrap();
        let eps = pair_epsilon(&ala, &ala);
        assert!((eps - (-0.07155 * KT_TO_KJ_MOL)).abs() < 1e-10);
    }

    #[test]
    fn epsilon_site_fallback() {
        let kh = KimHummer::new(HydrophobicScaling::NoScale);
        let ala = kh.params("ALA", BeadType::Residue).unwrap();
        let eps = pair_epsilon(&ala, &SITE);
        assert!((eps - SITE.epsilon).abs() < 1e-10);
    }

    #[test]
    fn pair_epsilon_follows_residue_not_type_name() {
        // Variant (CYX) and clustered (TYR1) type names must get the parent's MJ epsilon.
        let kh = KimHummer::new(HydrophobicScaling::NoScale);
        let cys = kh.params("CYS", BeadType::Residue).unwrap();
        let tyr = kh.params("TYR", BeadType::Titratable).unwrap();
        let yaml = kh.nonbonded_yaml(&[("CYX", 0.0, cys), ("TYR1", -0.1, tyr)]);
        let expected = format!(
            "[CYX, TYR1]:\n          - !KimHummer {{sigma: {}, epsilon: {:.4}}}",
            (cys.sigma + tyr.sigma) / 2.0,
            pair_epsilon(&cys, &tyr)
        );
        assert!(yaml.contains(&expected), "{yaml}");
        assert!((pair_epsilon(&cys, &tyr) - SITE.epsilon).abs() > 1e-6);
    }

    #[test]
    fn neutral_pairs_in_replace() {
        let kh = KimHummer::new(HydrophobicScaling::NoScale);
        let types = vec![
            ("ALA", 0.0, kh.params("ALA", BeadType::Residue).unwrap()),
            ("VAL", 0.0, kh.params("VAL", BeadType::Residue).unwrap()),
        ];
        let yaml = kh.nonbonded_yaml(&types);
        assert!(yaml.contains("replace:"));
        assert!(!yaml.contains("append:"));
    }

    #[test]
    fn charged_pairs_in_append() {
        let kh = KimHummer::new(HydrophobicScaling::NoScale);
        let types = vec![
            ("ALA", 0.0, kh.params("ALA", BeadType::Residue).unwrap()),
            ("ARG", 1.0, kh.params("ARG", BeadType::Residue).unwrap()),
        ];
        let yaml = kh.nonbonded_yaml(&types);
        // ALA-ALA neutral → replace, ALA-ARG and ARG-ARG charged → append
        assert!(yaml.contains("replace:"));
        assert!(yaml.contains("append:"));
    }

    #[test]
    fn hydrophobic_epsilon_scaling() {
        let kh = KimHummer::new(HydrophobicScaling::ScaleEpsilon(0.5));
        let types = vec![
            ("ALA", 0.0, kh.params("ALA", BeadType::Residue).unwrap()),
            ("GLU", -1.0, kh.params("GLU", BeadType::Residue).unwrap()),
        ];
        let yaml = kh.nonbonded_yaml(&types);
        // ALA-ALA is hydrophobic: epsilon should be scaled by 0.5
        let expected_eps = -0.07155 * KT_TO_KJ_MOL * 0.5;
        assert!(yaml.contains(&format!("{expected_eps:.4}")));
        // ALA-GLU is NOT hydrophobic: epsilon unscaled
        let unscaled = 0.12084 * KT_TO_KJ_MOL;
        assert!(yaml.contains(&format!("{unscaled:.4}")));
    }

    #[test]
    fn lambda_scaling_rejected() {
        // ScaleLambda should not affect KH (no lambda to scale)
        let kh = KimHummer::new(HydrophobicScaling::ScaleLambda(2.0));
        let types = vec![("ALA", 0.0, kh.params("ALA", BeadType::Residue).unwrap())];
        let yaml = kh.nonbonded_yaml(&types);
        // Epsilon should be unscaled
        let expected = -0.07155 * KT_TO_KJ_MOL;
        assert!(yaml.contains(&format!("{expected:.4}")));
    }
}
