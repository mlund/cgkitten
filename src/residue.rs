//! Amino acid definitions and titratable side-chain groups.
//!
//! Each amino acid is reduced to a single backbone bead at the center of mass.
//! Titratable residues get an additional side-chain bead placed at the charge center.

use crate::mmcif::AtomRecord;

/// The 20 standard amino acid three-letter codes.
pub const STANDARD_RESIDUES: &[&str] = &[
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
];

/// Returns true if residue name is a standard amino acid.
pub fn is_amino_acid(res_name: &str) -> bool {
    STANDARD_RESIDUES.contains(&res_name)
}

/// Known metal ions to keep as individual beads.
const METAL_IONS: &[&str] = &[
    "ZN", "FE", "MG", "CA", "MN", "CU", "CO", "NI", "CD", "NA", "K",
];

/// Returns true if the record is a coordinated metal ion.
pub fn is_metal_ion(record: &AtomRecord) -> bool {
    // Metal ions are typically HETATM with residue name == element symbol
    record.record_name == "HETATM" && METAL_IONS.contains(&record.res_name.as_str())
}

/// Atomic mass in Da from element symbol.
pub fn atomic_mass(element: &str) -> f64 {
    match element {
        "H" => 1.008,
        "C" => 12.011,
        "N" => 14.007,
        "O" => 15.999,
        "S" => 32.06,
        "SE" => 78.971,
        "P" => 30.974,
        // Both uppercase (mmCIF) and mixed case (PDB-style) element symbols
        "FE" | "Fe" => 55.845,
        "ZN" | "Zn" => 65.38,
        "MG" | "Mg" => 24.305,
        "CA" | "Ca" => 40.078,
        "MN" | "Mn" => 54.938,
        "CU" | "Cu" => 63.546,
        "CO" | "Co" => 58.933,
        "NI" | "Ni" => 58.693,
        "CD" | "Cd" => 112.414,
        "NA" | "Na" => 22.990,
        "K" => 39.098,
        _ => {
            log::warn!("Unknown element '{element}', using mass 0.0");
            0.0
        }
    }
}

/// A titratable group definition.
#[derive(Clone, Debug)]
pub struct TitratableGroup {
    /// Residue name (e.g., "ASP", "GLU", "HIS", "LYS", "ARG", "CYS", "TYR").
    pub res_name: &'static str,
    /// pKa value for Henderson-Hasselbalch calculation.
    pub pka: f64,
    /// Net charge of the protonated form.
    pub charge_protonated: f64,
    /// Net charge of the deprotonated form.
    pub charge_deprotonated: f64,
    /// Atom names that define the charge center (averaged for bead position).
    pub center_atoms: &'static [&'static str],
    /// Element symbol for the bead name (e.g., "O", "N", "S").
    pub element: &'static str,
}

/// Standard titratable groups with canonical pKa values.
pub static TITRATABLE_GROUPS: &[TitratableGroup] = &[
    TitratableGroup {
        res_name: "ASP",
        pka: 3.65,
        charge_protonated: 0.0,
        charge_deprotonated: -1.0,
        center_atoms: &["OD1", "OD2"],
        element: "O",
    },
    TitratableGroup {
        res_name: "GLU",
        pka: 4.25,
        charge_protonated: 0.0,
        charge_deprotonated: -1.0,
        center_atoms: &["OE1", "OE2"],
        element: "O",
    },
    TitratableGroup {
        res_name: "HIS",
        pka: 6.00,
        charge_protonated: 1.0,
        charge_deprotonated: 0.0,
        center_atoms: &["ND1", "NE2"],
        element: "N",
    },
    TitratableGroup {
        res_name: "CYS",
        pka: 8.18,
        charge_protonated: 0.0,
        charge_deprotonated: -1.0,
        center_atoms: &["SG"],
        element: "S",
    },
    TitratableGroup {
        res_name: "TYR",
        pka: 10.07,
        charge_protonated: 0.0,
        charge_deprotonated: -1.0,
        center_atoms: &["OH"],
        element: "O",
    },
    TitratableGroup {
        res_name: "LYS",
        pka: 10.54,
        charge_protonated: 1.0,
        charge_deprotonated: 0.0,
        center_atoms: &["NZ"],
        element: "N",
    },
    TitratableGroup {
        res_name: "ARG",
        pka: 12.48,
        charge_protonated: 1.0,
        charge_deprotonated: 0.0,
        center_atoms: &["NH1", "NH2"],
        element: "N",
    },
];

/// N-terminal group.
pub static NTERM_GROUP: TitratableGroup = TitratableGroup {
    res_name: "NTERM",
    pka: 8.0,
    charge_protonated: 1.0,
    charge_deprotonated: 0.0,
    center_atoms: &["N"],
    element: "N",
};

/// C-terminal group.
pub static CTERM_GROUP: TitratableGroup = TitratableGroup {
    res_name: "CTERM",
    pka: 3.1,
    charge_protonated: 0.0,
    charge_deprotonated: -1.0,
    center_atoms: &["OXT", "O"],
    element: "O",
};

/// Look up the titratable group for a residue, if any.
pub fn find_titratable_group(res_name: &str) -> Option<&'static TitratableGroup> {
    TITRATABLE_GROUPS.iter().find(|g| g.res_name == res_name)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_residues_recognized() {
        assert!(is_amino_acid("ALA"));
        assert!(is_amino_acid("GLU"));
        assert!(!is_amino_acid("HOH"));
        assert!(!is_amino_acid("ZN"));
    }

    #[test]
    fn titratable_groups_found() {
        assert!(find_titratable_group("ASP").is_some());
        assert!(find_titratable_group("GLU").is_some());
        assert!(find_titratable_group("HIS").is_some());
        assert!(find_titratable_group("LYS").is_some());
        assert!(find_titratable_group("ARG").is_some());
        assert!(find_titratable_group("CYS").is_some());
        assert!(find_titratable_group("TYR").is_some());
        assert!(find_titratable_group("ALA").is_none());
    }

    #[test]
    fn metal_ions_detected() {
        let zn = AtomRecord {
            record_name: "HETATM".into(),
            name: "ZN".into(),
            alt_loc: String::new(),
            res_name: "ZN".into(),
            chain_id: "A".into(),
            res_seq: 100,
            i_code: String::new(),
            x: 0.0,
            y: 0.0,
            z: 0.0,
            element: "ZN".into(),
        };
        assert!(is_metal_ion(&zn));
    }
}
