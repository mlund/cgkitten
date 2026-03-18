//! PDB file format parser.
//!
//! Parses ATOM/HETATM records and SSBOND disulfide bonds from PDB files.
//! Adapted from the voronota-ltr project (MIT License).

use std::io::BufRead;

use crate::mmcif::{
    AtomRecord, DisulfideBond, ParseOptions, ParsedMmcif, is_hydrogen, normalize_atom_name,
};

/// Extract a substring from fixed-width PDB columns (1-indexed, inclusive).
fn col(line: &str, start: usize, end: usize) -> &str {
    let len = line.len();
    let s = start.saturating_sub(1);
    let e = end.min(len);
    if s >= len {
        return "";
    }
    line.get(s..e).unwrap_or("").trim()
}

fn col_f64(line: &str, start: usize, end: usize) -> Option<f64> {
    let s = col(line, start, end);
    if s.is_empty() { None } else { s.parse().ok() }
}

fn col_i32(line: &str, start: usize, end: usize) -> Option<i32> {
    let s = col(line, start, end);
    if s.is_empty() { None } else { s.parse().ok() }
}

fn parse_atom_line(line: &str) -> Option<AtomRecord> {
    let record_name = col(line, 1, 6);
    if record_name != "ATOM" && record_name != "HETATM" {
        return None;
    }

    let x = col_f64(line, 31, 38)?;
    let y = col_f64(line, 39, 46)?;
    let z = col_f64(line, 47, 54)?;
    let res_seq = col_i32(line, 23, 26)?;

    let raw_name = col(line, 13, 16);
    if raw_name.is_empty() {
        return None;
    }

    let alt_loc = col(line, 17, 17).to_string();

    Some(AtomRecord {
        record_name: record_name.to_string(),
        name: normalize_atom_name(raw_name),
        alt_loc,
        res_name: col(line, 18, 20).to_string(),
        chain_id: col(line, 22, 22).to_string(),
        res_seq,
        i_code: col(line, 27, 27).to_string(),
        x,
        y,
        z,
        element: col(line, 77, 78).to_string(),
    })
}

fn parse_ssbond_line(line: &str) -> Option<DisulfideBond> {
    // SSBOND cols: chain1=16, seqNum1=18-21, chain2=30, seqNum2=32-35
    let chain_id_1 = col(line, 16, 16).to_string();
    let res_seq_1 = col_i32(line, 18, 21)?;
    let chain_id_2 = col(line, 30, 30).to_string();
    let res_seq_2 = col_i32(line, 32, 35)?;
    Some(DisulfideBond {
        chain_id_1,
        res_seq_1,
        chain_id_2,
        res_seq_2,
    })
}

fn is_acceptable(record: &AtomRecord, options: &ParseOptions) -> bool {
    // AltLoc: accept empty, "A", or "1"
    if !record.alt_loc.is_empty() && record.alt_loc != "A" && record.alt_loc != "1" {
        return false;
    }
    if !options.include_hydrogens && is_hydrogen(record) {
        return false;
    }
    true
}

/// Parse ATOM/HETATM records and SSBOND disulfide bonds from a PDB file.
/// Only the first MODEL is read.
pub fn parse_pdb<R: BufRead>(reader: R, options: &ParseOptions) -> ParsedMmcif {
    let mut atoms = Vec::new();
    let mut disulfide_bonds = Vec::new();

    for line in reader.lines().map_while(Result::ok) {
        let tag = col(&line, 1, 6);
        match tag {
            "ATOM" | "HETATM" => {
                if let Some(mut record) = parse_atom_line(&line)
                    && is_acceptable(&record, options)
                {
                    record.alt_loc.clear();
                    atoms.push(record);
                }
            }
            "SSBOND" => {
                if let Some(bond) = parse_ssbond_line(&line) {
                    disulfide_bonds.push(bond);
                }
            }
            "ENDMDL" | "END" => break,
            _ => {}
        }
    }

    ParsedMmcif {
        atoms,
        disulfide_bonds,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE_PDB: &str = "\
SSBOND   1 CYS A    6    CYS A  127
ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  MET A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   MET A   1       2.009   1.420   0.000  1.00  0.00           C
HETATM   99  ZN  ZN  A 100       5.000   5.000   5.000  1.00  0.00          ZN
END
";

    #[test]
    fn parse_atom_records() {
        let options = ParseOptions::default();
        let parsed = parse_pdb(SAMPLE_PDB.as_bytes(), &options);
        assert_eq!(parsed.atoms.len(), 4);
        assert_eq!(parsed.atoms[0].name, "N");
        assert_eq!(parsed.atoms[0].res_name, "MET");
        assert_eq!(parsed.atoms[1].name, "CA");
        assert_eq!(parsed.atoms[3].name, "ZN");
    }

    #[test]
    fn parse_ssbond() {
        let options = ParseOptions::default();
        let parsed = parse_pdb(SAMPLE_PDB.as_bytes(), &options);
        assert_eq!(parsed.disulfide_bonds.len(), 1);
        let bond = &parsed.disulfide_bonds[0];
        assert_eq!(bond.chain_id_1, "A");
        assert_eq!(bond.res_seq_1, 6);
        assert_eq!(bond.chain_id_2, "A");
        assert_eq!(bond.res_seq_2, 127);
    }

    #[test]
    fn normalize_numbered_atom() {
        assert_eq!(normalize_atom_name("1HG2"), "HG21");
        assert_eq!(normalize_atom_name("CA"), "CA");
        assert_eq!(normalize_atom_name("2HD1"), "HD12");
    }

    #[test]
    fn stops_at_endmdl() {
        let cif = "ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00  0.00           N\nENDMDL\nATOM      2  CA  MET A   1       1.458   0.000   0.000  1.00  0.00           C\n";
        let parsed = parse_pdb(cif.as_bytes(), &ParseOptions::default());
        assert_eq!(parsed.atoms.len(), 1);
    }
}
