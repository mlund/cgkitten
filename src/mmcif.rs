//! mmCIF file format parser.
//!
//! Parses `_atom_site` and `_struct_conn` loop blocks from mmCIF files.

use std::collections::HashMap;
use std::io::BufRead;

/// Parsed atom record with coordinates and metadata.
#[derive(Clone, Debug)]
pub struct AtomRecord {
    pub record_name: String,
    pub name: String,
    pub alt_loc: String,
    pub res_name: String,
    pub chain_id: String,
    pub res_seq: i32,
    pub i_code: String,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub element: String,
}

/// A disulfide bond between two residues, parsed from `_struct_conn`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DisulfideBond {
    pub chain_id_1: String,
    pub res_seq_1: i32,
    pub chain_id_2: String,
    pub res_seq_2: i32,
}

/// Result of parsing an mmCIF file.
#[derive(Clone, Debug)]
pub struct ParsedMmcif {
    pub atoms: Vec<AtomRecord>,
    pub disulfide_bonds: Vec<DisulfideBond>,
}

/// Options for parsing molecular files.
#[derive(Clone, Debug, Default)]
pub struct ParseOptions {
    pub include_hydrogens: bool,
}

// ---------------------------------------------------------------------------
// _atom_site fields
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum AtomField {
    ModelNum,
    GroupPdb,
    Id,
    AuthAtomId,
    LabelAtomId,
    LabelAltId,
    AuthCompId,
    LabelCompId,
    AuthAsymId,
    LabelAsymId,
    AuthSeqId,
    LabelSeqId,
    InsCode,
    CartnX,
    CartnY,
    CartnZ,
    TypeSymbol,
}

fn atom_field_from_name(name: &str) -> Option<AtomField> {
    match name {
        "pdbx_PDB_model_num" => Some(AtomField::ModelNum),
        "group_PDB" => Some(AtomField::GroupPdb),
        "id" => Some(AtomField::Id),
        "auth_atom_id" => Some(AtomField::AuthAtomId),
        "label_atom_id" => Some(AtomField::LabelAtomId),
        "label_alt_id" => Some(AtomField::LabelAltId),
        "auth_comp_id" => Some(AtomField::AuthCompId),
        "label_comp_id" => Some(AtomField::LabelCompId),
        "auth_asym_id" => Some(AtomField::AuthAsymId),
        "label_asym_id" => Some(AtomField::LabelAsymId),
        "auth_seq_id" => Some(AtomField::AuthSeqId),
        "label_seq_id" => Some(AtomField::LabelSeqId),
        "pdbx_PDB_ins_code" => Some(AtomField::InsCode),
        "Cartn_x" => Some(AtomField::CartnX),
        "Cartn_y" => Some(AtomField::CartnY),
        "Cartn_z" => Some(AtomField::CartnZ),
        "type_symbol" => Some(AtomField::TypeSymbol),
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// _struct_conn fields
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum ConnField {
    Type,
    Ptnr1AuthAsym,
    Ptnr1LabelAsym,
    Ptnr1AuthSeq,
    Ptnr1LabelSeq,
    Ptnr2AuthAsym,
    Ptnr2LabelAsym,
    Ptnr2AuthSeq,
    Ptnr2LabelSeq,
}

fn conn_field_from_name(name: &str) -> Option<ConnField> {
    match name {
        "conn_type_id" => Some(ConnField::Type),
        "ptnr1_auth_asym_id" => Some(ConnField::Ptnr1AuthAsym),
        "ptnr1_label_asym_id" => Some(ConnField::Ptnr1LabelAsym),
        "ptnr1_auth_seq_id" => Some(ConnField::Ptnr1AuthSeq),
        "ptnr1_label_seq_id" => Some(ConnField::Ptnr1LabelSeq),
        "ptnr2_auth_asym_id" => Some(ConnField::Ptnr2AuthAsym),
        "ptnr2_label_asym_id" => Some(ConnField::Ptnr2LabelAsym),
        "ptnr2_auth_seq_id" => Some(ConnField::Ptnr2AuthSeq),
        "ptnr2_label_seq_id" => Some(ConnField::Ptnr2LabelSeq),
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Tokenizer
// ---------------------------------------------------------------------------

struct Tokenizer<R> {
    reader: R,
    buffer: String,
    position: usize,
    eof: bool,
}

impl<R: BufRead> Tokenizer<R> {
    const fn new(reader: R) -> Self {
        Self {
            reader,
            buffer: String::new(),
            position: 0,
            eof: false,
        }
    }

    fn ensure_data(&mut self) -> bool {
        while self.position >= self.buffer.len() && !self.eof {
            self.buffer.clear();
            self.position = 0;
            match self.reader.read_line(&mut self.buffer) {
                Ok(0) => {
                    self.eof = true;
                    return false;
                }
                Err(e) => {
                    log::warn!("I/O error reading mmCIF input: {e}");
                    self.eof = true;
                    return false;
                }
                Ok(_) => {}
            }
        }
        !self.eof || self.position < self.buffer.len()
    }

    fn skip_whitespace_and_comments(&mut self) {
        loop {
            if !self.ensure_data() {
                return;
            }
            let bytes = self.buffer.as_bytes();
            while self.position < bytes.len() && bytes[self.position].is_ascii_whitespace() {
                self.position += 1;
            }
            if self.position >= bytes.len() {
                continue;
            }
            if bytes[self.position] == b'#' {
                self.position = bytes.len();
                continue;
            }
            break;
        }
    }

    fn next_token(&mut self) -> Option<String> {
        self.skip_whitespace_and_comments();
        if !self.ensure_data() {
            return None;
        }
        let bytes = self.buffer.as_bytes();
        if self.position >= bytes.len() {
            return None;
        }
        let ch = bytes[self.position];

        if ch == b'\'' || ch == b'"' {
            let quote = ch;
            self.position += 1;
            let start = self.position;
            while self.position < bytes.len() && bytes[self.position] != quote {
                self.position += 1;
            }
            let token = String::from_utf8_lossy(&bytes[start..self.position]).into_owned();
            if self.position < bytes.len() {
                self.position += 1;
            }
            return Some(token);
        }

        let start = self.position;
        while self.position < bytes.len() && !bytes[self.position].is_ascii_whitespace() {
            self.position += 1;
        }
        Some(String::from_utf8_lossy(&bytes[start..self.position]).into_owned())
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn fix_undefined(s: &str) -> String {
    if s == "." || s == "?" {
        String::new()
    } else {
        s.to_string()
    }
}

/// Check if token signals end of a loop data block.
fn is_stop_token(token: &str) -> bool {
    token.starts_with('_') || token.starts_with("data_") || token == "loop_"
}

pub(crate) fn normalize_atom_name(name: &str) -> String {
    if name.is_empty() {
        return String::new();
    }
    let first_letter = name.find(|c: char| !c.is_ascii_digit());
    match first_letter {
        Some(0) | None => name.to_string(),
        Some(pos) => format!("{}{}", &name[pos..], &name[..pos]),
    }
}

fn get_value<'a, K: Eq + std::hash::Hash>(
    row: &'a [String],
    field_map: &HashMap<K, usize>,
    primary: K,
    fallback: Option<K>,
) -> &'a str {
    if let Some(&idx) = field_map.get(&primary)
        && idx < row.len()
    {
        return &row[idx];
    }
    if let Some(fb) = fallback
        && let Some(&idx) = field_map.get(&fb)
        && idx < row.len()
    {
        return &row[idx];
    }
    ""
}

// ---------------------------------------------------------------------------
// _atom_site parsing
// ---------------------------------------------------------------------------

fn parse_atom_row(row: &[String], field_map: &HashMap<AtomField, usize>) -> Option<AtomRecord> {
    let x: f64 = get_value(row, field_map, AtomField::CartnX, None)
        .parse()
        .ok()?;
    let y: f64 = get_value(row, field_map, AtomField::CartnY, None)
        .parse()
        .ok()?;
    let z: f64 = get_value(row, field_map, AtomField::CartnZ, None)
        .parse()
        .ok()?;
    let res_seq: i32 = get_value(
        row,
        field_map,
        AtomField::AuthSeqId,
        Some(AtomField::LabelSeqId),
    )
    .parse()
    .ok()?;

    let raw_name = get_value(
        row,
        field_map,
        AtomField::AuthAtomId,
        Some(AtomField::LabelAtomId),
    );
    if raw_name.is_empty() {
        return None;
    }
    let record_name = get_value(row, field_map, AtomField::GroupPdb, None);
    if record_name.is_empty() {
        return None;
    }

    Some(AtomRecord {
        record_name: record_name.to_string(),
        name: normalize_atom_name(raw_name),
        alt_loc: fix_undefined(get_value(row, field_map, AtomField::LabelAltId, None)),
        res_name: fix_undefined(get_value(
            row,
            field_map,
            AtomField::AuthCompId,
            Some(AtomField::LabelCompId),
        )),
        chain_id: fix_undefined(get_value(
            row,
            field_map,
            AtomField::AuthAsymId,
            Some(AtomField::LabelAsymId),
        )),
        res_seq,
        i_code: fix_undefined(get_value(row, field_map, AtomField::InsCode, None)),
        x,
        y,
        z,
        element: fix_undefined(get_value(row, field_map, AtomField::TypeSymbol, None)),
    })
}

pub(crate) fn is_hydrogen(record: &AtomRecord) -> bool {
    record.name.starts_with('H') || record.element == "H" || record.element == "D"
}

fn is_acceptable(record: &AtomRecord, options: &ParseOptions) -> bool {
    if !record.alt_loc.is_empty()
        && record.alt_loc != "A"
        && record.alt_loc != "1"
        && record.alt_loc != "."
    {
        return false;
    }
    if !options.include_hydrogens && is_hydrogen(record) {
        return false;
    }
    if record.res_name == "HOH" || record.res_name == "WAT" || record.res_name == "DOD" {
        return false;
    }
    true
}

/// Returns the stop token (if any) that ended the data block.
fn read_atom_site_data<R: BufRead>(
    tokenizer: &mut Tokenizer<R>,
    field_map: &HashMap<AtomField, usize>,
    num_cols: usize,
    first_token: &str,
    options: &ParseOptions,
    records: &mut Vec<AtomRecord>,
) -> Option<String> {
    let mut first_model_id: Option<String> = None;
    let mut row = Vec::with_capacity(num_cols);
    row.push(first_token.to_string());

    loop {
        // Check for end-of-loop before filling the row
        if is_stop_token(&row[0]) {
            return Some(row.swap_remove(0));
        }

        while row.len() < num_cols {
            match tokenizer.next_token() {
                Some(tok) => row.push(tok),
                None => return None,
            }
        }

        let model_id = get_value(&row, field_map, AtomField::ModelNum, None);
        let first = first_model_id.get_or_insert_with(|| model_id.to_string());
        let should_process = first == model_id;

        if should_process
            && let Some(mut record) = parse_atom_row(&row, field_map)
            && is_acceptable(&record, options)
        {
            record.alt_loc.clear();
            records.push(record);
        }

        row.clear();
        match tokenizer.next_token() {
            Some(tok) => row.push(tok),
            None => return None,
        }
    }
}

// ---------------------------------------------------------------------------
// _struct_conn parsing
// ---------------------------------------------------------------------------

fn parse_disulfide_row(
    row: &[String],
    field_map: &HashMap<ConnField, usize>,
) -> Option<DisulfideBond> {
    let conn_type = get_value(row, field_map, ConnField::Type, None);
    if conn_type != "disulf" {
        return None;
    }

    let chain1 = fix_undefined(get_value(
        row,
        field_map,
        ConnField::Ptnr1AuthAsym,
        Some(ConnField::Ptnr1LabelAsym),
    ));
    let seq1: i32 = get_value(
        row,
        field_map,
        ConnField::Ptnr1AuthSeq,
        Some(ConnField::Ptnr1LabelSeq),
    )
    .parse()
    .ok()?;

    let chain2 = fix_undefined(get_value(
        row,
        field_map,
        ConnField::Ptnr2AuthAsym,
        Some(ConnField::Ptnr2LabelAsym),
    ));
    let seq2: i32 = get_value(
        row,
        field_map,
        ConnField::Ptnr2AuthSeq,
        Some(ConnField::Ptnr2LabelSeq),
    )
    .parse()
    .ok()?;

    Some(DisulfideBond {
        chain_id_1: chain1,
        res_seq_1: seq1,
        chain_id_2: chain2,
        res_seq_2: seq2,
    })
}

/// Returns the stop token (if any) that ended the data block.
fn read_struct_conn_data<R: BufRead>(
    tokenizer: &mut Tokenizer<R>,
    field_map: &HashMap<ConnField, usize>,
    num_cols: usize,
    first_token: &str,
    bonds: &mut Vec<DisulfideBond>,
) -> Option<String> {
    let mut row = Vec::with_capacity(num_cols);
    row.push(first_token.to_string());

    loop {
        // Check for end-of-loop before filling the row
        if is_stop_token(&row[0]) {
            return Some(row.swap_remove(0));
        }

        while row.len() < num_cols {
            match tokenizer.next_token() {
                Some(tok) => row.push(tok),
                None => return None,
            }
        }

        if let Some(bond) = parse_disulfide_row(&row, field_map) {
            bonds.push(bond);
        }

        row.clear();
        match tokenizer.next_token() {
            Some(tok) => row.push(tok),
            None => return None,
        }
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Parse mmCIF format, extracting atom records and disulfide bonds.
pub fn parse_mmcif<R: BufRead>(reader: R, options: &ParseOptions) -> ParsedMmcif {
    let mut tokenizer = Tokenizer::new(reader);
    let mut atoms = Vec::new();
    let mut disulfide_bonds = Vec::new();
    let mut pending_token: Option<String> = None;

    loop {
        let token = if let Some(pt) = pending_token.take() {
            pt
        } else {
            match tokenizer.next_token() {
                Some(t) => t,
                None => break,
            }
        };

        if token != "loop_" {
            continue;
        }

        let Some(first_field) = tokenizer.next_token() else {
            break;
        };

        if first_field.starts_with("_atom_site.") {
            let (field_map, header_len, first_data) = read_loop_header(
                &mut tokenizer,
                &first_field,
                "_atom_site.",
                atom_field_from_name,
            );
            if let Some(first_data) = first_data {
                pending_token = read_atom_site_data(
                    &mut tokenizer,
                    &field_map,
                    header_len,
                    &first_data,
                    options,
                    &mut atoms,
                );
            }
        } else if first_field.starts_with("_struct_conn.") {
            let (field_map, header_len, first_data) = read_loop_header(
                &mut tokenizer,
                &first_field,
                "_struct_conn.",
                conn_field_from_name,
            );
            if let Some(first_data) = first_data {
                pending_token = read_struct_conn_data(
                    &mut tokenizer,
                    &field_map,
                    header_len,
                    &first_data,
                    &mut disulfide_bonds,
                );
            }
        }
    }

    ParsedMmcif {
        atoms,
        disulfide_bonds,
    }
}

/// Read loop header fields. Returns (field_map, column_count, first_data_token).
fn read_loop_header<R: BufRead, K: Eq + std::hash::Hash>(
    tokenizer: &mut Tokenizer<R>,
    first_field: &str,
    prefix: &str,
    field_parser: fn(&str) -> Option<K>,
) -> (HashMap<K, usize>, usize, Option<String>) {
    let mut field_map = HashMap::new();
    let mut count = 1;

    if let Some(field) = field_parser(&first_field[prefix.len()..]) {
        field_map.insert(field, 0);
    }

    loop {
        let Some(tok) = tokenizer.next_token() else {
            return (field_map, count, None);
        };
        if !tok.starts_with(prefix) {
            return (field_map, count, Some(tok));
        }
        if let Some(field) = field_parser(&tok[prefix.len()..]) {
            field_map.insert(field, count);
        }
        count += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE_MMCIF: &str = r#"
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
HETATM 99 ZN ZN A 100 5.000 5.000 5.000 ZN . 1
HETATM 100 O HOH B 200 10.000 10.000 10.000 O . 1
"#;

    #[test]
    fn parse_atom_records() {
        let options = ParseOptions::default();
        let parsed = parse_mmcif(SAMPLE_MMCIF.as_bytes(), &options);
        assert!(parsed.atoms.iter().all(|r| r.res_name != "HOH"));
        assert!(parsed.atoms.iter().any(|r| r.element == "ZN"));
    }

    #[test]
    fn water_is_removed() {
        let options = ParseOptions::default();
        let parsed = parse_mmcif(SAMPLE_MMCIF.as_bytes(), &options);
        assert!(!parsed.atoms.iter().any(|r| r.res_name == "HOH"));
    }

    #[test]
    fn parse_struct_conn_disulfides() {
        let cif = r#"
data_test
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_seq_id
disulf1 disulf A 6 A 127
disulf2 disulf A 30 A 115
covale1 covale A 50 B 1
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
ATOM 1 CA ALA A 1 0.0 0.0 0.0 C . 1
"#;
        let parsed = parse_mmcif(cif.as_bytes(), &ParseOptions::default());
        assert_eq!(parsed.disulfide_bonds.len(), 2);
        assert_eq!(
            parsed.disulfide_bonds[0],
            DisulfideBond {
                chain_id_1: "A".into(),
                res_seq_1: 6,
                chain_id_2: "A".into(),
                res_seq_2: 127,
            }
        );
        assert_eq!(
            parsed.disulfide_bonds[1],
            DisulfideBond {
                chain_id_1: "A".into(),
                res_seq_1: 30,
                chain_id_2: "A".into(),
                res_seq_2: 115,
            }
        );
        assert_eq!(parsed.atoms.len(), 1);
    }
}
