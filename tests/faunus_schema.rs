//! Generated topologies are written by hand-formatted strings, so nothing ties them to
//! faunus' input schema. Parse them with faunus here so schema changes fail CI, not users.

use cgkitten::forcefield::HydrophobicScaling;
use faunus::energy::NonbondedMatrix;
use faunus::interatomic::coulomb::{Medium, Salt};
use faunus::topology::Topology;
use std::path::Path;

/// Hydrophobic and charged residues, so `replace:` pair blocks are emitted too.
const CIF: &str = r#"
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
ATOM 1 N LEU A 1 0.000 0.000 0.000 N . 1
ATOM 2 CA LEU A 1 1.458 0.000 0.000 C . 1
ATOM 3 C LEU A 1 2.400 -1.100 0.000 C . 1
ATOM 4 O LEU A 1 2.200 -2.300 0.000 O . 1
ATOM 5 CB LEU A 1 2.000 1.000 0.000 C . 1
ATOM 6 N PHE A 2 3.600 -0.500 0.000 N . 1
ATOM 7 CA PHE A 2 4.900 -1.100 0.000 C . 1
ATOM 8 C PHE A 2 5.800 -2.200 0.000 C . 1
ATOM 9 O PHE A 2 5.600 -3.400 0.000 O . 1
ATOM 10 CB PHE A 2 5.500 0.100 0.800 C . 1
ATOM 11 N LYS A 3 7.000 -1.700 0.000 N . 1
ATOM 12 CA LYS A 3 8.300 -2.300 0.000 C . 1
ATOM 13 C LYS A 3 9.200 -3.400 0.000 C . 1
ATOM 14 O LYS A 3 9.000 -4.600 0.000 O . 1
ATOM 15 NZ LYS A 3 9.000 -1.000 1.500 N . 1
ATOM 16 N GLU A 4 10.400 -2.900 0.000 N . 1
ATOM 17 CA GLU A 4 11.700 -3.500 0.000 C . 1
ATOM 18 C GLU A 4 12.600 -4.600 0.000 C . 1
ATOM 19 O GLU A 4 12.400 -5.800 0.000 O . 1
ATOM 20 OE1 GLU A 4 12.300 -2.300 1.200 O . 1
"#;

fn assert_faunus_parses(model: &str, multi_bead: bool) {
    let policy: &dyn cgkitten::CoarseGrain = if multi_bead {
        &cgkitten::MultiBead
    } else {
        &cgkitten::SingleBead
    };
    // No MC sweeps: deterministic charges are enough to exercise the schema.
    let yaml = cgkitten::coarse_grain_to_strings(
        CIF.as_bytes(),
        false,
        &cgkitten::ChargeCalc::new().mc(0),
        model,
        policy,
        0.02,
        HydrophobicScaling::NoScale,
        &[] as &[&str],
    )
    .unwrap()
    .topology;

    let ctx = format!("{model} (multi_bead={multi_bead}):\n{yaml}");
    let mut topology = Topology::from_str_partial(&yaml).expect(&ctx);
    topology.finalize_atoms().expect(&ctx);
    topology.finalize_molecules(Path::new(".")).expect(&ctx);
    let medium = Medium::salt_water(298.15, Salt::SodiumChloride, 0.1);
    NonbondedMatrix::from_str(&yaml, &topology, Some(medium)).expect(&ctx);
}

#[test]
fn generated_topologies_parse_in_faunus() {
    for model in ["calvados3", "kimhummer", "pasquier"] {
        for multi_bead in [false, true] {
            assert_faunus_parses(model, multi_bead);
        }
    }
}
