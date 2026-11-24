//! Generated topologies are written by hand-formatted strings, so nothing ties them to
//! faunus' input schema. Parse them with faunus here so schema changes fail CI, not users.

use cgkitten::forcefield::{self, ForceField, HydrophobicScaling};
use faunus::energy::{NonbondedMatrix, PairPot};
use faunus::interatomic::{
    Cutoff,
    coulomb::{Medium, Salt},
    twobody::IsotropicTwobodyEnergy,
};
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

fn topology_yaml(model: &str, multi_bead: bool) -> String {
    let policy: &dyn cgkitten::CoarseGrain = if multi_bead {
        &cgkitten::MultiBead
    } else {
        &cgkitten::SingleBead
    };
    // No MC sweeps: deterministic charges are enough to exercise the schema.
    cgkitten::coarse_grain_to_strings(
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
    .topology
}

fn nonbonded(yaml: &str, medium: Option<Medium>) -> NonbondedMatrix {
    let mut topology = Topology::from_str_partial(yaml).expect(yaml);
    topology.finalize_atoms().expect(yaml);
    topology.finalize_molecules(Path::new(".")).expect(yaml);
    NonbondedMatrix::from_str(yaml, &topology, medium).expect(yaml)
}

fn assert_faunus_parses_yaml(yaml: &str) {
    nonbonded(
        yaml,
        Some(Medium::salt_water(298.15, Salt::SodiumChloride, 0.1)),
    );
}

/// Every pair of the model, without the medium's Coulomb term.
fn short_range_pairs(model: &str) -> Vec<PairPot> {
    let matrix = nonbonded(&topology_yaml(model, false), None);
    matrix.get_potentials().iter().cloned().collect()
}

fn range_of(model: &str) -> f64 {
    let ff: Box<dyn ForceField> = forcefield::from_name(model, HydrophobicScaling::NoScale)
        .unwrap()
        .unwrap();
    ff.range()
}

#[test]
fn generated_topologies_parse_in_faunus() {
    for model in ["calvados3", "kimhummer", "pasquier"] {
        for multi_bead in [false, true] {
            println!("{model} (multi_bead={multi_bead})");
            assert_faunus_parses_yaml(&topology_yaml(model, multi_bead));
        }
    }
}

/// Cross pairs exist only if both molecules share one topology.
#[test]
fn joint_topology_has_cross_pairs() {
    let other = CIF.replace("LEU", "ALA");
    let system = cgkitten::coarse_grain_jointly(
        [(CIF.as_bytes(), false), (other.as_bytes(), false)],
        &cgkitten::ChargeCalc::new().mc(0),
        "kimhummer",
        &cgkitten::SingleBead,
        0.02,
        HydrophobicScaling::NoScale,
        &[] as &[&str],
    )
    .unwrap();
    let yaml = &system.topology;
    assert!(
        yaml.contains("[LEU, ALA]") || yaml.contains("[ALA, LEU]"),
        "{yaml}"
    );
    assert!(!system.xyz[1].contains("LEU"), "{}", system.xyz[1]);
    assert!(system.xyz[1].contains("ALA"), "{}", system.xyz[1]);
    assert_faunus_parses_yaml(yaml);
}

/// C1: pairs are cut and shifted at the model's range, so consumers splining over
/// it lose nothing and see no jump.
#[test]
fn pairs_vanish_continuously_at_the_model_range() {
    for model in ["calvados3", "kimhummer", "pasquier"] {
        let range = range_of(model);
        for pair in short_range_pairs(model) {
            assert_eq!(pair.cutoff(), range, "{model}: {pair:?}");
            assert!(!pair.has_coulomb(), "{model}: {pair:?}");
            let u = pair.isotropic_twobody_energy((range - 1e-6).powi(2));
            assert!(u.abs() < 1e-9, "{model}: u = {u} at the range, {pair:?}");
            assert_eq!(pair.isotropic_twobody_energy(range.powi(2)), 0.0);
        }
    }
}
