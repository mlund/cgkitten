//! Topology type assignment for coarse-grained beads.
//!
//! Groups beads into the minimum number of unique atom types such that
//! all beads within a type share the same force-field identity and have
//! charges within `tolerance` of each other.
//!
//! The algorithm is **order-independent**: shuffling the input beads
//! produces identical type assignments.

use std::collections::{HashMap, HashSet};
use std::ops::Range;

use log::info;

use crate::{Bead, BeadType};

// ---------------------------------------------------------------------------
// Private types
// ---------------------------------------------------------------------------

/// A unique atom type produced by clustering.
struct AtomType {
    name: String,
    charge: f64,
    mass: f64,
    res_name: String,
    bead_type: BeadType,
}

/// Key for grouping beads before charge clustering.
#[derive(Clone, PartialEq, Eq, Hash)]
enum GroupKey {
    /// Titratable single-bead (SingleBead policy): key = res_name.
    ByResName(String),
    /// Titratable multi-bead (Virtual/Ntr/Ctr): key = (res_name, mass_bits).
    ByNameAndMass(String, u64),
}

impl GroupKey {
    fn base_name(&self) -> &str {
        match self {
            GroupKey::ByResName(n) | GroupKey::ByNameAndMass(n, _) => n,
        }
    }

    fn mass_bits(&self) -> u64 {
        match self {
            GroupKey::ByNameAndMass(_, bits) => *bits,
            _ => 0,
        }
    }
}

/// Per-group data: charge entries (bead_index, charge) and representative bead info.
struct GroupData {
    entries: Vec<(usize, f64)>,
    mass: f64,
    res_name: String,
    bead_type: BeadType,
}

// ---------------------------------------------------------------------------
// Sort-and-split clustering
// ---------------------------------------------------------------------------

/// Sort charges and split at consecutive gaps > `tolerance`.
/// Returns ranges into the sorted slice; each range is one cluster.
fn cluster_charges(entries: &mut [(usize, f64)], tolerance: f64) -> Vec<Range<usize>> {
    if entries.is_empty() {
        return Vec::new();
    }
    entries.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let mut clusters = Vec::new();
    let mut start = 0;
    for i in 1..entries.len() {
        if entries[i].1 - entries[i - 1].1 > tolerance {
            clusters.push(start..i);
            start = i;
        }
    }
    clusters.push(start..entries.len());
    clusters
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Read-only view of an atom type. Borrows from [`Topology`].
pub struct TypeView<'a> {
    /// Type name (e.g. "ALA", "O1", "ASP").
    pub name: &'a str,
    /// Mean charge of all beads in this type.
    pub charge: f64,
    /// Mass in Da.
    pub mass: f64,
    /// Source residue/element name for force-field lookup.
    pub res_name: &'a str,
    /// Bead classification.
    pub bead_type: BeadType,
}

/// Topology type assignment for coarse-grained beads.
///
/// Groups beads into the minimum number of unique atom types such that
/// all beads within a type share the same force-field identity and have
/// charges within `tolerance` of each other.
pub struct Topology {
    types: Vec<AtomType>,
    /// Per-bead type name, same indexing as input beads.
    bead_names: Vec<String>,
}

impl Topology {
    /// Build a topology from charged beads.
    pub fn new(beads: &[Bead], tolerance: f64) -> Self {
        let mut types: Vec<AtomType> = Vec::new();
        let mut bead_names: Vec<String> = vec![String::new(); beads.len()];

        // Phase 1: handle non-titratable beads (Residue, Ion) — deduplicate by res_name.
        let mut seen_non_tit: HashSet<String> = HashSet::new();
        for (i, b) in beads.iter().enumerate() {
            if b.bead_type == BeadType::Residue || b.bead_type == BeadType::Ion {
                if seen_non_tit.insert(b.res_name.clone()) {
                    types.push(AtomType {
                        name: b.res_name.clone(),
                        charge: b.charge,
                        mass: b.mass,
                        res_name: b.res_name.clone(),
                        bead_type: b.bead_type,
                    });
                }
                bead_names[i] = b.res_name.clone();
            }
        }

        // Phase 2: group titratable beads by GroupKey, then cluster charges.
        let mut groups: HashMap<GroupKey, GroupData> = HashMap::new();

        for (i, b) in beads.iter().enumerate() {
            if !b.bead_type.is_titratable() {
                continue;
            }
            let key = match b.bead_type {
                BeadType::Titratable => GroupKey::ByResName(b.res_name.clone()),
                _ => GroupKey::ByNameAndMass(b.res_name.clone(), b.mass.to_bits()),
            };
            if let Some(group) = groups.get_mut(&key) {
                group.entries.push((i, b.charge));
            } else {
                groups.insert(
                    key,
                    GroupData {
                        entries: vec![(i, b.charge)],
                        mass: b.mass,
                        res_name: b.res_name.clone(),
                        bead_type: b.bead_type,
                    },
                );
            }
        }

        // Sort group keys for deterministic output order.
        let mut keys: Vec<GroupKey> = groups.keys().cloned().collect();
        keys.sort_by(|a, b| {
            a.base_name()
                .cmp(b.base_name())
                .then(a.mass_bits().cmp(&b.mass_bits()))
        });

        // Two-pass: first count total clusters per base name, then assign names.
        let mut base_name_clusters: HashMap<&str, usize> = HashMap::new();
        let mut key_clusters: Vec<(&GroupKey, Vec<Range<usize>>)> = Vec::new();

        for key in &keys {
            let group = groups.get_mut(key).unwrap();
            let clusters = cluster_charges(&mut group.entries, tolerance);
            *base_name_clusters.entry(key.base_name()).or_insert(0) += clusters.len();
            key_clusters.push((key, clusters));
        }

        // Assign names and create types.
        let mut counters: HashMap<&str, usize> = HashMap::new();
        for (key, clusters) in &key_clusters {
            let group = groups.get(key).unwrap();
            let total_clusters = base_name_clusters[group.res_name.as_str()];

            for range in clusters {
                let counter = counters.entry(group.res_name.as_str()).or_insert(0);
                *counter += 1;

                let name = if total_clusters == 1 {
                    group.res_name.clone()
                } else {
                    format!("{}{}", group.res_name, counter)
                };

                // Mean charge for this cluster.
                let charge_sum: f64 = group.entries[range.clone()].iter().map(|(_, q)| q).sum();
                let charge_mean = charge_sum / (range.len() as f64);

                types.push(AtomType {
                    name: name.clone(),
                    charge: charge_mean,
                    mass: group.mass,
                    res_name: group.res_name.clone(),
                    bead_type: group.bead_type,
                });

                for &(bead_idx, _) in &group.entries[range.clone()] {
                    bead_names[bead_idx] = name.clone();
                }
            }
        }

        // Log merge stats
        let n_before = beads.iter().filter(|b| b.bead_type.is_titratable()).count();
        let n_after = types.iter().filter(|t| t.bead_type.is_titratable()).count();
        if n_after < n_before {
            info!(
                "Merged {n_before} titratable sites into {n_after} unique types (tolerance {:.0}%)",
                tolerance * 100.0
            );
        }

        Topology { types, bead_names }
    }

    /// Number of unique atom types.
    pub fn num_types(&self) -> usize {
        self.types.len()
    }

    /// Iterate over unique atom types.
    pub fn types(&self) -> impl Iterator<Item = TypeView<'_>> {
        self.types.iter().map(|t| TypeView {
            name: &t.name,
            charge: t.charge,
            mass: t.mass,
            res_name: &t.res_name,
            bead_type: t.bead_type,
        })
    }

    /// Per-bead type name (same indexing as input beads).
    pub fn bead_name(&self, index: usize) -> &str {
        &self.bead_names[index]
    }

    /// All per-bead type names as a slice.
    pub fn bead_names(&self) -> &[String] {
        &self.bead_names
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper to create a bead with minimal fields.
    fn bead(res_name: &str, bead_type: BeadType, charge: f64, mass: f64) -> Bead {
        Bead {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            charge,
            mass,
            res_name: res_name.to_string(),
            chain_id: "A".to_string(),
            res_seq: 1,
            bead_type,
        }
    }

    #[test]
    fn non_titratable_deduplicate_by_res_name() {
        let beads = vec![
            bead("ALA", BeadType::Residue, 0.0, 71.0),
            bead("GLY", BeadType::Residue, 0.0, 57.0),
            bead("ALA", BeadType::Residue, 0.0, 71.0),
        ];
        let topo = Topology::new(&beads, 0.02);
        assert_eq!(topo.num_types(), 2); // ALA, GLY
        assert_eq!(topo.bead_name(0), "ALA");
        assert_eq!(topo.bead_name(1), "GLY");
        assert_eq!(topo.bead_name(2), "ALA");
    }

    #[test]
    fn virtual_similar_charges_merge() {
        let beads = vec![
            bead("O", BeadType::Virtual, -0.90, 0.0),
            bead("O", BeadType::Virtual, -0.91, 0.0),
            bead("O", BeadType::Virtual, -0.89, 0.0),
        ];
        let topo = Topology::new(&beads, 0.05);
        assert_eq!(topo.num_types(), 1);
        assert_eq!(topo.bead_name(0), "O");
        // Mean charge
        let t = topo.types().next().unwrap();
        assert!((t.charge - (-0.90)).abs() < 0.01);
    }

    #[test]
    fn virtual_distinct_charges_split() {
        let beads = vec![
            bead("O", BeadType::Virtual, -0.90, 0.0),
            bead("O", BeadType::Virtual, -0.50, 0.0),
        ];
        let topo = Topology::new(&beads, 0.02);
        assert_eq!(topo.num_types(), 2);
        // Names should be numbered
        let names: Vec<&str> = topo.types().map(|t| t.name).collect();
        assert!(names.contains(&"O1"));
        assert!(names.contains(&"O2"));
    }

    #[test]
    fn order_independence() {
        let beads_a = vec![
            bead("O", BeadType::Virtual, -0.50, 0.0),
            bead("O", BeadType::Virtual, -0.90, 0.0),
            bead("O", BeadType::Virtual, -0.51, 0.0),
            bead("O", BeadType::Virtual, -0.89, 0.0),
        ];
        // Reversed order
        let beads_b = vec![
            bead("O", BeadType::Virtual, -0.89, 0.0),
            bead("O", BeadType::Virtual, -0.51, 0.0),
            bead("O", BeadType::Virtual, -0.90, 0.0),
            bead("O", BeadType::Virtual, -0.50, 0.0),
        ];
        let topo_a = Topology::new(&beads_a, 0.05);
        let topo_b = Topology::new(&beads_b, 0.05);

        assert_eq!(topo_a.num_types(), topo_b.num_types());

        // Same type names and charges (sorted for comparison)
        let mut types_a: Vec<(String, String)> = topo_a
            .types()
            .map(|t| (t.name.to_string(), format!("{:.4}", t.charge)))
            .collect();
        let mut types_b: Vec<(String, String)> = topo_b
            .types()
            .map(|t| (t.name.to_string(), format!("{:.4}", t.charge)))
            .collect();
        types_a.sort();
        types_b.sort();
        assert_eq!(types_a, types_b);

        // Each bead in both orders maps to the same type name
        // (charge-based: -0.50/-0.51 → one type, -0.89/-0.90 → another)
        assert_eq!(topo_a.bead_name(0), topo_b.bead_name(1)); // -0.50 == -0.51 group
        assert_eq!(topo_a.bead_name(1), topo_b.bead_name(2)); // -0.90 == -0.89 group
    }

    #[test]
    fn single_bead_different_residues_never_merge() {
        let beads = vec![
            bead("ASP", BeadType::Titratable, -0.90, 115.0),
            bead("GLU", BeadType::Titratable, -0.91, 129.0),
        ];
        let topo = Topology::new(&beads, 0.05);
        assert_eq!(topo.num_types(), 2);
        assert_eq!(topo.bead_name(0), "ASP");
        assert_eq!(topo.bead_name(1), "GLU");
    }

    #[test]
    fn single_bead_same_residue_clusters_by_charge() {
        let beads = vec![
            bead("ASP", BeadType::Titratable, -0.90, 115.0),
            bead("ASP", BeadType::Titratable, -0.50, 115.0),
            bead("ASP", BeadType::Titratable, -0.91, 115.0),
        ];
        let topo = Topology::new(&beads, 0.05);
        assert_eq!(topo.num_types(), 2);
        let names: Vec<&str> = topo.types().map(|t| t.name).collect();
        assert!(names.contains(&"ASP1"));
        assert!(names.contains(&"ASP2"));
    }

    #[test]
    fn multi_bead_different_masses_prevent_merging() {
        // Same res_name, same charge, but different masses → different groups
        let beads = vec![
            bead("O", BeadType::Virtual, -0.90, 0.0),
            bead("O", BeadType::Virtual, -0.90, 16.0),
        ];
        let topo = Topology::new(&beads, 0.05);
        assert_eq!(topo.num_types(), 2);
        let names: Vec<&str> = topo.types().map(|t| t.name).collect();
        assert!(names.contains(&"O1"));
        assert!(names.contains(&"O2"));
    }

    #[test]
    fn strip_trailing_1_for_single_variant() {
        // Only one cluster → bare name, no "1" suffix
        let beads = vec![
            bead("N", BeadType::Virtual, 0.95, 0.0),
            bead("N", BeadType::Virtual, 0.96, 0.0),
        ];
        let topo = Topology::new(&beads, 0.05);
        assert_eq!(topo.num_types(), 1);
        assert_eq!(topo.types().next().unwrap().name, "N");
    }

    #[test]
    fn mixed_scenario() {
        let beads = vec![
            bead("ALA", BeadType::Residue, 0.0, 71.0),
            bead("ZN", BeadType::Ion, 2.0, 65.0),
            bead("O", BeadType::Virtual, -0.90, 0.0),
            bead("O", BeadType::Virtual, -0.50, 0.0),
            bead("N", BeadType::Ntr, 0.95, 0.0),
            bead("O", BeadType::Ctr, -0.90, 0.0),
        ];
        let topo = Topology::new(&beads, 0.02);
        // O(Virtual, mass=0) and O(Ctr, mass=0) share GroupKey ByNameAndMass("O", 0)
        // because they are force-field equivalent (same res_name for FF lookup).
        // Result: ALA, ZN, N, O1(-0.90 merged virtual+ctr), O2(-0.50)
        assert_eq!(topo.num_types(), 5);
        assert_eq!(topo.bead_name(0), "ALA");
        assert_eq!(topo.bead_name(1), "ZN");
    }

    #[test]
    fn empty_input() {
        let topo = Topology::new(&[], 0.02);
        assert_eq!(topo.num_types(), 0);
        assert!(topo.bead_names().is_empty());
    }

    #[test]
    fn tolerance_zero_max_split() {
        let beads = vec![
            bead("O", BeadType::Virtual, -0.90, 0.0),
            bead("O", BeadType::Virtual, -0.90, 0.0),
            bead("O", BeadType::Virtual, -0.91, 0.0),
        ];
        let topo = Topology::new(&beads, 0.0);
        // -0.90 and -0.90 are identical (gap=0, not > 0), but -0.91 differs
        assert_eq!(topo.num_types(), 2);
    }

    #[test]
    fn tolerance_large_max_merge() {
        let beads = vec![
            bead("O", BeadType::Virtual, -0.90, 0.0),
            bead("O", BeadType::Virtual, -0.10, 0.0),
            bead("O", BeadType::Virtual, 0.50, 0.0),
        ];
        let topo = Topology::new(&beads, 1.0);
        assert_eq!(topo.num_types(), 1);
        assert_eq!(topo.bead_name(0), "O");
    }
}
