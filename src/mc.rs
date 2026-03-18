//! Metropolis Monte Carlo titration using Yukawa (screened Coulomb) electrostatics.
//!
//! Titratable sites are randomly toggled between protonated and deprotonated
//! states. The acceptance criterion combines the screened Coulomb energy change
//! with the intrinsic protonation free energy:
//!
//! ΔU/kT = Δq · φ_i + ln(10) · (pKa − pH)   (for deprotonation)
//! ΔU/kT = Δq · φ_i + ln(10) · (pH − pKa)    (for protonation)
//!
//! where φ_i = λ_B · Σ_{j≠i} q_j · exp(−r_ij / λ_D) / r_ij

use std::f64::consts::PI;

use log::debug;
use rand::Rng;

use crate::Bead;
use crate::charge::Conditions;

const LN_10: f64 = std::f64::consts::LN_10;

// Physical constants (SI)
const ELEMENTARY_CHARGE: f64 = 1.602_176_634e-19; // C
const VACUUM_PERMITTIVITY: f64 = 8.854_187_812_8e-12; // F/m
const BOLTZMANN_CONSTANT: f64 = 1.380_649e-23; // J/K
const AVOGADRO: f64 = 6.022_140_76e23; // mol⁻¹
const WATER_DIELECTRIC: f64 = 78.4;

/// Bjerrum length in Ångströms.
///
/// λ_B = e² / (4π ε₀ ε_r k T)
pub fn bjerrum_length(temperature: f64) -> f64 {
    let lb_m = ELEMENTARY_CHARGE.powi(2)
        / (4.0 * PI * VACUUM_PERMITTIVITY * WATER_DIELECTRIC * BOLTZMANN_CONSTANT * temperature);
    lb_m * 1e10 // convert m → Å
}

/// Debye screening length in Ångströms.
///
/// λ_D = 1/κ where κ² = 2 N_A e² I · 1000 / (ε₀ ε_r k T)
pub fn debye_length(temperature: f64, ionic_strength: f64) -> f64 {
    if ionic_strength <= 0.0 {
        return f64::INFINITY;
    }
    // Factor 1000 converts mol/L to mol/m³
    let kappa_sq = 2.0 * AVOGADRO * ELEMENTARY_CHARGE.powi(2) * ionic_strength * 1000.0
        / (VACUUM_PERMITTIVITY * WATER_DIELECTRIC * BOLTZMANN_CONSTANT * temperature);
    (1.0 / kappa_sq.sqrt()) * 1e10 // convert m → Å
}

/// Convert a length to its inverse, treating infinity as zero.
/// When ionic strength is zero, λ_D → ∞ and 1/λ_D → 0, giving unscreened Coulomb.
fn inv_length(l: f64) -> f64 {
    if l.is_infinite() { 0.0 } else { 1.0 / l }
}

/// A titratable site tracked during MC.
pub(crate) struct TitratableSite {
    /// Index into the original beads array (for mapping results back).
    pub bead_index: usize,
    /// Current protonation state.
    pub is_protonated: bool,
    /// Charge change for deprotonation (q_deprot - q_prot).
    dq_deprot: f64,
    /// Chemical energy for deprotonation: ln(10) · (pKa - pH).
    du_chem_deprot: f64,
    /// Chemical energy for protonation: ln(10) · (pH - pKa).
    du_chem_prot: f64,
}

/// Result of an MC titration run.
pub(crate) struct TitrationResult {
    /// Average charge per bead (same indexing as input beads).
    pub charges: Vec<f64>,
    /// Ensemble-averaged net charge ⟨Z⟩.
    pub mean_charge: f64,
    /// Ensemble-averaged ⟨Z²⟩.
    pub mean_charge_sq: f64,
    /// Ensemble-averaged dipole moment magnitude ⟨μ⟩ in e·Å.
    pub mean_dipole: f64,
    /// Ensemble-averaged ⟨μ²⟩ in (e·Å)².
    pub mean_dipole_sq: f64,
}

/// Metropolis Monte Carlo titration engine.
///
/// Only titratable sites are stored. Precomputes the pairwise Yukawa kernel
/// between titratable sites, then performs MC moves toggling protonation states.
/// Metal ions contribute a constant background potential included in φ
/// but never updated during sweeps.
pub(crate) struct TitrationMC {
    /// Charges of titratable sites only (length = n_sites).
    charges: Vec<f64>,
    /// Positions of titratable sites only (length = n_sites).
    positions: Vec<[f64; 3]>,
    sites: Vec<TitratableSite>,
    /// Precomputed Yukawa kernel between titratable sites: n_sites × n_sites.
    kernel: Vec<f64>,
    /// Cached electrostatic potential at each site (kT/e), updated incrementally.
    phi: Vec<f64>,
    /// Precomputed background potential from fixed-charge beads, stored so
    /// recompute_phi (test helper) can correctly restore φ after parameter changes.
    #[cfg(test)]
    phi_fixed: Vec<f64>,
    lambda_b: f64,
    /// Running net charge (updated incrementally on accepted moves).
    z: f64,
    /// Running dipole moment components (updated incrementally).
    mx: f64,
    my: f64,
    mz: f64,
    /// Total number of beads in the original input (for result mapping).
    n_beads: usize,
}

impl TitrationMC {
    /// Create a new MC titration from coarse-grained beads.
    ///
    /// Only titratable sites are kept. The Yukawa kernel is computed between
    /// titratable sites only. Initial protonation states are set from
    /// Henderson-Hasselbalch.
    pub fn new(beads: &[Bead], conditions: &Conditions) -> Self {
        let lambda_b = bjerrum_length(conditions.temperature);
        let lambda_d = debye_length(conditions.temperature, conditions.ionic_strength);
        let inv_debye = inv_length(lambda_d);

        debug!("λ_B = {lambda_b:.2} Å, λ_D = {lambda_d:.2} Å");

        let rn = crate::res_name_map(beads);

        // Identify titratable sites and set initial discrete protonation states
        let mut sites = Vec::new();
        let mut charges = Vec::new();
        let mut positions = Vec::new();

        for (i, b) in beads.iter().enumerate() {
            if let Some(group) = crate::titratable_group_for_bead(b, &rn) {
                // Use HH fraction to set initial discrete state; MC will equilibrate from here
                let f_deprot = 1.0 / (1.0 + 10_f64.powf(group.pka - conditions.ph));
                let is_protonated = f_deprot < 0.5;
                charges.push(if is_protonated {
                    group.charge_protonated
                } else {
                    group.charge_deprotonated
                });
                positions.push([b.x, b.y, b.z]);
                sites.push(TitratableSite {
                    bead_index: i,
                    is_protonated,
                    dq_deprot: group.charge_deprotonated - group.charge_protonated,
                    du_chem_deprot: LN_10 * (group.pka - conditions.ph),
                    du_chem_prot: LN_10 * (conditions.ph - group.pka),
                });
            }
        }

        let ns = sites.len();
        debug!("{ns} titratable sites");

        // A nearby Zn²⁺ (or other metal ion) raises the local electrostatic potential and
        // shifts pKa values of titratable residues — a real physical effect that must be
        // included. Because metal ion charges never change during the simulation, their
        // Yukawa contribution to φ_i is a constant offset: compute it once here and add it
        // to φ. NTR/CTR already titrate and appear in titratable_indices; backbone beads
        // carry zero charge; so in practice only metal ions contribute.
        let titratable_indices: std::collections::HashSet<usize> =
            sites.iter().map(|s| s.bead_index).collect();
        let phi_fixed: Vec<f64> = positions
            .iter()
            .map(|&pos_i| {
                lambda_b
                    * beads
                        .iter()
                        .enumerate()
                        .filter(|(j, b)| !titratable_indices.contains(j) && b.charge != 0.0)
                        .map(|(_, b)| {
                            let dx = pos_i[0] - b.x;
                            let dy = pos_i[1] - b.y;
                            let dz = pos_i[2] - b.z;
                            let r = dz.mul_add(dz, dx.mul_add(dx, dy * dy)).sqrt();
                            if r > 0.0 {
                                b.charge * (-r * inv_debye).exp() / r
                            } else {
                                0.0
                            }
                        })
                        .sum::<f64>()
            })
            .collect();

        // ΔU for a protonation move only involves charge changes at titratable sites, so the
        // kernel only needs titratable-titratable pairs: O(n_sites²) instead of O(n_beads²).
        // Fixed-charge beads cancel out of ΔU and are handled by phi_fixed above.
        let mut kernel = vec![0.0; ns * ns];
        for i in 0..ns {
            for j in (i + 1)..ns {
                let dx = positions[i][0] - positions[j][0];
                let dy = positions[i][1] - positions[j][1];
                let dz = positions[i][2] - positions[j][2];
                let r = dz.mul_add(dz, dx.mul_add(dx, dy * dy)).sqrt();
                if r > 0.0 {
                    let k = (-r * inv_debye).exp() / r;
                    kernel[i * ns + j] = k;
                    kernel[j * ns + i] = k;
                }
            }
        }

        // Compute initial cached potentials from titratable-titratable interactions
        // (diagonal is zero, no branch needed), then add the fixed-charge background
        // so φ_i correctly reflects the full electrostatic environment from the start.
        let phi = (0..ns)
            .map(|i| {
                let row = &kernel[i * ns..(i + 1) * ns];
                let phi_tit = lambda_b
                    * row
                        .iter()
                        .zip(charges.iter())
                        .map(|(&k, &q)| q * k)
                        .sum::<f64>();
                phi_tit + phi_fixed[i]
            })
            .collect();

        // Initial running multipole from current charges
        let mut z = 0.0;
        let (mut mx, mut my, mut mz_val) = (0.0, 0.0, 0.0);
        for (i, &q) in charges.iter().enumerate() {
            z += q;
            let [x, y, zz] = positions[i];
            mx += q * x;
            my += q * y;
            mz_val += q * zz;
        }

        Self {
            charges,
            positions,
            sites,
            kernel,
            #[cfg(test)]
            phi_fixed,
            phi,
            lambda_b,
            z,
            mx,
            my,
            mz: mz_val,
            n_beads: beads.len(),
        }
    }

    /// Override the Bjerrum length (e.g., set to 0 to disable electrostatics).
    #[cfg(test)]
    pub fn set_bjerrum_length(&mut self, lambda_b: f64) {
        self.lambda_b = lambda_b;
        self.recompute_phi();
    }

    /// Recompute all cached potentials from scratch.
    #[cfg(test)]
    fn recompute_phi(&mut self) {
        let ns = self.sites.len();
        for i in 0..ns {
            let row = &self.kernel[i * ns..(i + 1) * ns];
            self.phi[i] = self.lambda_b
                * row
                    .iter()
                    .zip(self.charges.iter())
                    .map(|(&k, &q)| q * k)
                    .sum::<f64>()
                + self.phi_fixed[i];
        }
    }

    /// Number of titratable sites.
    pub fn num_sites(&self) -> usize {
        self.sites.len()
    }

    /// Electrostatic energy of a titratable site due to all other sites (in kT).
    #[cfg(test)]
    pub fn site_energy(&self, site_idx: usize) -> f64 {
        self.charges[site_idx] * self.phi[site_idx]
    }

    /// Perform one MC sweep (N random trial moves, where N = number of sites).
    ///
    /// Returns the number of accepted moves.
    pub fn sweep<R: Rng>(&mut self, rng: &mut R) -> usize {
        let n_sites = self.sites.len();
        if n_sites == 0 {
            return 0;
        }

        let mut accepted = 0;

        for _ in 0..n_sites {
            let site_idx = rng.r#gen_range(0..n_sites);
            let site = &self.sites[site_idx];

            let (dq, du_chem) = if site.is_protonated {
                (site.dq_deprot, site.du_chem_deprot)
            } else {
                (-site.dq_deprot, site.du_chem_prot)
            };

            let du = dq.mul_add(self.phi[site_idx], du_chem);

            if du < 0.0 || rng.r#gen::<f64>() < (-du).exp() {
                self.charges[site_idx] += dq;
                self.sites[site_idx].is_protonated = !self.sites[site_idx].is_protonated;
                accepted += 1;

                // Incrementally update running multipole
                self.z += dq;
                let [x, y, zz] = self.positions[site_idx];
                self.mx += dq * x;
                self.my += dq * y;
                self.mz += dq * zz;

                // Update all φ_j incrementally: O(n_sites) instead of O(n_sites²) full recompute.
                // Diagonal kernel[i][i] = 0, so the self-potential stays correct without branching.
                let dq_lb = dq * self.lambda_b;
                let row = &self.kernel[site_idx * n_sites..(site_idx + 1) * n_sites];
                for (phi, &k) in self.phi.iter_mut().zip(row.iter()) {
                    *phi += dq_lb * k;
                }
            }
        }

        accepted
    }

    /// Run MC titration and return ensemble-averaged results.
    ///
    /// Performs `sweeps` sweeps (each sweep = one trial move per titratable site),
    /// collecting per-site charge averages, ⟨Z⟩, and ⟨μ⟩ from every sweep.
    pub fn run<R: Rng>(&mut self, sweeps: usize, rng: &mut R) -> TitrationResult {
        let ns = self.sites.len();
        let mut site_sum = vec![0.0; ns];
        let mut charge_sum = 0.0;
        let mut charge_sq_sum = 0.0;
        let mut dipole_sum = 0.0;
        let mut dipole_sq_sum = 0.0;

        for _ in 0..sweeps {
            self.sweep(rng);

            // Sample running multipole (updated incrementally in sweep)
            for (i, &q) in self.charges.iter().enumerate() {
                site_sum[i] += q;
            }
            charge_sum += self.z;
            charge_sq_sum += self.z * self.z;
            let mu_sq = self
                .mz
                .mul_add(self.mz, self.mx.mul_add(self.mx, self.my * self.my));
            dipole_sum += mu_sq.sqrt();
            dipole_sq_sum += mu_sq;
        }

        // Map site-averaged charges back to full bead array
        let mut charges = vec![0.0; self.n_beads];
        for (si, site) in self.sites.iter().enumerate() {
            charges[site.bead_index] = site_sum[si] / sweeps as f64;
        }

        let s = sweeps as f64;
        TitrationResult {
            charges,
            mean_charge: charge_sum / s,
            mean_charge_sq: charge_sq_sum / s,
            mean_dipole: dipole_sum / s,
            mean_dipole_sq: dipole_sq_sum / s,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BeadType, coarse_grain};
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn bjerrum_length_water_298k() {
        let lb = bjerrum_length(298.15);
        assert!((lb - 7.12).abs() < 0.1, "λ_B = {lb} Å");
    }

    #[test]
    fn debye_length_0_1m() {
        let ld = debye_length(298.15, 0.1);
        assert!((ld - 9.6).abs() < 0.5, "λ_D = {ld} Å");
    }

    #[test]
    fn debye_length_zero_ionic_strength() {
        let ld = debye_length(298.15, 0.0);
        assert!(ld.is_infinite());
    }

    #[test]
    fn isolated_asp_at_ph7() {
        let cif = r#"
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
ATOM 1 N ASP A 1 0.0 0.0 0.0 N . 1
ATOM 2 CA ASP A 1 1.5 0.0 0.0 C . 1
ATOM 3 C ASP A 1 2.4 -1.1 0.0 C . 1
ATOM 4 O ASP A 1 2.2 -2.3 0.0 O . 1
ATOM 5 OD1 ASP A 1 100.0 100.0 100.0 O . 1
ATOM 6 OD2 ASP A 1 100.0 100.0 102.0 O . 1
"#;
        let beads = coarse_grain(cif.as_bytes());
        assert_eq!(beads.len(), 4);

        let conditions = Conditions::default();
        let mut mc = TitrationMC::new(&beads, &conditions);
        assert_eq!(mc.num_sites(), 3);

        let mut rng = StdRng::seed_from_u64(42);
        let result = mc.run(5000, &mut rng);

        let sc = beads
            .iter()
            .position(|b| b.bead_type == BeadType::Virtual)
            .unwrap();
        assert!(
            (result.charges[sc] + 1.0).abs() < 0.05,
            "ASP avg charge at pH 7: {}",
            result.charges[sc]
        );
    }

    #[test]
    fn isolated_asp_at_pka() {
        let cif = r#"
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
ATOM 1 N ASP A 1 0.0 0.0 0.0 N . 1
ATOM 2 CA ASP A 1 1.5 0.0 0.0 C . 1
ATOM 3 C ASP A 1 2.4 -1.1 0.0 C . 1
ATOM 4 O ASP A 1 2.2 -2.3 0.0 O . 1
ATOM 5 OD1 ASP A 1 100.0 100.0 100.0 O . 1
ATOM 6 OD2 ASP A 1 100.0 100.0 102.0 O . 1
"#;
        let beads = coarse_grain(cif.as_bytes());

        let conditions = Conditions {
            ph: 3.65,
            ..Default::default()
        };
        let mut mc = TitrationMC::new(&beads, &conditions);
        let mut rng = StdRng::seed_from_u64(123);
        let result = mc.run(20000, &mut rng);

        let sc = beads
            .iter()
            .position(|b| b.bead_type == BeadType::Virtual)
            .unwrap();
        assert!(
            (result.charges[sc] + 0.5).abs() < 0.05,
            "ASP avg charge at pKa: {}",
            result.charges[sc]
        );
    }

    #[test]
    fn mc_without_electrostatics_matches_hh() {
        use crate::charge::{ChargeCalculator, HendersonHasselbalch};
        use crate::residue::TITRATABLE_GROUPS;

        let ph = 7.0;
        let conditions = Conditions {
            ph,
            ..Default::default()
        };
        let hh = HendersonHasselbalch;

        let mut beads = Vec::new();
        for (i, group) in TITRATABLE_GROUPS.iter().enumerate() {
            let x = (i as f64) * 1000.0;
            beads.push(Bead {
                x,
                y: 0.0,
                z: 0.0,
                charge: 0.0,
                res_name: group.res_name.to_string(),
                chain_id: "A".into(),
                res_seq: (i + 1) as i32,
                mass: 0.0,
                bead_type: BeadType::Residue,
            });
            beads.push(Bead {
                x,
                y: 10.0,
                z: 0.0,
                charge: hh.charge(group, &conditions),
                res_name: group.element.to_string(),
                chain_id: "A".into(),
                res_seq: (i + 1) as i32,
                mass: 0.0,
                bead_type: BeadType::Virtual,
            });
        }

        let mut mc = TitrationMC::new(&beads, &conditions);
        mc.set_bjerrum_length(0.0);

        let mut rng = StdRng::seed_from_u64(42);
        let result = mc.run(50000, &mut rng);

        for (i, group) in TITRATABLE_GROUPS.iter().enumerate() {
            let hh_charge = hh.charge(group, &conditions);
            let mc_charge = result.charges[2 * i + 1];
            assert!(
                (mc_charge - hh_charge).abs() < 0.03,
                "{}: MC={mc_charge:.4}, HH={hh_charge:.4}",
                group.res_name,
            );
        }
    }

    #[test]
    fn site_energy_sign() {
        let beads = vec![
            Bead {
                x: 0.0,
                y: 0.0,
                z: 0.0,
                charge: 0.0,
                res_name: "ASP".into(),
                chain_id: "A".into(),
                res_seq: 1,
                mass: 0.0,
                bead_type: BeadType::Residue,
            },
            Bead {
                x: 5.0,
                y: 0.0,
                z: 0.0,
                charge: -1.0,
                res_name: "O".into(),
                chain_id: "A".into(),
                res_seq: 1,
                mass: 0.0,
                bead_type: BeadType::Virtual,
            },
            Bead {
                x: 10.0,
                y: 0.0,
                z: 0.0,
                charge: 0.0,
                res_name: "GLU".into(),
                chain_id: "A".into(),
                res_seq: 2,
                mass: 0.0,
                bead_type: BeadType::Residue,
            },
            Bead {
                x: 12.0,
                y: 0.0,
                z: 0.0,
                charge: -1.0,
                res_name: "O".into(),
                chain_id: "A".into(),
                res_seq: 2,
                mass: 0.0,
                bead_type: BeadType::Virtual,
            },
        ];

        let conditions = Conditions::default();
        let mc = TitrationMC::new(&beads, &conditions);

        let e0 = mc.site_energy(0);
        let e1 = mc.site_energy(1);
        assert!(e0 > 0.0, "repulsive energy site 0: {e0}");
        assert!(e1 > 0.0, "repulsive energy site 1: {e1}");
    }
}
