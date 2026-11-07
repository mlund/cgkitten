//! Screened Coulomb (Debye–Hückel) potential of bead charges on a surface, the same
//! model the titration uses, so colours and titration agree.

use crate::mc::yukawa;
use crate::surface::Mesh;
use std::collections::HashMap;

/// Water probe radius (Å). The potential is read where a solvent molecule's centre
/// can be, rather than on the surface, where nearby charges spike.
const PROBE_RADIUS: f64 = 1.4;

/// Screened charges beyond this many Debye lengths add < e⁻³ of a near one, invisible
/// in a colour map; skipping them makes large assemblies tractable.
const CUTOFF_DEBYE_LENGTHS: f64 = 3.0;

/// A bead's charge and its titration fluctuation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointCharge {
    /// Å.
    pub position: [f64; 3],
    /// Mean charge, e.
    pub charge: f64,
    /// Charge variance, e²; zero for fixed charges.
    pub variance: f64,
}

/// Mean potential φ and its standard deviation σ_φ, sites independent.
pub struct SurfacePotential {
    bjerrum: f64,
    inv_debye: f64,
    /// Only charges that can contribute.
    sources: Vec<PointCharge>,
    /// With screening, sources by cell of edge = cutoff; without, every charge matters.
    grid: Option<CellGrid>,
}

impl SurfacePotential {
    /// `bjerrum` and `debye` in Å. An infinite `debye` means no salt: nothing is then
    /// screened, so every point sums every charge.
    pub fn new(sources: &[PointCharge], bjerrum: f64, debye: f64) -> Self {
        let inv_debye = crate::mc::inv_length(debye);
        let sources: Vec<_> = sources
            .iter()
            .filter(|s| s.charge != 0.0 || s.variance > 0.0)
            .copied()
            .collect();
        let grid = (inv_debye > 0.0).then(|| {
            CellGrid::new(
                sources.iter().map(|s| s.position),
                CUTOFF_DEBYE_LENGTHS / inv_debye,
            )
        });
        Self {
            bjerrum,
            inv_debye,
            sources,
            grid,
        }
    }

    /// (φ, σ_φ) in kT/e one water probe radius outside each vertex.
    pub fn on_mesh(&self, mesh: &Mesh) -> Vec<(f64, f64)> {
        let probes: Vec<_> = mesh
            .positions
            .iter()
            .zip(&mesh.normals)
            .map(|(p, n)| [0, 1, 2].map(|k| n[k].mul_add(PROBE_RADIUS, p[k])))
            .collect();
        self.at(&probes)
    }

    /// (φ, σ_φ) in kT/e at each point.
    pub fn at(&self, points: &[[f64; 3]]) -> Vec<(f64, f64)> {
        points.iter().map(|&p| self.at_point(p)).collect()
    }

    fn at_point(&self, p: [f64; 3]) -> (f64, f64) {
        match &self.grid {
            None => self.sum(p, self.sources.iter()),
            Some(grid) => {
                let cutoff_sq = grid.size * grid.size;
                let near = grid
                    .near(p)
                    .map(|i| &self.sources[i])
                    .filter(|s| distance_sq(p, s.position) <= cutoff_sq);
                self.sum(p, near)
            }
        }
    }

    fn sum<'a>(&self, p: [f64; 3], sources: impl Iterator<Item = &'a PointCharge>) -> (f64, f64) {
        let (mut mean, mut var) = (0.0, 0.0);
        for s in sources {
            let k = self.bjerrum * yukawa(distance_sq(p, s.position).sqrt(), self.inv_debye);
            mean += s.charge * k;
            var += s.variance * k * k;
        }
        (mean, var.sqrt())
    }
}

/// Uniform cell grid for neighbour queries within one fixed cutoff: linear in the
/// number of points, and simpler than a tree for a single cutoff.
struct CellGrid {
    size: f64,
    cells: HashMap<[i64; 3], Vec<usize>>,
}

impl CellGrid {
    /// Cells of edge `size`, which must be at least the largest query distance.
    fn new(points: impl IntoIterator<Item = [f64; 3]>, size: f64) -> Self {
        let mut cells: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
        for (i, p) in points.into_iter().enumerate() {
            cells.entry(cell(p, size)).or_default().push(i);
        }
        Self { size, cells }
    }

    /// Indices in the 27 cells around `p`: a superset of those within `size`.
    fn near(&self, p: [f64; 3]) -> impl Iterator<Item = usize> + '_ {
        let [x, y, z] = cell(p, self.size);
        (-1..=1)
            .flat_map(move |dx| {
                (-1..=1).flat_map(move |dy| (-1..=1).map(move |dz| [x + dx, y + dy, z + dz]))
            })
            .filter_map(|c| self.cells.get(&c))
            .flatten()
            .copied()
    }
}

fn cell(p: [f64; 3], size: f64) -> [i64; 3] {
    p.map(|c| (c / size).floor() as i64)
}

fn distance_sq(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).powi(2)).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    const BJERRUM: f64 = 7.0;

    fn charge(position: [f64; 3], charge: f64, variance: f64) -> PointCharge {
        PointCharge {
            position,
            charge,
            variance,
        }
    }

    #[test]
    fn one_charge_is_screened_coulomb() {
        let q = [charge([1.0, 1.0, 1.0], -2.0, 0.0)];
        let p = [[1.0, 1.0, 6.0]]; // r = 5 Å
        let (phi, sigma) = SurfacePotential::new(&q, BJERRUM, 10.0).at(&p)[0];
        // −2 · 7 · e^{−0.5} / 5
        assert!((phi - (-1.698_286)).abs() < 1e-6, "φ = {phi}");
        assert_eq!(sigma, 0.0);
        // No salt: −2 · 7 / 5.
        let (phi, _) = SurfacePotential::new(&q, BJERRUM, f64::INFINITY).at(&p)[0];
        assert!((phi - (-2.8)).abs() < 1e-12, "φ = {phi}");
    }

    #[test]
    fn fluctuating_bead_spreads_the_potential() {
        let q = [charge([0.0; 3], 0.0, 0.25)];
        let (phi, sigma) = SurfacePotential::new(&q, BJERRUM, 10.0).at(&[[5.0, 0.0, 0.0]])[0];
        assert_eq!(phi, 0.0);
        // √0.25 · 7 · e^{−0.5} / 5
        assert!((sigma - 0.424_571).abs() < 1e-6, "σ_φ = {sigma}");
    }

    /// Screened charges beyond 3 λ_D are skipped: invisible in a colour map, and
    /// skipping them keeps large assemblies tractable.
    #[test]
    fn charges_beyond_three_debye_lengths_are_skipped() {
        let near = charge([0.0, 0.0, 5.0], 1.0, 0.0);
        let far = charge([0.0, 0.0, 35.0], 1.0, 0.0); // 3.5 λ_D away
        let p = [[0.0; 3]];
        let alone = SurfacePotential::new(&[near], BJERRUM, 10.0).at(&p)[0];
        let both = SurfacePotential::new(&[near, far], BJERRUM, 10.0).at(&p)[0];
        assert_eq!(both, alone);
    }

    /// The cutoff error is bounded by what the skipped charges could add at 3 λ_D.
    #[test]
    fn cutoff_matches_full_sum_within_e_minus_3() {
        const DEBYE: f64 = 4.0;
        // Deterministic cloud spanning several Debye lengths.
        let sources: Vec<_> = (0..300)
            .map(|i| {
                let t = f64::from(i);
                let position = [
                    (t * 0.37).sin() * 20.0,
                    (t * 0.91).cos() * 20.0,
                    (t * 0.13).sin() * 20.0,
                ];
                charge(position, if i % 3 == 0 { -1.0 } else { 0.5 }, 0.0)
            })
            .collect();
        let points = [[0.0; 3], [7.0, -3.0, 2.0], [-15.0, 10.0, 4.0]];
        let cut = SurfacePotential::new(&sources, BJERRUM, DEBYE).at(&points);
        let total_charge: f64 = sources.iter().map(|s| s.charge.abs()).sum();
        let bound = total_charge * BJERRUM * (-3.0_f64).exp() / (3.0 * DEBYE);
        for (p, (phi, _)) in points.iter().zip(cut) {
            let full: f64 = sources
                .iter()
                .map(|s| {
                    let r = distance_sq(*p, s.position).sqrt();
                    s.charge * BJERRUM * (-r / DEBYE).exp() / r
                })
                .sum();
            assert!((phi - full).abs() <= bound, "φ = {phi}, full = {full}");
            assert!(phi != full || bound == 0.0, "nothing was cut at {p:?}");
        }
    }
}
