// Mesher vendored from bemtzmann (https://github.com/mlund/bemtzmann, MIT OR Apache-2.0),
// src/geometry/molecular_surface at e3e3e0fa61cee83c67809133944ae97f05b6ca92, via purrette.
// Changes across the vendored files: errors as String, serial density sampling (no
// rayon on wasm), upstream unit tests and the test-only gradient dropped.

//! Gaussian molecular surface as a triangle mesh with vertex normals.
//!
//! Density ρ(r) = Σᵢ exp(d·(1 − |r − rᵢ|²/Rᵢ²)) is sampled on a grid (serially: no
//! threads on wasm), polygonised at ρ = 1 by marching cubes and Taubin-smoothed.
//! With d = 1 a lone sphere's surface lies on its radius.

mod density;
mod grid;
mod marching_cubes;
mod smoothing;

use grid::Grid;

type Result<T> = std::result::Result<T, String>;

/// Closed, outward-wound triangle mesh.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Mesh {
    /// Å.
    pub positions: Vec<[f64; 3]>,
    /// Outward unit normals, one per position.
    pub normals: Vec<[f64; 3]>,
    /// Counter-clockwise seen from outside.
    pub triangles: Vec<[u32; 3]>,
}

/// Surface around spheres; spacing coarsens with size so large assemblies still
/// mesh in about a second.
pub fn gaussian_surface(centers: &[[f64; 3]], radii: &[f64]) -> Result<Mesh> {
    let spacing = match centers.len() {
        0..5_000 => 1.0,
        5_000..50_000 => 1.5,
        _ => 2.5,
    };
    let (positions, triangles) = build_mesh(centers, radii, spacing)?;
    let normals = vertex_normals(&positions, &triangles);
    Ok(Mesh {
        positions,
        normals,
        triangles,
    })
}

/// Area-weighted face normals averaged per vertex: cheap and smooth on fine meshes.
fn vertex_normals(positions: &[[f64; 3]], triangles: &[[u32; 3]]) -> Vec<[f64; 3]> {
    let mut normals = vec![[0.0; 3]; positions.len()];
    for t in triangles {
        let [a, b, c] = t.map(|i| positions[i as usize]);
        let n = cross(sub(b, a), sub(c, a));
        for &i in t {
            let acc = &mut normals[i as usize];
            (0..3).for_each(|k| acc[k] += n[k]);
        }
    }
    normals.into_iter().map(normalize_or_zero).collect()
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalize_or_zero(v: [f64; 3]) -> [f64; 3] {
    let norm = v.iter().map(|c| c * c).sum::<f64>().sqrt();
    if norm > 0.0 {
        v.map(|c| c / norm)
    } else {
        [0.0; 3]
    }
}

/// Raw mesh arrays: vertices and outward-wound triangles.
type MeshArrays = (Vec<[f64; 3]>, Vec<[u32; 3]>);

/// Mesher knobs, fixed at TMSmesh-style defaults.
struct GaussianMeshParams {
    /// Decay coefficient `d` in `exp(d · (1 − r²/R²))`. Range [0.1, 1.0]:
    /// smaller is smoother / more inflated. TMSmesh's default is 1.0.
    decay: f64,
    /// Isolevel ρ₀ defining the surface. With d = 1, ρ₀ = 1 puts the surface
    /// approximately on the van-der-Waals shell.
    isolevel: f64,
    /// Extra margin beyond `max(R) + cutoff_radius` on every axis, in Å. Stops
    /// the marching cubes step from emitting holes at the volume boundary.
    padding: f64,
    /// Atoms whose Gaussian contribution at a query point is below this are
    /// dropped from the density sum (via the KD-tree radius query).
    cutoff_eps: f64,
    /// Number of Taubin smoothing iterations to apply post-meshing.
    ///
    /// Default 3: marching cubes on a perfectly axis-aligned grid produces
    /// triangles with high-symmetry coplanarity (many panel pairs share
    /// edges that pass through other panel centroids), which triggers a
    /// log-divergence in BEMtzmann's WRG quadrature for the BEM matrix.
    /// A few Taubin iterations move vertices off the lattice without
    /// changing topology, breaking those alignments.
    ///
    /// Mirrors the design choice in TABI-PB
    /// (<https://github.com/Treecodes/TABI-PB>), which sets
    /// `Smooth_Mesh = true` in its NanoShaper configuration
    /// (`src/elements.cpp`) for the same BEM-conditioning reason.
    smoothing_iters: u32,
}

impl Default for GaussianMeshParams {
    fn default() -> Self {
        Self {
            decay: 1.0,
            isolevel: 1.0,
            padding: 2.0,
            cutoff_eps: 1e-6,
            smoothing_iters: 3,
        }
    }
}

/// Raw `(vertices, faces)`; manifold and outward-wound by construction.
fn build_mesh(positions: &[[f64; 3]], radii: &[f64], grid_spacing: f64) -> Result<MeshArrays> {
    let params = GaussianMeshParams::default();
    validate_inputs(positions, radii, grid_spacing)?;

    // why: precompute squared radii because `density_at` runs on every grid
    // corner (millions for a moderate protein) and looks up `r * r` per atom
    // hit on each call. Squaring once up front amortises that.
    let radii_sq: Vec<f64> = radii.iter().map(|r| r * r).collect();
    let max_radius = radii.iter().copied().fold(0.0_f64, f64::max);
    let max_cutoff_sq = density::cutoff_radius_sq(max_radius, params.decay, params.cutoff_eps);
    let max_cutoff = max_cutoff_sq.sqrt();

    let grid = Grid::fit(positions, radii, grid_spacing, max_cutoff, params.padding);
    let mut densities = sample_density_grid(positions, &radii_sq, &params, &grid, max_cutoff_sq)?;

    // why: corners landing exactly on the isolevel cause edge classification
    // to be ambiguous and let two distinct edges produce vertices at the same
    // point — degenerate triangles. Nudge them strictly above so every grid
    // corner is unambiguously inside (>= iso) or outside (< iso). Done once
    // here so every cube that shares the corner reads the same value.
    perturb_isolevel_ties(&mut densities, params.isolevel);

    let (mut vertices, faces) = marching_cubes::polygonize(&grid, &densities, params.isolevel);

    if vertices.is_empty() || faces.is_empty() {
        return Err(
            "marching cubes produced no triangles — check grid_spacing/decay vs atom radii"
                .to_string(),
        );
    }

    if params.smoothing_iters > 0 {
        smoothing::taubin_smooth(&mut vertices, &faces, params.smoothing_iters);
    }

    Ok((vertices, faces))
}

fn perturb_isolevel_ties(densities: &mut [f64], isolevel: f64) {
    const EPS: f64 = 1e-9;
    for d in densities {
        if (*d - isolevel).abs() < EPS {
            *d = isolevel + EPS;
        }
    }
}

fn validate_inputs(positions: &[[f64; 3]], radii: &[f64], grid_spacing: f64) -> Result<()> {
    if positions.len() != radii.len() {
        return Err(format!(
            "{} positions but {} radii",
            positions.len(),
            radii.len()
        ));
    }
    if positions.is_empty() {
        return Err("atoms list is empty".to_string());
    }
    if !grid_spacing.is_finite() || grid_spacing <= 0.0 {
        return Err(format!(
            "grid_spacing must be positive and finite, got {grid_spacing}"
        ));
    }
    for (i, &r) in radii.iter().enumerate() {
        if !r.is_finite() || r <= 0.0 {
            return Err(format!("atom {i} has non-positive radius {r}"));
        }
    }
    // why: the kd-tree's `add`/`within` return `NonFiniteCoordinate` on
    // NaN/Inf. Catch it here with a structured error rather than letting
    // the downstream call propagate a generic kd-tree failure.
    for (i, p) in positions.iter().enumerate() {
        if !p.iter().all(|c| c.is_finite()) {
            return Err(format!("atom {i} has non-finite position {p:?}"));
        }
    }
    Ok(())
}

fn sample_density_grid(
    positions: &[[f64; 3]],
    radii_sq: &[f64],
    params: &GaussianMeshParams,
    grid: &Grid,
    query_cutoff_sq: f64,
) -> Result<Vec<f64>> {
    let tree = density::build_atom_tree(positions)?;
    let plane = grid.nx * grid.ny;

    let mut densities = vec![0.0_f64; grid.n_samples()];
    for (k, slab) in densities.chunks_mut(plane).enumerate() {
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                let p = grid.coord(i, j, k);
                slab[i + grid.nx * j] =
                    density::density_at(p, radii_sq, params.decay, query_cutoff_sq, &tree)?;
            }
        }
    }
    Ok(densities)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn sphere_is_closed_and_on_its_radius() {
        const R: f64 = 5.0;
        let center = [1.0, -2.0, 3.0];
        let mesh = gaussian_surface(&[center], &[R]).unwrap();
        let half_step = 0.5; // spacing is 1 Å below 5000 spheres
        for p in &mesh.positions {
            let r = p
                .iter()
                .zip(center)
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt();
            assert!((r - R).abs() < half_step, "vertex at r = {r}");
        }
        let edges: HashSet<(u32, u32)> = mesh
            .triangles
            .iter()
            .flat_map(|&[a, b, c]| [(a, b), (b, c), (c, a)])
            .map(|(a, b)| (a.min(b), a.max(b)))
            .collect();
        let (v, e, f) = (mesh.positions.len(), edges.len(), mesh.triangles.len());
        assert_eq!(v as i64 - e as i64 + f as i64, 2, "Euler characteristic");
        assert_eq!(mesh.normals.len(), v);
    }
}
