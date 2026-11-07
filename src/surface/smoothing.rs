// Vendored from bemtzmann (https://github.com/mlund/bemtzmann, MIT OR Apache-2.0),
// src/geometry/molecular_surface/smoothing.rs at e3e3e0fa61cee83c67809133944ae97f05b6ca92.
// Changes across the vendored files: errors as String, serial density sampling (no
// rayon on wasm), unit tests and the test-only gradient dropped (covered upstream).

//! Taubin λ/μ mesh smoothing — internal helper, not part of the public API.
//!
//! Two-pass per iteration: shrink (λ > 0) then unshrink (μ < 0, |μ| > λ).
//! Defaults `λ = 0.33`, `μ = -0.34`. Faces (connectivity) unchanged; only
//! vertex positions move along the umbrella Laplacian.

// why: Taubin (1995) showed that a single Laplacian relaxation with positive
// weight λ shrinks the mesh; alternating with a slightly larger negative
// weight μ (so |μ| > λ) gives a low-pass filter that smooths without net
// shrinkage. 0.33 / -0.34 is the canonical pair from the original paper.
const LAMBDA: f64 = 0.33;
const MU: f64 = -0.34;

/// Apply `iters` Taubin iterations to `vertices` in place. No-op when
/// `iters == 0`. The neighbour set per vertex is built once from `faces`.
pub(super) fn taubin_smooth(vertices: &mut [[f64; 3]], faces: &[[u32; 3]], iters: u32) {
    if iters == 0 || vertices.len() < 3 {
        return;
    }
    let neighbours = vertex_neighbours(vertices.len(), faces);
    // why: `relax` reads each vertex's neighbours and writes a new position
    // for that vertex. Updating in place would let later vertices see
    // already-updated neighbours, biasing the result. Scratch holds the new
    // positions until the pass is complete, then a single `copy_from_slice`
    // commits them.
    let mut scratch = vec![[0.0_f64; 3]; vertices.len()];
    for _ in 0..iters {
        relax(vertices, &neighbours, &mut scratch, LAMBDA);
        relax(vertices, &neighbours, &mut scratch, MU);
    }
}

fn relax(
    vertices: &mut [[f64; 3]],
    neighbours: &[Vec<u32>],
    scratch: &mut [[f64; 3]],
    weight: f64,
) {
    for (i, nbrs) in neighbours.iter().enumerate() {
        if nbrs.is_empty() {
            scratch[i] = vertices[i];
            continue;
        }
        let inv_n = 1.0 / nbrs.len() as f64;
        let mut avg = [0.0_f64; 3];
        for &j in nbrs {
            let q = vertices[j as usize];
            avg[0] += q[0];
            avg[1] += q[1];
            avg[2] += q[2];
        }
        for v in &mut avg {
            *v *= inv_n;
        }
        let p = vertices[i];
        scratch[i] = [
            p[0] + weight * (avg[0] - p[0]),
            p[1] + weight * (avg[1] - p[1]),
            p[2] + weight * (avg[2] - p[2]),
        ];
    }
    vertices.copy_from_slice(scratch);
}

fn vertex_neighbours(n_verts: usize, faces: &[[u32; 3]]) -> Vec<Vec<u32>> {
    // Build adjacency directly into Vec<Vec<u32>>; each vertex shows up in
    // 6 face/edge slots on average, so duplicates are bounded. Sort+dedup
    // per vertex avoids the temporary hash-table memory of a HashSet build.
    let mut adjacency: Vec<Vec<u32>> = (0..n_verts).map(|_| Vec::new()).collect();
    for &[a, b, c] in faces {
        adjacency[a as usize].extend([b, c]);
        adjacency[b as usize].extend([a, c]);
        adjacency[c as usize].extend([a, b]);
    }
    for nbrs in &mut adjacency {
        nbrs.sort_unstable();
        nbrs.dedup();
    }
    adjacency
}
