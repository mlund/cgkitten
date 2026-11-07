// Vendored from bemtzmann (https://github.com/mlund/bemtzmann, MIT OR Apache-2.0),
// src/geometry/molecular_surface/grid.rs at e3e3e0fa61cee83c67809133944ae97f05b6ca92.
// Changes across the vendored files: errors as String, serial density sampling (no
// rayon on wasm), unit tests and the test-only gradient dropped (covered upstream).

//! Uniform 3D grid covering a padded bounding box around a set of weighted
//! spheres. Holds the geometry only — densities are stored separately in a
//! flat `Vec<f64>` indexed by `Grid::flat_idx`.

/// Axis-aligned uniform sampling grid.
pub(super) struct Grid {
    pub origin: [f64; 3],
    pub spacing: f64,
    /// Number of *samples* (corners) along each axis. The number of cubes
    /// (cells) is `nx-1`, `ny-1`, `nz-1`.
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
}

impl Grid {
    /// Build a grid that covers all atom spheres plus enough padding for the
    /// Gaussian density to fall below `cutoff_eps` at every boundary cell.
    ///
    /// `padding` is the per-axis margin to add beyond `max(R) + cutoff_radius`
    /// — a small extra buffer so the marching-cubes step never hits a cube
    /// with non-zero density at a boundary face (which would emit holes at the
    /// edge of the volume).
    pub(super) fn fit(
        positions: &[[f64; 3]],
        radii: &[f64],
        spacing: f64,
        max_cutoff: f64,
        padding: f64,
    ) -> Self {
        debug_assert!(spacing > 0.0 && padding >= 0.0 && max_cutoff > 0.0);
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for (p, &r) in positions.iter().zip(radii) {
            for axis in 0..3 {
                lo[axis] = lo[axis].min(p[axis] - r);
                hi[axis] = hi[axis].max(p[axis] + r);
            }
        }
        let total_pad = max_cutoff + padding;
        for (l, h) in lo.iter_mut().zip(hi.iter_mut()) {
            *l -= total_pad;
            *h += total_pad;
        }
        // why: snap origin down to a multiple of `spacing` so refinement
        // studies (halving the spacing) put the new grid corners on the same
        // analytical surface intersections as the coarser one — keeps
        // convergence-test residuals from being dominated by sub-cell phase.
        for l in &mut lo {
            *l = (*l / spacing).floor() * spacing;
        }
        // why: `+ 1` because we count *corners* (sample points), not *cells*.
        // Marching cubes consumes 8 corners per cell, so n_cells = n_corners − 1
        // along each axis; we need at least 2 corners per axis to form one cell.
        let nx = ((hi[0] - lo[0]) / spacing).ceil() as usize + 1;
        let ny = ((hi[1] - lo[1]) / spacing).ceil() as usize + 1;
        let nz = ((hi[2] - lo[2]) / spacing).ceil() as usize + 1;
        Self {
            origin: lo,
            spacing,
            nx,
            ny,
            nz,
        }
    }

    pub(super) const fn n_samples(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    pub(super) const fn flat_idx(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.nx * (j + self.ny * k)
    }

    pub(super) const fn coord(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        [
            self.origin[0] + (i as f64) * self.spacing,
            self.origin[1] + (j as f64) * self.spacing,
            self.origin[2] + (k as f64) * self.spacing,
        ]
    }
}
