//! Running charge moments about a fixed reference point.
//!
//! Dipole and quadrupole of a charged molecule depend on the origin; the centre of
//! mass is used (as in faunus) so values do not depend on the PDB frame.

use crate::Bead;

/// Mass-weighted centre, falling back to the geometric centre for massless input.
pub(crate) fn center_of_mass(beads: &[Bead]) -> [f64; 3] {
    let total: f64 = beads.iter().map(|b| b.mass).sum();
    let massless = total <= 0.0;
    let norm = if massless {
        beads.len().max(1) as f64
    } else {
        total
    };
    let mut c = [0.0; 3];
    for b in beads {
        let wi = if massless { 1.0 } else { b.mass };
        c[0] += wi * b.x;
        c[1] += wi * b.y;
        c[2] += wi * b.z;
    }
    c.map(|v| v / norm)
}

/// Net charge, dipole vector and primitive second moment Σ q d dᵀ.
///
/// Additive in charge, so an accepted MC move updates it in O(1).
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Moments {
    pub z: f64,
    mu: [f64; 3],
    /// Upper triangle of Σ q d dᵀ: xx, xy, xz, yy, yz, zz.
    p: [f64; 6],
}

impl Moments {
    /// Moments of `charges` at `positions`, both relative to the reference.
    pub fn from_charges<'a>(items: impl IntoIterator<Item = (f64, &'a [f64; 3])>) -> Self {
        let mut m = Self::default();
        for (q, d) in items {
            m.add(q, d);
        }
        m
    }

    pub fn add(&mut self, q: f64, d: &[f64; 3]) {
        let [x, y, z] = *d;
        self.z += q;
        self.mu[0] += q * x;
        self.mu[1] += q * y;
        self.mu[2] += q * z;
        self.p[0] += q * x * x;
        self.p[1] += q * x * y;
        self.p[2] += q * x * z;
        self.p[3] += q * y * y;
        self.p[4] += q * y * z;
        self.p[5] += q * z * z;
    }

    pub fn dipole_sq(&self) -> f64 {
        self.mu.iter().map(|m| m * m).sum()
    }

    /// Squared Frobenius norm of the traceless (Buckingham) quadrupole
    /// Θ = ½(3P − tr P·I); traceless so isotropic charge shells give zero.
    pub fn quadrupole_sq(&self) -> f64 {
        let [xx, xy, xz, yy, yz, zz] = self.p;
        let tr = xx + yy + zz;
        let diag: f64 = [xx, yy, zz]
            .iter()
            .map(|&a| (0.5 * 3.0f64.mul_add(a, -tr)).powi(2))
            .sum();
        let off: f64 = [xy, xz, yz].iter().map(|&a| (1.5 * a).powi(2)).sum();
        2.0f64.mul_add(off, diag)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_quadrupole_matches_buckingham() {
        // +q at ±a on z, −2q at origin: Θ_zz = 2qa², Θ_xx = Θ_yy = −qa²
        let (q, a) = (1.0, 2.0);
        let m = Moments::from_charges([
            (q, &[0.0, 0.0, a]),
            (q, &[0.0, 0.0, -a]),
            (-2.0 * q, &[0.0, 0.0, 0.0]),
        ]);
        let expected = 6.0 * (q * a * a).powi(2);
        assert!((m.quadrupole_sq() - expected).abs() < 1e-12);
        assert!(m.dipole_sq().abs() < 1e-12);
        assert!(m.z.abs() < 1e-12);
    }

    #[test]
    fn isotropic_shell_has_no_quadrupole() {
        let m = Moments::from_charges(
            [
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
            ]
            .iter()
            .map(|d| (1.0, d)),
        );
        assert!(m.quadrupole_sq().abs() < 1e-12);
    }
}
