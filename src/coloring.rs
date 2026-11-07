//! Colour maps for charge and potential, shared by viewers so they read alike.

/// Display (sRGB) RGBA, 0–1, as egui and the viewers' render targets take it.
pub type Rgba = [f32; 4];

const NEUTRAL: [f32; 3] = [0.80, 0.80, 0.80];
const POSITIVE: [f32; 3] = [0.17, 0.32, 0.86];
const NEGATIVE: [f32; 3] = [0.84, 0.15, 0.15];

/// Red ← grey → blue for `t` in −1..1 (clamped). Grey rather than white so neutral
/// regions keep their shading.
pub fn diverging(t: f64) -> Rgba {
    let t = t.clamp(-1.0, 1.0) as f32;
    let end = if t >= 0.0 { POSITIVE } else { NEGATIVE };
    mix(NEUTRAL, end, t.abs())
}

/// Grey → gold → magenta for `t` in 0..1 (clamped); distinct from the red/blue map so
/// a fluctuation isn't read as a sign.
pub fn sequential(t: f64) -> Rgba {
    const STOPS: [[f32; 3]; 3] = [[0.80, 0.80, 0.80], [0.96, 0.76, 0.18], [0.78, 0.13, 0.62]];
    let t = 2.0 * t.clamp(0.0, 1.0) as f32;
    if t <= 1.0 {
        mix(STOPS[0], STOPS[1], t)
    } else {
        mix(STOPS[1], STOPS[2], t - 1.0)
    }
}

/// The `quantile` of `values`, at least `floor`: the value that saturates a colour
/// map so most of a surface uses its range, while a near-neutral one stays pale.
pub fn percentile_scale(values: impl IntoIterator<Item = f64>, quantile: f64, floor: f64) -> f64 {
    let mut v: Vec<f64> = values.into_iter().collect();
    v.sort_by(f64::total_cmp);
    let i = ((v.len() as f64 * quantile) as usize).min(v.len().saturating_sub(1));
    v.get(i).copied().unwrap_or(0.0).max(floor)
}

/// (1 − f)a + fb rather than a + f(b − a), so f = 1 gives `b` exactly.
fn mix(a: [f32; 3], b: [f32; 3], f: f32) -> Rgba {
    let c = |i: usize| (1.0 - f) * a[i] + f * b[i];
    [c(0), c(1), c(2), 1.0]
}

#[cfg(test)]
mod tests {
    use super::*;

    const fn opaque([r, g, b]: [f32; 3]) -> Rgba {
        [r, g, b, 1.0]
    }
    const GREY: Rgba = opaque(NEUTRAL);
    const BLUE: Rgba = opaque(POSITIVE);
    const RED: Rgba = opaque(NEGATIVE);

    #[test]
    fn diverging_is_grey_at_zero_and_clamps_at_the_ends() {
        assert_eq!(diverging(0.0), GREY);
        assert_eq!(diverging(1.0), BLUE);
        assert_eq!(diverging(-1.0), RED);
        assert_eq!(diverging(7.0), BLUE);
        assert_eq!(diverging(-7.0), RED);
    }
}
