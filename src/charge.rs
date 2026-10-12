//! Solution conditions and Henderson-Hasselbalch site charges.

use crate::residue::{PkaOverrides, TitratableGroup};

/// Solution conditions for charge calculation.
#[derive(Clone, Debug)]
pub struct Conditions {
    pub ph: f64,
    /// Temperature in Kelvin.
    pub temperature: f64,
    /// Ionic strength in mol/L.
    pub ionic_strength: f64,
    /// pKa overrides; groups not listed use their defaults.
    pub pka_overrides: PkaOverrides,
}

impl Conditions {
    /// Henderson-Hasselbalch deprotonated fraction of `group`, honouring overrides.
    pub fn fraction_deprotonated(&self, group: &TitratableGroup) -> f64 {
        1.0 / (1.0 + 10_f64.powf(self.pka_overrides.effective(group) - self.ph))
    }

    /// Henderson-Hasselbalch mean charge of `group`.
    pub fn mean_charge(&self, group: &TitratableGroup) -> f64 {
        let f = self.fraction_deprotonated(group);
        f.mul_add(
            group.charge_deprotonated - group.charge_protonated,
            group.charge_protonated,
        )
    }

    /// Henderson-Hasselbalch charge variance of `group`, (Δq)² f(1−f).
    pub fn charge_variance(&self, group: &TitratableGroup) -> f64 {
        let f = self.fraction_deprotonated(group);
        (group.charge_deprotonated - group.charge_protonated).powi(2) * f * (1.0 - f)
    }
}

impl Default for Conditions {
    fn default() -> Self {
        Self {
            ph: 7.0,
            temperature: 298.15,
            ionic_strength: 0.1,
            pka_overrides: PkaOverrides::default(),
        }
    }
}

impl From<&crate::ChargeCalc> for Conditions {
    fn from(calc: &crate::ChargeCalc) -> Self {
        Self {
            ph: calc.ph,
            temperature: calc.temperature,
            ionic_strength: calc.ionic_strength,
            pka_overrides: calc.pka_overrides.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::residue::TITRATABLE_GROUPS;

    fn conditions_at_ph(ph: f64) -> Conditions {
        Conditions {
            ph,
            ..Default::default()
        }
    }

    #[test]
    fn at_pka_charge_is_halfway() {
        for group in TITRATABLE_GROUPS {
            let charge = conditions_at_ph(group.pka).mean_charge(group);
            let expected = (group.charge_protonated + group.charge_deprotonated) / 2.0;
            assert!(
                (charge - expected).abs() < 1e-10,
                "{}: charge at pKa = {charge}, expected {expected}",
                group.res_name,
            );
        }
    }

    #[test]
    fn pka_override_shifts_charge() {
        let his = crate::residue::find_titratable_group("HIS").unwrap();
        let mut conditions = conditions_at_ph(7.0);
        let default = conditions.mean_charge(his);
        conditions.pka_overrides.insert("HIS", 7.0).unwrap();
        assert!((conditions.mean_charge(his) - 0.5).abs() < 1e-10);
        assert!(default < 0.5);
        // Overrides for other groups leave HIS untouched
        conditions.pka_overrides = PkaOverrides::default();
        conditions.pka_overrides.insert("LYS", 5.0).unwrap();
        assert_eq!(conditions.mean_charge(his), default);
    }

    #[test]
    fn asp_charges_at_extremes() {
        let asp = TITRATABLE_GROUPS
            .iter()
            .find(|g| g.res_name == "ASP")
            .unwrap();

        let charge_low = conditions_at_ph(1.0).mean_charge(asp);
        assert!(charge_low.abs() < 0.01, "ASP at pH 1: {charge_low}");

        let charge_high = conditions_at_ph(10.0).mean_charge(asp);
        assert!(
            (charge_high + 1.0).abs() < 0.01,
            "ASP at pH 10: {charge_high}"
        );
    }

    #[test]
    fn lys_charges_at_extremes() {
        let lys = TITRATABLE_GROUPS
            .iter()
            .find(|g| g.res_name == "LYS")
            .unwrap();

        let charge = conditions_at_ph(7.0).mean_charge(lys);
        assert!((charge - 1.0).abs() < 0.01, "LYS at pH 7: {charge}");

        let charge_high = conditions_at_ph(14.0).mean_charge(lys);
        assert!(charge_high.abs() < 0.01, "LYS at pH 14: {charge_high}");
    }
}
