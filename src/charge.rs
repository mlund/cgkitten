//! Charge calculation using the Henderson-Hasselbalch equation.
//!
//! Provides a [`ChargeCalculator`] trait for extensibility, with a default
//! implementation using Henderson-Hasselbalch.

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

/// Trait for computing the charge of a titratable group under given conditions.
pub trait ChargeCalculator {
    fn charge(&self, group: &TitratableGroup, conditions: &Conditions) -> f64;
}

/// Henderson-Hasselbalch charge calculator.
///
/// For an acid HA <-> H+ + A-:
///   fraction_deprotonated = 1 / (1 + 10^(pKa - pH))
///   charge = f_deprot * q_deprot + (1 - f_deprot) * q_prot
#[derive(Clone, Debug, Default)]
pub struct HendersonHasselbalch;

impl ChargeCalculator for HendersonHasselbalch {
    fn charge(&self, group: &TitratableGroup, conditions: &Conditions) -> f64 {
        let fraction_deprotonated = conditions.fraction_deprotonated(group);
        fraction_deprotonated * group.charge_deprotonated
            + (1.0 - fraction_deprotonated) * group.charge_protonated
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
        let hh = HendersonHasselbalch;
        for group in TITRATABLE_GROUPS {
            let charge = hh.charge(group, &conditions_at_ph(group.pka));
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
        let hh = HendersonHasselbalch;
        let his = crate::residue::find_titratable_group("HIS").unwrap();
        let mut conditions = conditions_at_ph(7.0);
        let default = hh.charge(his, &conditions);
        conditions.pka_overrides.insert("HIS", 7.0).unwrap();
        assert!((hh.charge(his, &conditions) - 0.5).abs() < 1e-10);
        assert!(default < 0.5);
        // Overrides for other groups leave HIS untouched
        conditions.pka_overrides = PkaOverrides::default();
        conditions.pka_overrides.insert("LYS", 5.0).unwrap();
        assert_eq!(hh.charge(his, &conditions), default);
    }

    #[test]
    fn asp_charges_at_extremes() {
        let hh = HendersonHasselbalch;
        let asp = TITRATABLE_GROUPS
            .iter()
            .find(|g| g.res_name == "ASP")
            .unwrap();

        let charge_low = hh.charge(asp, &conditions_at_ph(1.0));
        assert!(charge_low.abs() < 0.01, "ASP at pH 1: {charge_low}");

        let charge_high = hh.charge(asp, &conditions_at_ph(10.0));
        assert!(
            (charge_high + 1.0).abs() < 0.01,
            "ASP at pH 10: {charge_high}"
        );
    }

    #[test]
    fn lys_charges_at_extremes() {
        let hh = HendersonHasselbalch;
        let lys = TITRATABLE_GROUPS
            .iter()
            .find(|g| g.res_name == "LYS")
            .unwrap();

        let charge = hh.charge(lys, &conditions_at_ph(7.0));
        assert!((charge - 1.0).abs() < 0.01, "LYS at pH 7: {charge}");

        let charge_high = hh.charge(lys, &conditions_at_ph(14.0));
        assert!(charge_high.abs() < 0.01, "LYS at pH 14: {charge_high}");
    }
}
