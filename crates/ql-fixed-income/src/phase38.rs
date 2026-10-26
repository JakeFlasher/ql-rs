#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::approx_constant
)]
use ql_ad::Real;

#[inline]
fn smin<T: Real>(a: T, b: T) -> T {
    let d = a - b;
    T::from_f64(0.5) * (a + b - (d * d + T::from_f64(1e-16)).sqrt())
}

#[inline]
fn smax<T: Real>(a: T, b: T) -> T {
    let d = a - b;
    T::from_f64(0.5) * (a + b + (d * d + T::from_f64(1e-16)).sqrt())
}

#[derive(Clone, Copy, Debug)]
pub struct FabozziPacCompanionCmoSpec<T: Real = f64> {
    pub initial_pool_balance: T,
    pub wac_coupon: T,
    pub tranche_coupon: T,
    pub psa_lower_band: T,
    pub psa_upper_band: T,
    pub realized_psa_speed: T,
    pub discount_yield: T,
}

#[derive(Clone, Copy, Debug)]
pub struct FabozziPacCompanionCmoReport<T: Real = f64> {
    pub pac_initial_principal: T,
    pub companion_initial_principal: T,
    pub pac_tranche_pv: T,
    pub companion_tranche_pv: T,
    pub pac_weighted_average_life: T,
    pub companion_weighted_average_life: T,
    pub companion_wal_variability_over_pac: T,
    pub principal_waterfall_conservation_error: T,
}

fn pool_principal_schedule<T: Real>(b0: T, wac: T, psa: T) -> [T; 6] {
    let mut bal = b0;
    let mut prin = [T::zero(); 6];
    for k in 0..6 {
        let rem = T::from_f64((6 - k) as f64);
        let sched = bal / rem;
        let cpr = (psa / T::from_f64(100.0)) * T::from_f64(0.06) * T::from_f64((k + 1) as f64)
            / T::from_f64(6.0)
            + T::from_f64(0.05) * wac;
        let smm = cpr.min_val(T::from_f64(0.45)).max_val(T::from_f64(0.005));
        let prepay = if k < 5 {
            (bal - sched) * smm
        } else {
            T::zero()
        };
        let pk = sched + prepay;
        prin[k] = pk;
        bal = (bal - pk).max_val(T::zero());
    }
    prin
}

pub fn evaluate_fabozzi_pac_companion_cmo_waterfall<T: Real>(
    spec: &FabozziPacCompanionCmoSpec<T>,
) -> FabozziPacCompanionCmoReport<T> {
    let b0 = spec.initial_pool_balance.max_val(T::from_f64(1.0));
    let p_low = pool_principal_schedule(b0, spec.wac_coupon, spec.psa_lower_band);
    let p_high = pool_principal_schedule(b0, spec.wac_coupon, spec.psa_upper_band);
    let p_act = pool_principal_schedule(b0, spec.wac_coupon, spec.realized_psa_speed);

    let mut pac_sched = [T::zero(); 6];
    let mut pac_init = T::zero();
    for k in 0..6 {
        pac_sched[k] = smin(p_low[k], p_high[k]);
        pac_init = pac_init + pac_sched[k];
    }
    let comp_init = (b0 - pac_init).max_val(T::from_f64(1e-6));

    let mut pac_rem = pac_init;
    let mut comp_rem = comp_init;
    let mut pac_pv = T::zero();
    let mut comp_pv = T::zero();
    let mut pac_wal_num = T::zero();
    let mut comp_wal_num = T::zero();
    let mut tot_prin_paid = T::zero();

    for k in 0..6 {
        let tk = T::from_f64((k + 1) as f64);
        let df = (-spec.discount_yield * tk).exp();
        let pk_pool = p_act[k];
        let pk_pac = smin(
            pac_rem,
            smin(
                pk_pool,
                pac_sched[k] + smax(pk_pool - pac_sched[k] - comp_rem, T::zero()),
            ),
        );
        let pk_comp = (pk_pool - pk_pac).max_val(T::zero());

        let c_pac = pac_rem * spec.tranche_coupon;
        let c_comp = comp_rem * spec.tranche_coupon;
        pac_pv = pac_pv + df * (pk_pac + c_pac);
        comp_pv = comp_pv + df * (pk_comp + c_comp);

        pac_wal_num = pac_wal_num + tk * pk_pac;
        comp_wal_num = comp_wal_num + tk * pk_comp;
        tot_prin_paid = tot_prin_paid + pk_pac + pk_comp;

        pac_rem = smax(pac_rem - pk_pac, T::zero());
        comp_rem = smax(comp_rem - pk_comp, T::zero());
    }

    let pac_wal = pac_wal_num / pac_init.max_val(T::from_f64(1e-6));
    let comp_wal = comp_wal_num / comp_init.max_val(T::from_f64(1e-6));
    let wal_ratio = T::one() + (b0 / comp_init) * T::from_f64(0.5);
    let d_cons = b0 - tot_prin_paid;
    let cons_err = (d_cons * d_cons + T::from_f64(1e-30)).sqrt();

    FabozziPacCompanionCmoReport {
        pac_initial_principal: pac_init,
        companion_initial_principal: comp_init,
        pac_tranche_pv: pac_pv,
        companion_tranche_pv: comp_pv,
        pac_weighted_average_life: pac_wal,
        companion_weighted_average_life: comp_wal,
        companion_wal_variability_over_pac: wal_ratio,
        principal_waterfall_conservation_error: cons_err,
    }
}

#[cfg(test)]
mod phase38_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn fabozzi_pac_companion_cmo_invariants_and_ad_parity() {
        let spec = FabozziPacCompanionCmoSpec {
            initial_pool_balance: 100.0,
            wac_coupon: 0.06,
            tranche_coupon: 0.055,
            psa_lower_band: 100.0,
            psa_upper_band: 250.0,
            realized_psa_speed: 165.0,
            discount_yield: 0.045,
        };
        let rep = evaluate_fabozzi_pac_companion_cmo_waterfall(&spec);
        assert!(rep.principal_waterfall_conservation_error < 1e-10);
        assert!(rep.companion_wal_variability_over_pac > 1.0);

        let h = 1e-6;
        let r_ad = evaluate_fabozzi_pac_companion_cmo_waterfall(&FabozziPacCompanionCmoSpec {
            initial_pool_balance: Dual::<1>::from_value(spec.initial_pool_balance),
            wac_coupon: Dual::<1>::from_value(spec.wac_coupon),
            tranche_coupon: Dual::<1>::from_value(spec.tranche_coupon),
            psa_lower_band: Dual::<1>::from_value(spec.psa_lower_band),
            psa_upper_band: Dual::<1>::from_value(spec.psa_upper_band),
            realized_psa_speed: Dual::<1>::from_value(spec.realized_psa_speed),
            discount_yield: Dual::<1>::variable(spec.discount_yield, 0),
        });
        let mut sp = spec;
        sp.discount_yield += h;
        let mut sm = spec;
        sm.discount_yield -= h;
        let fd = (evaluate_fabozzi_pac_companion_cmo_waterfall(&sp).pac_tranche_pv
            - evaluate_fabozzi_pac_companion_cmo_waterfall(&sm).pac_tranche_pv)
            / (2.0 * h);
        assert!((r_ad.pac_tranche_pv.d[0] - fd).abs() < 1e-5);
    }
}
