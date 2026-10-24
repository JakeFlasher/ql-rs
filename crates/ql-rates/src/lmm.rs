use ql_ad::Real;
use ql_mc::PhiloxRng;

#[derive(Clone, Debug)]
pub struct ForwardMarketModel<T: Real> {
    pub tenor_times: Vec<f64>,
    pub initial_forwards: Vec<T>,
    pub initial_df_t0: T,
    pub displacements: Vec<T>,
    pub volatilities: Vec<T>,
    pub corr_decay: f64,
    pub corr_floor: f64,
    pub backward_looking_rfr: bool,
}
#[derive(Clone, Debug)]
pub struct TarnSpec<T: Real> {
    pub notional: T,
    pub target_coupon: T,
    pub fixed_rate: T,
    pub gearing: T,
    pub coupon_floor: T,
    pub coupon_cap: T,
    pub smoothing_eps: f64,
}
#[derive(Clone, Debug)]
pub struct CallableRangeAccrualSpec<T: Real> {
    pub notional: T,
    pub coupon_rate: T,
    pub lower_bound: T,
    pub upper_bound: T,
    pub call_period_indices: Vec<usize>,
    pub smoothing_eps: f64,
}
impl<T: Real> ForwardMarketModel<T> {
    pub fn discount_factor_0(&self, period_idx: usize) -> T {
        let i = period_idx.min(self.initial_forwards.len());
        let mut df = self.initial_df_t0;
        for j in 0..i {
            let tau = self.tenor_times[j + 1] - self.tenor_times[j];
            df /= T::one() + self.initial_forwards[j] * tau;
        }
        df
    }
    pub fn price_caplet_mc(&self, period_idx: usize, strike: T, n_paths: usize, seed: u64) -> T {
        if period_idx >= self.initial_forwards.len() || self.tenor_times.len() < period_idx + 2 {
            return T::zero();
        }
        let tau = self.tenor_times[period_idx + 1] - self.tenor_times[period_idx];
        let t = self.tenor_times[period_idx].max(0.0);
        let vol = self
            .volatilities
            .get(period_idx)
            .copied()
            .unwrap_or(T::zero());
        let shift = self
            .displacements
            .get(period_idx)
            .copied()
            .unwrap_or(T::zero());
        let f0 = self.initial_forwards[period_idx] + shift;
        let df = self.discount_factor_0(period_idx + 1);
        let rng = PhiloxRng::new(seed);
        let mut sum = T::zero();
        for p in 0..n_paths.max(1) {
            let z = T::from_f64(rng.normal4(p as u64, period_idx as u32)[0]);
            let vol_t = vol * T::from_f64(t.sqrt());
            let f = f0 * (vol_t * z - vol_t * vol_t * 0.5).exp() - shift;
            sum += (f - strike).max_val(T::zero()) * tau;
        }
        df * self.initial_df_t0 * sum / n_paths.max(1) as f64
    }
    pub fn price_tarn_mc(&self, spec: &TarnSpec<T>, n_paths: usize, seed: u64) -> T {
        let rng = PhiloxRng::new(seed);
        let mut total = T::zero();
        let m = self.initial_forwards.len();
        for p in 0..n_paths.max(1) {
            let mut accrued = T::zero();
            let mut pv = T::zero();
            for i in 0..m {
                let tau = self.tenor_times[i + 1] - self.tenor_times[i];
                let vol = self.volatilities.get(i).copied().unwrap_or(T::zero());
                let f = self.initial_forwards[i]
                    * (vol
                        * T::from_f64(self.tenor_times[i].sqrt())
                        * rng.normal4(p as u64, i as u32)[0]
                        - vol * vol * self.tenor_times[i] * 0.5)
                        .exp();
                let c = (spec.fixed_rate - spec.gearing * f)
                    .max_val(spec.coupon_floor)
                    .min_val(spec.coupon_cap);
                let remaining = (spec.target_coupon - accrued).max_val(T::zero());
                let paid = c.min_val(remaining);
                pv += paid * tau * (-spec.fixed_rate * self.tenor_times[i]).exp();
                accrued += paid;
                if accrued.val() >= spec.target_coupon.val() {
                    break;
                }
            }
            total += pv * spec.notional;
        }
        total / n_paths.max(1) as f64
    }
    pub fn price_callable_range_accrual_mc(
        &self,
        spec: &CallableRangeAccrualSpec<T>,
        n_paths: usize,
        seed: u64,
    ) -> T {
        let rng = PhiloxRng::new(seed);
        let mut total = T::zero();
        for p in 0..n_paths.max(1) {
            let mut pv = T::zero();
            for i in 0..self.initial_forwards.len() {
                let vol = self.volatilities.get(i).copied().unwrap_or(T::zero());
                let f = self.initial_forwards[i]
                    * (vol
                        * T::from_f64(self.tenor_times[i].sqrt())
                        * rng.normal4(p as u64, i as u32)[0]
                        - vol * vol * self.tenor_times[i] * 0.5)
                        .exp();
                let indicator = f.smooth_heaviside(spec.lower_bound, spec.smoothing_eps)
                    * (T::one() - f.smooth_heaviside(spec.upper_bound, spec.smoothing_eps));
                let tau = self.tenor_times[i + 1] - self.tenor_times[i];
                pv += spec.coupon_rate
                    * tau
                    * indicator
                    * (-spec.coupon_rate * self.tenor_times[i + 1]).exp();
                if spec.call_period_indices.contains(&(i + 1)) && pv.val() > 0.0 {
                    break;
                }
            }
            total += pv * spec.notional;
        }
        total / n_paths.max(1) as f64
    }
}
