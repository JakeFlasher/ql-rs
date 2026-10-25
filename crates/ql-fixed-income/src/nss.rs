use ql_ad::Real;
use ql_curves::YieldTermStructure;
use ql_time::{Date, DayCounter};
#[derive(Clone, Copy, Debug)]
pub struct NelsonSiegelSvensson<T: Real> {
    pub reference_date: Date,
    pub beta0: T,
    pub beta1: T,
    pub beta2: T,
    pub beta3: T,
    pub tau1: T,
    pub tau2: T,
}
impl<T: Real> NelsonSiegelSvensson<T> {
    fn h1(x: T) -> T {
        if x.val().abs() < 1e-6 {
            T::one() - x * 0.5 + x * x / 6. - x * x * x / 24.
        } else {
            (T::one() - (-x).exp()) / x
        }
    }
    pub fn zero_rate(&self, t: T) -> T {
        let tau1 = self.tau1.max_val(T::from_f64(1e-12));
        let tau2 = self.tau2.max_val(T::from_f64(1e-12));
        let x1 = t / tau1;
        let x2 = t / tau2;
        let h_1 = Self::h1(x1);
        let h_2 = h_1 - (-x1).exp();
        self.beta0 + self.beta1 * h_1 + self.beta2 * h_2 + self.beta3 * (Self::h1(x2) - (-x2).exp())
    }
    pub fn forward_rate(&self, t: T) -> T {
        let tau1 = self.tau1.max_val(T::from_f64(1e-12));
        let tau2 = self.tau2.max_val(T::from_f64(1e-12));
        let e1 = (-t / tau1).exp();
        let e2 = (-t / tau2).exp();
        self.beta0 + self.beta1 * e1 + self.beta2 * (t / tau1) * e1 + self.beta3 * (t / tau2) * e2
    }
    pub fn discount_factor(&self, t: T) -> T {
        (-self.zero_rate(t) * t).exp()
    }
    pub fn fit_ols_for_taus(times: &[f64], rates: &[f64], tau1: f64, tau2: f64) -> Self {
        use nalgebra::{DMatrix, DVector};
        assert_eq!(
            times.len(),
            rates.len(),
            "NSS times/rates dimension mismatch"
        );
        assert!(!times.is_empty(), "NSS fit requires observations");
        assert!(
            times.iter().chain(rates.iter()).all(|x| x.is_finite()),
            "NSS fit inputs must be finite"
        );
        let tau1 = tau1.max(1e-12);
        let tau2 = tau2.max(1e-12);
        let mut a = DMatrix::zeros(times.len(), 4);
        for (i, &t) in times.iter().enumerate() {
            let x1 = t / tau1;
            let x2 = t / tau2;
            // Keep the OLS design matrix finite at a zero/very short
            // maturity, matching the AD-safe `h1` used by `zero_rate`.
            let h1 = if x1.abs() < 1e-6 {
                1. - x1 * 0.5 + x1 * x1 / 6. - x1 * x1 * x1 / 24.
            } else {
                (1. - (-x1).exp()) / x1
            };
            let h2 = h1 - (-x1).exp();
            let h1_2 = if x2.abs() < 1e-6 {
                1. - x2 * 0.5 + x2 * x2 / 6. - x2 * x2 * x2 / 24.
            } else {
                (1. - (-x2).exp()) / x2
            };
            let h3 = h1_2 - (-x2).exp();
            a[(i, 0)] = 1.;
            a[(i, 1)] = h1;
            a[(i, 2)] = h2;
            a[(i, 3)] = h3;
        }
        let y = DVector::from_column_slice(rates);
        let b = a
            .svd(true, true)
            .solve(&y, 1e-12)
            .expect("NSS least-squares solve");
        Self {
            reference_date: Date(0),
            beta0: T::from_f64(b[0]),
            beta1: T::from_f64(b[1]),
            beta2: T::from_f64(b[2]),
            beta3: T::from_f64(b[3]),
            tau1: T::from_f64(tau1),
            tau2: T::from_f64(tau2),
        }
    }

    pub fn with_reference_date(mut self, reference_date: Date) -> Self {
        self.reference_date = reference_date;
        self
    }
}

impl<T: Real> YieldTermStructure<T> for NelsonSiegelSvensson<T> {
    fn reference_date(&self) -> Date {
        self.reference_date
    }
    fn discount_factor(&self, t: f64) -> T {
        NelsonSiegelSvensson::discount_factor(self, T::from_f64(t))
    }
    fn zero_rate(&self, t: f64) -> T {
        NelsonSiegelSvensson::zero_rate(self, T::from_f64(t))
    }
    fn instantaneous_forward(&self, t: f64) -> T {
        NelsonSiegelSvensson::forward_rate(self, T::from_f64(t))
    }
    fn forward_rate(&self, d1: Date, d2: Date, dc: DayCounter) -> T {
        let t1 = dc.year_fraction(self.reference_date, d1, None, None);
        let t2 = dc.year_fraction(self.reference_date, d2, None, None);
        let tau = dc.year_fraction(d1, d2, None, None);
        (<Self as YieldTermStructure<T>>::discount_factor(self, t1)
            / <Self as YieldTermStructure<T>>::discount_factor(self, t2)
            - T::one())
            / tau
    }
}

#[derive(Clone, Debug)]
pub struct TuckmanKeyRateDurationSpec<T: Real = f64> {
    pub nss_curve: NelsonSiegelSvensson<T>,
    pub cashflow_times: Vec<T>,
    pub cashflow_amounts: Vec<T>,
    pub key_rate_tenors: Vec<T>,
    pub pca_decay_lambda: T,
    pub benchmark_hedge_tenors: [T; 3],
}
#[derive(Clone, Debug)]
pub struct TuckmanKeyRateDurationReport<T: Real = f64> {
    pub dirty_price: T,
    pub effective_duration: T,
    pub effective_convexity: T,
    pub key_rate_durations: Vec<T>,
    pub key_rate_convexities: Vec<T>,
    pub key_rate_diagonal_convexities: Vec<T>,
    pub key_rate_dv01s: Vec<T>,
    pub krd_sum_partition_error: T,
    pub pca_level_exposure: T,
    pub pca_slope_exposure: T,
    pub pca_curvature_exposure: T,
    pub pca_butterfly_hedge_notionals: [T; 3],
    pub hedged_pca_residual_l2: T,
}
fn tuck_tent<T: Real>(t: T, k: usize, tenors: &[T]) -> T {
    let m = tenors.len();
    if m == 0 {
        return T::zero();
    }
    if k == 0 {
        if t.val() <= tenors[0].val() {
            T::one()
        } else if t.val() < tenors[1].val() {
            (tenors[1] - t) / (tenors[1] - tenors[0])
        } else {
            T::zero()
        }
    } else if k + 1 == m {
        if t.val() <= tenors[m - 2].val() {
            T::zero()
        } else if t.val() < tenors[m - 1].val() {
            (t - tenors[m - 2]) / (tenors[m - 1] - tenors[m - 2])
        } else {
            T::one()
        }
    } else if t.val() <= tenors[k - 1].val() || t.val() >= tenors[k + 1].val() {
        T::zero()
    } else if t.val() <= tenors[k].val() {
        (t - tenors[k - 1]) / (tenors[k] - tenors[k - 1])
    } else {
        (tenors[k + 1] - t) / (tenors[k + 1] - tenors[k])
    }
}
fn solve3<T: Real>(mut a: [[T; 3]; 3], mut b: [T; 3]) -> [T; 3] {
    let safe = |x: T| {
        if x.abs().val() < 1e-14 {
            if x.val() < 0.0 {
                T::from_f64(-1e-14)
            } else {
                T::from_f64(1e-14)
            }
        } else {
            x
        }
    };
    for i in 0..3 {
        let mut p = i;
        for r in i + 1..3 {
            if a[r][i].abs().val() > a[p][i].abs().val() {
                p = r;
            }
        }
        if p != i {
            a.swap(i, p);
            b.swap(i, p);
        }
        let d = safe(a[i][i]);
        for r in i + 1..3 {
            let f = a[r][i] / d;
            for c in i..3 {
                a[r][c] -= f * a[i][c];
            }
            b[r] -= f * b[i];
        }
    }
    let mut x = [T::zero(); 3];
    for i in (0..3).rev() {
        let mut s = b[i];
        for j in i + 1..3 {
            s -= a[i][j] * x[j];
        }
        x[i] = s / safe(a[i][i]);
    }
    x
}
pub fn analyze_tuckman_key_rate_durations_and_pca_hedge<T: Real>(
    s: &TuckmanKeyRateDurationSpec<T>,
) -> TuckmanKeyRateDurationReport<T> {
    let n = s.cashflow_times.len().min(s.cashflow_amounts.len());
    let m = s.key_rate_tenors.len();
    assert!(
        m >= 2
            && s.key_rate_tenors
                .windows(2)
                .all(|w| w[1].val() > w[0].val())
    );
    let mut pv = T::zero();
    let mut d = T::zero();
    let mut c = T::zero();
    let mut kr = vec![T::zero(); m];
    let mut kc = kr.clone();
    let mut kd = kr.clone();
    for i in 0..n {
        let t = s.cashflow_times[i];
        let p = s.cashflow_amounts[i] * (-s.nss_curve.zero_rate(t) * t).exp();
        pv += p;
        d += t * p;
        c += t * t * p;
        for k in 0..m {
            let w = tuck_tent(t, k, &s.key_rate_tenors);
            kr[k] += t * w * p;
            kc[k] += t * t * w * p;
            kd[k] += t * t * w * w * p;
        }
    }
    let pv_safe = pv.max_val(T::from_f64(1e-14));
    let effd = d / pv_safe;
    let effc = c / pv_safe;
    let dv: Vec<T> = kr.iter().map(|x| *x * T::from_f64(1e-4)).collect();
    let mut beta = [T::zero(); 3];
    let mut h = [[T::zero(); 3]; 3];
    for k in 0..m {
        let tau = s.key_rate_tenors[k];
        let l = [
            T::one(),
            (-s.pca_decay_lambda * tau).exp() - 0.5,
            tau * (-s.pca_decay_lambda * tau).exp() - 0.25,
        ];
        for f in 0..3 {
            beta[f] += l[f] * dv[k];
            for j in 0..3 {
                let th = s.benchmark_hedge_tenors[j];
                let ph = (-s.nss_curve.zero_rate(th) * th).exp();
                h[f][j] +=
                    l[f] * (T::from_f64(1e-4) * ph * th * tuck_tent(th, k, &s.key_rate_tenors));
            }
        }
    }
    // Normalize columns before elimination: short-tenor DV01 columns can be
    // several orders of magnitude smaller than long-tenor columns.
    let mut hs = h;
    let mut scales = [T::one(); 3];
    for j in 0..3 {
        scales[j] = (0..3)
            .map(|f| h[f][j].abs())
            .fold(T::zero(), |a, b| a.max_val(b))
            .max_val(T::from_f64(1e-14));
        for row in &mut hs {
            row[j] /= scales[j];
        }
    }
    let ys = solve3(hs, [-beta[0], -beta[1], -beta[2]]);
    let hn = [ys[0] / scales[0], ys[1] / scales[1], ys[2] / scales[2]];
    let mut res = T::zero();
    for f in 0..3 {
        let mut x = beta[f];
        for j in 0..3 {
            x += h[f][j] * hn[j];
        }
        res += x * x;
    }
    TuckmanKeyRateDurationReport {
        dirty_price: pv,
        effective_duration: effd,
        effective_convexity: effc,
        key_rate_durations: kr.iter().map(|x| *x / pv_safe).collect(),
        key_rate_convexities: kc.iter().map(|x| *x / pv_safe).collect(),
        key_rate_diagonal_convexities: kd.iter().map(|x| *x / pv_safe).collect(),
        key_rate_dv01s: dv,
        krd_sum_partition_error: (kr.iter().copied().fold(T::zero(), |a, b| a + b) / pv_safe
            - effd)
            .abs(),
        pca_level_exposure: beta[0],
        pca_slope_exposure: beta[1],
        pca_curvature_exposure: beta[2],
        pca_butterfly_hedge_notionals: hn,
        hedged_pca_residual_l2: res.sqrt(),
    }
}

#[cfg(test)]
mod phase22_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn tuckman_dual_fd_parity_and_hedge_invariants() {
        let make = |b0| TuckmanKeyRateDurationSpec {
            nss_curve: NelsonSiegelSvensson {
                reference_date: Date::from_ymd(2024, 1, 1).unwrap(),
                beta0: b0,
                beta1: -0.01,
                beta2: 0.02,
                beta3: -0.005,
                tau1: 1.5,
                tau2: 4.0,
            },
            cashflow_times: vec![1.0, 2.0, 3.0, 5.0],
            cashflow_amounts: vec![2.0, 2.0, 2.0, 102.0],
            key_rate_tenors: vec![1.0, 2.0, 5.0, 10.0],
            pca_decay_lambda: 0.25,
            benchmark_hedge_tenors: [2.0, 5.0, 10.0],
        };
        let r = analyze_tuckman_key_rate_durations_and_pca_hedge(&make(0.03));
        assert!(r.krd_sum_partition_error.val() < 1e-12);
        assert!(r.hedged_pca_residual_l2.val() < 1e-10);
        let h = 1e-5;
        let fd = (analyze_tuckman_key_rate_durations_and_pca_hedge(&make(0.03 + h)).dirty_price
            - analyze_tuckman_key_rate_durations_and_pca_hedge(&make(0.03 - h)).dirty_price)
            / (2.0 * h);
        let ad = analyze_tuckman_key_rate_durations_and_pca_hedge(&TuckmanKeyRateDurationSpec {
            nss_curve: NelsonSiegelSvensson {
                reference_date: Date::from_ymd(2024, 1, 1).unwrap(),
                beta0: Dual::<1>::variable(0.03, 0),
                beta1: Dual::<1>::from_value(-0.01),
                beta2: Dual::<1>::from_value(0.02),
                beta3: Dual::<1>::from_value(-0.005),
                tau1: Dual::<1>::from_value(1.5),
                tau2: Dual::<1>::from_value(4.0),
            },
            cashflow_times: vec![
                Dual::<1>::from_value(1.0),
                Dual::<1>::from_value(2.0),
                Dual::<1>::from_value(3.0),
                Dual::<1>::from_value(5.0),
            ],
            cashflow_amounts: vec![
                Dual::<1>::from_value(2.0),
                Dual::<1>::from_value(2.0),
                Dual::<1>::from_value(2.0),
                Dual::<1>::from_value(102.0),
            ],
            key_rate_tenors: vec![
                Dual::<1>::from_value(1.0),
                Dual::<1>::from_value(2.0),
                Dual::<1>::from_value(5.0),
                Dual::<1>::from_value(10.0),
            ],
            pca_decay_lambda: Dual::<1>::from_value(0.25),
            benchmark_hedge_tenors: [
                Dual::<1>::from_value(2.0),
                Dual::<1>::from_value(5.0),
                Dual::<1>::from_value(10.0),
            ],
        });
        assert!((ad.dirty_price.d[0] - fd).abs() < 1e-5);
    }
}
