use ql_ad::Real;
use ql_curves::YieldTermStructure;
use ql_models::HullWhite1F;
#[derive(Clone, Debug)]
pub struct NettingSetSwap {
    pub notional: f64,
    pub fixed_rate: f64,
    pub payment_times: Vec<f64>,
    pub accruals: Vec<f64>,
    pub is_payer: bool,
}
impl NettingSetSwap {
    pub fn npv_at<T: Real>(
        &self,
        model: &HullWhite1F<T>,
        curve: &impl YieldTermStructure<T>,
        t: f64,
        r_t: T,
    ) -> T {
        let mut fixed = T::zero();
        let mut last = None;
        for (i, &tm) in self.payment_times.iter().enumerate() {
            if tm > t {
                fixed += T::from_f64(self.fixed_rate * self.accruals.get(i).copied().unwrap_or(1.))
                    * model.zcb_price(curve, t, tm, r_t);
                last = Some(tm);
            }
        }
        let Some(last) = last else { return T::zero() };
        let floating = T::one() - model.zcb_price(curve, t, last, r_t);
        let v = if self.is_payer {
            floating - fixed
        } else {
            fixed - floating
        };
        T::from_f64(self.notional) * v
    }
}
#[derive(Clone, Debug)]
pub struct NpvCube<T: Real = f64> {
    pub n_netting_sets: usize,
    pub times: Vec<f64>,
    pub n_paths: usize,
    pub data: Vec<T>,
}
impl<T: Real> NpvCube<T> {
    pub fn new(ns: usize, times: Vec<f64>, paths: usize) -> Self {
        let nt = times.len();
        Self {
            n_netting_sets: ns,
            times,
            n_paths: paths,
            data: vec![T::zero(); ns * nt * paths],
        }
    }
    pub fn get(&self, ns: usize, time: usize, path: usize) -> T {
        assert!(ns < self.n_netting_sets, "netting-set index out of bounds");
        assert!(time < self.times.len(), "time index out of bounds");
        assert!(path < self.n_paths, "path index out of bounds");
        self.data[(ns * self.times.len() + time) * self.n_paths + path]
    }
    pub fn set(&mut self, ns: usize, time: usize, path: usize, v: T) {
        assert!(ns < self.n_netting_sets, "netting-set index out of bounds");
        assert!(time < self.times.len(), "time index out of bounds");
        assert!(path < self.n_paths, "path index out of bounds");
        let k = (ns * self.times.len() + time) * self.n_paths + path;
        self.data[k] = v
    }
}
pub fn simulate_hw1f_exposure_cube<T: Real>(
    model: &HullWhite1F<T>,
    curve: &impl YieldTermStructure<T>,
    swaps: &[NettingSetSwap],
    times: &[f64],
    n_paths: usize,
    seed: u64,
) -> NpvCube<T> {
    // The slice is one legal netting set.  Its trades are netted before the
    // exposure profile is written to the cube.
    let mut c = NpvCube::new(1, times.to_vec(), n_paths);
    let rng = ql_mc::PhiloxRng::new(seed);
    for pth in 0..n_paths {
        let mut x = T::zero();
        let mut prev = 0.;
        let sign = if pth % 2 == 0 { 1. } else { -1. };
        let pid = (pth / 2) as u64;
        for (k, &tm) in times.iter().enumerate() {
            let dt = tm - prev;
            let z = T::from_f64(rng.normal4(pid, k as u32)[0] * sign);
            let e = (-model.a * dt).exp();
            let var = if model.a.val().abs() < 1e-10 {
                T::from_f64(dt)
            } else {
                (T::one() - (-model.a * T::from_f64(2. * dt)).exp()) / (model.a * T::from_f64(2.))
            };
            x = e * x + model.sigma * var.sqrt() * z;
            let b = model.b_factor(0., tm);
            let alpha = curve.instantaneous_forward(tm)
                + T::from_f64(0.5) * model.sigma * model.sigma * b * b;
            let r = alpha + x;
            let net = swaps
                .iter()
                .map(|sw| sw.npv_at(model, curve, tm, r))
                .fold(T::zero(), |a, b| a + b);
            c.set(0, k, pth, net);
            prev = tm;
        }
    }
    c
}

/// Simulate several legally distinct netting sets.  Each inner vector is
/// netted independently and occupies one cube slice.
pub fn simulate_hw1f_multi_netting_set_cube<T: Real>(
    model: &HullWhite1F<T>,
    curve: &impl YieldTermStructure<T>,
    netting_sets: &[Vec<NettingSetSwap>],
    times: &[f64],
    n_paths: usize,
    seed: u64,
) -> NpvCube<T> {
    let mut c = NpvCube::new(netting_sets.len(), times.to_vec(), n_paths);
    let rng = ql_mc::PhiloxRng::new(seed);
    for pth in 0..n_paths {
        let mut x = T::zero();
        let mut prev = 0.0;
        let sign = if pth % 2 == 0 { 1.0 } else { -1.0 };
        let pid = (pth / 2) as u64;
        for (k, &tm) in times.iter().enumerate() {
            let dt = tm - prev;
            let z = T::from_f64(rng.normal4(pid, k as u32)[0] * sign);
            let e = (-model.a * dt).exp();
            let var = if model.a.val().abs() < 1e-10 {
                T::from_f64(dt)
            } else {
                (T::one() - (-model.a * (2.0 * dt)).exp()) / (model.a * 2.0)
            };
            x = e * x + model.sigma * var.sqrt() * z;
            let b = model.b_factor(0.0, tm);
            let alpha = curve.instantaneous_forward(tm) + model.sigma * model.sigma * b * b * 0.5;
            let r = alpha + x;
            for (ns, set) in netting_sets.iter().enumerate() {
                let net = set
                    .iter()
                    .map(|sw| sw.npv_at(model, curve, tm, r))
                    .fold(T::zero(), |a, b| a + b);
                c.set(ns, k, pth, net);
            }
            prev = tm;
        }
    }
    c
}
