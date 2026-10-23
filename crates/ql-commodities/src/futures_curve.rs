use ql_ad::{Dual, Real};

#[derive(Clone, Copy, Debug)]
pub struct SeasonalParams<T: Real> {
    pub a1: T,
    pub b1: T,
    pub a2: T,
    pub b2: T,
}
impl<T: Real> SeasonalParams<T> {
    pub fn factor(&self, t: T) -> T {
        self.a1 * (T::from_f64(2.0 * std::f64::consts::PI) * t).cos()
            + self.b1 * (T::from_f64(2.0 * std::f64::consts::PI) * t).sin()
            + self.a2 * (T::from_f64(4.0 * std::f64::consts::PI) * t).cos()
            + self.b2 * (T::from_f64(4.0 * std::f64::consts::PI) * t).sin()
    }
}

#[derive(Clone, Debug)]
pub struct CommodityFuturesCurve<T: Real> {
    pub times: Vec<T>,
    pub log_base_forwards: Vec<T>,
    pub seasonal: SeasonalParams<T>,
}
impl<T: Real> CommodityFuturesCurve<T> {
    pub fn from_futures_pillars(
        times: Vec<T>,
        market_futures: &[T],
        seasonal: SeasonalParams<T>,
    ) -> Self {
        assert_eq!(times.len(), market_futures.len());
        assert!(!times.is_empty());
        let log_base_forwards = times
            .iter()
            .zip(market_futures)
            .map(|(t, f)| f.ln() - seasonal.factor(*t))
            .collect();
        Self {
            times,
            log_base_forwards,
            seasonal,
        }
    }
    pub fn futures_price(&self, t: T) -> T {
        let n = self.times.len();
        let i = if t.val() <= self.times[0].val() {
            0
        } else if t.val() >= self.times[n - 1].val() {
            n - 1
        } else {
            self.times
                .partition_point(|x| x.val() <= t.val())
                .saturating_sub(1)
        };
        let base = if i + 1 >= n {
            self.log_base_forwards[i]
        } else {
            let w = (t - self.times[i]) / (self.times[i + 1] - self.times[i]);
            self.log_base_forwards[i]
                + w * (self.log_base_forwards[i + 1] - self.log_base_forwards[i])
        };
        (base + self.seasonal.factor(t)).exp()
    }
    pub fn implied_convenience_yield(&self, spot: T, rate: T, storage_cost: T, t: T) -> T {
        rate + storage_cost - (self.futures_price(t) / spot).ln() / t
    }
}

#[derive(Clone, Debug)]
pub struct CalibratedCommodityCurve {
    pub times: Vec<f64>,
    pub market_futures: Vec<f64>,
    pub log_base_forwards: Vec<f64>,
    pub seasonal: SeasonalParams<f64>,
}
impl CalibratedCommodityCurve {
    pub fn new(times: Vec<f64>, market_futures: Vec<f64>, seasonal: SeasonalParams<f64>) -> Self {
        assert_eq!(times.len(), market_futures.len());
        let log_base_forwards = times
            .iter()
            .zip(&market_futures)
            .map(|(&t, &f)| f.ln() - seasonal.factor(t))
            .collect();
        Self {
            times,
            market_futures,
            log_base_forwards,
            seasonal,
        }
    }
    pub fn curve(&self) -> CommodityFuturesCurve<f64> {
        CommodityFuturesCurve::from_futures_pillars(
            self.times.clone(),
            &self.market_futures,
            self.seasonal,
        )
    }
    pub fn dp_dfutures<F>(&self, pricer: F) -> (f64, Vec<f64>)
    where
        F: Fn(&CommodityFuturesCurve<Dual<1>>) -> Dual<1>,
    {
        let mut out = Vec::with_capacity(self.market_futures.len());
        for j in 0..self.market_futures.len() {
            let times = self
                .times
                .iter()
                .map(|&x| Dual::from_value(x))
                .collect::<Vec<_>>();
            let theta = self
                .log_base_forwards
                .iter()
                .enumerate()
                .map(|(i, &x)| {
                    if i == j {
                        Dual::variable(x, 0)
                    } else {
                        Dual::from_value(x)
                    }
                })
                .collect::<Vec<_>>();
            let seasonal = SeasonalParams {
                a1: Dual::from_value(self.seasonal.a1),
                b1: Dual::from_value(self.seasonal.b1),
                a2: Dual::from_value(self.seasonal.a2),
                b2: Dual::from_value(self.seasonal.b2),
            };
            out.push(
                pricer(&CommodityFuturesCurve {
                    times,
                    log_base_forwards: theta,
                    seasonal,
                })
                .d[0]
                    / self.market_futures[j],
            );
        }
        let curve = CommodityFuturesCurve {
            times: self.times.iter().map(|&x| Dual::from_value(x)).collect(),
            log_base_forwards: self
                .log_base_forwards
                .iter()
                .map(|&x| Dual::from_value(x))
                .collect(),
            seasonal: SeasonalParams {
                a1: Dual::from_value(self.seasonal.a1),
                b1: Dual::from_value(self.seasonal.b1),
                a2: Dual::from_value(self.seasonal.a2),
                b2: Dual::from_value(self.seasonal.b2),
            },
        };
        (pricer(&curve).val, out)
    }
}
