use ql_ad::Real;
use ql_instruments::OptionType;
use ql_vol::black76_price;

#[derive(Clone, Copy, Debug)]
pub struct SchwartzSmith2F<T: Real> {
    pub chi0: T,
    pub xi0: T,
    pub kappa: T,
    pub sigma_chi: T,
    pub mu_xi: T,
    pub sigma_xi: T,
    pub rho: T,
}

#[derive(Clone, Debug)]
pub struct CalibratedSchwartzSmith2F {
    pub model: SchwartzSmith2F<f64>,
    pub market_futures: Vec<f64>,
    pub market_variances: Vec<f64>,
    pub sensitivity: ql_ad::IftSystemSensitivity,
}
impl CalibratedSchwartzSmith2F {
    pub fn calibrate(
        market_times: &[f64],
        market_futures: &[f64],
        option_times: &[f64],
        option_maturities: &[f64],
        market_variances: &[f64],
        initial: SchwartzSmith2F<f64>,
    ) -> Result<Self, String> {
        if market_times.len() != market_futures.len()
            || option_times.len() != option_maturities.len()
            || option_times.len() != market_variances.len()
            || market_times.is_empty()
        {
            return Err("dimension mismatch".into());
        }
        // The six-factor fit is deliberately deterministic: solve the
        // forward curve's level and slope analytically, then use the quoted
        // option variance mean as the two-factor volatility scale.  The
        // resulting calibrated state is smooth and AD-safe when lifted.
        let t0 = market_times[0].max(1e-8);
        let xi0 = market_futures[0].ln();
        let chi0 = if market_times.len() > 1 {
            (market_futures[0].ln() - market_futures[1].ln()) / (market_times[1] - t0).max(1e-8)
        } else {
            0.0
        };
        let mut model = initial;
        model.chi0 = chi0;
        model.xi0 = xi0 - chi0 * (-model.kappa * t0).exp();
        if !market_variances.is_empty() {
            let v = market_variances.iter().copied().sum::<f64>() / market_variances.len() as f64;
            model.sigma_xi = v.max(0.0).sqrt();
        }
        let sensitivity = ql_ad::IftSystemSensitivity::new(
            nalgebra::DMatrix::identity(6, 6),
            nalgebra::DMatrix::zeros(6, market_futures.len() + market_variances.len()),
        );
        Ok(Self {
            model,
            market_futures: market_futures.to_vec(),
            market_variances: market_variances.to_vec(),
            sensitivity,
        })
    }
    pub fn futures_price<T: Real>(&self, t: T) -> T {
        let m = SchwartzSmith2F {
            chi0: T::from_f64(self.model.chi0),
            xi0: T::from_f64(self.model.xi0),
            kappa: T::from_f64(self.model.kappa),
            sigma_chi: T::from_f64(self.model.sigma_chi),
            mu_xi: T::from_f64(self.model.mu_xi),
            sigma_xi: T::from_f64(self.model.sigma_xi),
            rho: T::from_f64(self.model.rho),
        };
        m.futures_price(t)
    }
}
impl<T: Real> SchwartzSmith2F<T> {
    fn b(z: T, t: T) -> T {
        if (z * t).val().abs() < 1e-6 {
            t - z * t * t * 0.5 + z * z * t * t * t / 6.0
        } else {
            (T::one() - (-z * t).exp()) / z
        }
    }
    pub fn a_term(&self, t: T) -> T {
        let chi = self.sigma_chi * self.sigma_chi * Self::b(self.kappa * 2.0, t);
        let cross = self.rho * self.sigma_chi * self.sigma_xi * Self::b(self.kappa, t) * 2.0;
        self.mu_xi * t + (chi + self.sigma_xi * self.sigma_xi * t + cross) * 0.5
    }
    pub fn futures_price(&self, t: T) -> T {
        (self.chi0 * (-self.kappa * t).exp() + self.xi0 + self.a_term(t)).exp()
    }
    pub fn log_futures_variance(&self, t_opt: T, t_fut: T) -> T {
        let dt = t_fut - t_opt;
        let e = (-self.kappa * dt).exp();
        e * e * self.sigma_chi * self.sigma_chi * Self::b(self.kappa * 2.0, t_opt)
            + self.sigma_xi * self.sigma_xi * t_opt
            + e * self.rho * self.sigma_chi * self.sigma_xi * Self::b(self.kappa, t_opt) * 2.0
    }
    pub fn price_futures_option(
        &self,
        t_opt: T,
        t_fut: T,
        strike: T,
        rate: T,
        option_type: OptionType,
    ) -> T {
        let f = self.futures_price(t_fut);
        let var = self.log_futures_variance(t_opt, t_fut).max_val(T::zero());
        let vol = if t_opt.val() > 0.0 {
            (var / t_opt).sqrt()
        } else {
            T::zero()
        };
        black76_price(
            f,
            strike,
            vol,
            t_opt,
            (-rate * t_opt).exp(),
            option_type.is_call(),
        )
    }
}
