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
impl<T: Real> SchwartzSmith2F<T> {
    fn b(z: T, t: T) -> T {
        if z.val().abs() < 1e-12 {
            t
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
