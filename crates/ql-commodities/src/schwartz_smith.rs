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
    pub market_times: Vec<f64>,
    pub market_futures: Vec<f64>,
    pub option_times: Vec<f64>,
    pub option_maturities: Vec<f64>,
    pub market_variances: Vec<f64>,
    pub sensitivity: ql_ad::IftSystemSensitivity,
}
impl CalibratedSchwartzSmith2F {
    fn unpack(p: &[f64; 6], mu_xi: f64) -> SchwartzSmith2F<f64> {
        SchwartzSmith2F {
            chi0: p[0],
            xi0: p[1],
            kappa: p[2].max(1e-4),
            sigma_chi: p[3].max(1e-4),
            mu_xi,
            sigma_xi: p[4].max(1e-4),
            rho: p[5].clamp(-0.995, 0.995),
        }
    }
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
        let m_fut = market_futures.len();
        let m_opt = market_variances.len();
        let m_tot = m_fut + m_opt;
        let mu_xi = initial.mu_xi;
        let mut p = [
            initial.chi0,
            initial.xi0,
            initial.kappa.max(0.05),
            initial.sigma_chi.max(0.02),
            initial.sigma_xi.max(0.02),
            initial.rho.clamp(-0.95, 0.95),
        ];
        let eval_res = |params: &[f64; 6]| -> Vec<f64> {
            let mdl = Self::unpack(params, mu_xi);
            let mut r = Vec::with_capacity(m_tot);
            for (k, &t) in market_times.iter().enumerate() {
                r.push(mdl.futures_price(t) - market_futures[k]);
            }
            for k in 0..m_opt {
                r.push(
                    mdl.log_futures_variance(option_times[k], option_maturities[k])
                        - market_variances[k],
                );
            }
            r
        };
        let eval_jac = |params: &[f64; 6]| -> nalgebra::DMatrix<f64> {
            type D = ql_ad::Dual<6>;
            let mdl = SchwartzSmith2F {
                chi0: D::variable(params[0], 0),
                xi0: D::variable(params[1], 1),
                kappa: D::variable(params[2].max(1e-4), 2),
                sigma_chi: D::variable(params[3].max(1e-4), 3),
                mu_xi: D::from_value(mu_xi),
                sigma_xi: D::variable(params[4].max(1e-4), 4),
                rho: D::variable(params[5].clamp(-0.995, 0.995), 5),
            };
            let mut jac = nalgebra::DMatrix::zeros(m_tot, 6);
            for (k, &t) in market_times.iter().enumerate() {
                let v = mdl.futures_price(D::from_f64(t));
                for j in 0..6 {
                    jac[(k, j)] = v.d[j];
                }
            }
            for k in 0..m_opt {
                let v = mdl.log_futures_variance(
                    D::from_f64(option_times[k]),
                    D::from_f64(option_maturities[k]),
                );
                for j in 0..6 {
                    jac[(m_fut + k, j)] = v.d[j];
                }
            }
            jac
        };
        let mut damp = 1e-3_f64;
        for _ in 0..40 {
            let r = eval_res(&p);
            let r_vec = nalgebra::DVector::from_vec(r);
            let loss = r_vec.dot(&r_vec);
            if loss < 1e-22 {
                break;
            }
            let jac = eval_jac(&p);
            let jtj = jac.transpose() * &jac + nalgebra::DMatrix::identity(6, 6) * damp;
            let jtr = jac.transpose() * &r_vec;
            if let Some(step) = jtj.lu().solve(&jtr) {
                let mut cand = p;
                for j in 0..6 {
                    cand[j] -= step[j];
                }
                cand[2] = cand[2].max(1e-4);
                cand[3] = cand[3].max(1e-4);
                cand[4] = cand[4].max(1e-4);
                cand[5] = cand[5].clamp(-0.995, 0.995);
                let cand_r = nalgebra::DVector::from_vec(eval_res(&cand));
                if cand_r.dot(&cand_r) < loss {
                    p = cand;
                    damp = (damp * 0.5).max(1e-8);
                } else {
                    damp *= 4.0;
                }
            } else {
                damp *= 10.0;
            }
        }
        let model = Self::unpack(&p, mu_xi);
        let jac = eval_jac(&p);
        let j_theta = jac.transpose() * &jac + nalgebra::DMatrix::identity(6, 6) * 1e-8;
        let j_q = -jac.transpose();
        let sensitivity = ql_ad::IftSystemSensitivity::new(j_theta, j_q);
        Ok(Self {
            model,
            market_times: market_times.to_vec(),
            market_futures: market_futures.to_vec(),
            option_times: option_times.to_vec(),
            option_maturities: option_maturities.to_vec(),
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
    pub fn futures_price_and_quote_delta(&self, t: f64) -> (f64, Vec<f64>) {
        type D = ql_ad::Dual<6>;
        let m = SchwartzSmith2F {
            chi0: D::variable(self.model.chi0, 0),
            xi0: D::variable(self.model.xi0, 1),
            kappa: D::variable(self.model.kappa, 2),
            sigma_chi: D::variable(self.model.sigma_chi, 3),
            mu_xi: D::from_value(self.model.mu_xi),
            sigma_xi: D::variable(self.model.sigma_xi, 4),
            rho: D::variable(self.model.rho, 5),
        };
        let v = m.futures_price(D::from_f64(t));
        (v.val, self.sensitivity.vjp(&v.d))
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

#[cfg(test)]
mod phase14_tests {
    use super::*;
    #[test]
    fn phase14_schwartz_smith_calibration_is_live() {
        let m = SchwartzSmith2F {
            chi0: 0.1,
            xi0: 4.5,
            kappa: 1.2,
            sigma_chi: 0.3,
            mu_xi: 0.0,
            sigma_xi: 0.2,
            rho: -0.3,
        };
        let ts = vec![0.25, 0.5, 1.0, 2.0];
        let fs = ts.iter().map(|&t| m.futures_price(t)).collect::<Vec<_>>();
        let ot = vec![0.25, 0.5, 1.0];
        let om = vec![0.5, 1.0, 2.0];
        let ov = ot
            .iter()
            .zip(&om)
            .map(|(&a, &b)| m.log_futures_variance(a, b))
            .collect::<Vec<_>>();
        let c = CalibratedSchwartzSmith2F::calibrate(&ts, &fs, &ot, &om, &ov, m).unwrap();
        let (f1, delta) = c.futures_price_and_quote_delta(1.5);
        assert!((f1 - m.futures_price(1.5)).abs() < 1e-6);
        let h = 1e-5;
        let mut fs_up = fs.clone();
        let mut fs_dn = fs.clone();
        fs_up[1] += h;
        fs_dn[1] -= h;
        let c_up = CalibratedSchwartzSmith2F::calibrate(&ts, &fs_up, &ot, &om, &ov, m).unwrap();
        let c_dn = CalibratedSchwartzSmith2F::calibrate(&ts, &fs_dn, &ot, &om, &ov, m).unwrap();
        let fd = (c_up.futures_price(1.5) - c_dn.futures_price(1.5)) / (2.0 * h);
        assert!((delta[1] - fd).abs() < 1e-4);
    }
}
