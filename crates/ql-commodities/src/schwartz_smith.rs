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

#[derive(Copy, Clone, Debug)]
pub struct GibsonSchwartz2F<T: Real = f64> {
    pub spot: T,
    pub convenience_yield0: T,
    pub rate: T,
    pub kappa: T,
    pub alpha_rn: T,
    pub sigma_spot: T,
    pub sigma_cy: T,
    pub rho: T,
}
impl<T: Real> GibsonSchwartz2F<T> {
    fn b(&self, t: T) -> T {
        let at = self.kappa * t;
        if at.val().abs() < 1e-6 {
            t * (T::one() - at * 0.5 + at * at / 6.0)
        } else {
            (T::one() - (-at).exp()) / self.kappa
        }
    }
    pub fn futures_price(&self, maturity: T) -> T {
        let t = maturity.max_val(T::zero());
        let b = self.b(t);
        let nodes = [
            -0.9602898565,
            -0.7966664774,
            -0.5255324099,
            -0.1834346425,
            0.1834346425,
            0.5255324099,
            0.7966664774,
            0.9602898565,
        ];
        let weights = [
            0.1012285363,
            0.2223810345,
            0.3137066459,
            0.3626837834,
            0.3626837834,
            0.3137066459,
            0.2223810345,
            0.1012285363,
        ];
        let mut var = T::zero();
        for (&x, &w) in nodes.iter().zip(weights.iter()) {
            let u = t * (x + 1.0) * 0.5;
            let bu = self.b(u);
            var += (self.sigma_cy * self.sigma_cy * bu * bu
                - self.rho * self.sigma_spot * self.sigma_cy * bu * 2.0)
                * (t * 0.5 * w);
        }
        self.spot
            * (-b * self.convenience_yield0
                + (self.rate - self.alpha_rn) * t
                + b * self.alpha_rn
                + var * 0.5)
                .exp()
    }
    pub fn log_futures_variance(&self, t_opt: T, t_fut: T) -> T {
        let t = t_opt.max_val(T::zero());
        let nodes = [
            -0.9602898565,
            -0.7966664774,
            -0.5255324099,
            -0.1834346425,
            0.1834346425,
            0.5255324099,
            0.7966664774,
            0.9602898565,
        ];
        let weights = [
            0.1012285363,
            0.2223810345,
            0.3137066459,
            0.3626837834,
            0.3626837834,
            0.3137066459,
            0.2223810345,
            0.1012285363,
        ];
        let mut v = T::zero();
        for (&x, &w) in nodes.iter().zip(weights.iter()) {
            let u = t * (x + 1.0) * 0.5;
            let b = self.b(t_fut - u);
            v += (self.sigma_spot * self.sigma_spot
                - self.rho * self.sigma_spot * self.sigma_cy * b * 2.0
                + self.sigma_cy * self.sigma_cy * b * b)
                * (t * 0.5 * w);
        }
        v.max_val(T::zero())
    }
    pub fn price_futures_option(
        &self,
        t_opt: T,
        t_fut: T,
        strike: T,
        option_type: ql_instruments::OptionType,
    ) -> T {
        let var = self.log_futures_variance(t_opt, t_fut);
        let vol = if t_opt.val() > 1e-14 {
            (var / t_opt).sqrt()
        } else {
            T::zero()
        };
        ql_vol::black76_price(
            self.futures_price(t_fut),
            strike,
            vol,
            t_opt,
            (-self.rate * t_opt).exp(),
            option_type.is_call(),
        )
    }
}

#[derive(Clone, Debug)]
pub struct ClewlowStricklandSwaptionSpec<T: Real = f64> {
    pub option_expiry: T,
    pub delivery_times: Vec<f64>,
    pub forward_prices: Vec<T>,
    pub delivery_weights: Vec<f64>,
    pub strike: T,
    pub rate: T,
    pub sigma_short: T,
    pub alpha_short: T,
    pub sigma_long: T,
    pub rho_factors: T,
    pub seasonal_vol_amp: T,
    pub is_payer_call: bool,
}

#[derive(Copy, Clone, Debug)]
pub struct ClewlowStricklandValuation<T: Real> {
    pub npv: T,
    pub forward_swap_price: T,
    pub effective_strip_vol: T,
    pub first_contract_samuelson_vol: T,
    pub last_contract_samuelson_vol: T,
}

pub fn price_clewlow_strickland_strip_swaption<T: Real>(
    spec: &ClewlowStricklandSwaptionSpec<T>,
) -> ClewlowStricklandValuation<T> {
    let m = spec
        .delivery_times
        .len()
        .min(spec.forward_prices.len())
        .min(spec.delivery_weights.len());
    let te = spec.option_expiry.max_val(T::from_f64(1e-8));
    if m == 0 {
        return ClewlowStricklandValuation {
            npv: T::zero(),
            forward_swap_price: T::zero(),
            effective_strip_vol: T::zero(),
            first_contract_samuelson_vol: T::zero(),
            last_contract_samuelson_vol: T::zero(),
        };
    }

    let w_sum: f64 = spec.delivery_weights[..m].iter().sum::<f64>().max(1e-12);
    let mut w_norm = vec![0.0; m];
    let mut f_strip = T::zero();
    for i in 0..m {
        w_norm[i] = spec.delivery_weights[i] / w_sum;
        f_strip += spec.forward_prices[i] * w_norm[i];
    }

    let a = spec.alpha_short.max_val(T::from_f64(1e-6));
    let int_ss = ((a * te * 2.0).exp() - T::one()) / (a * 2.0);
    let int_sl = ((a * te).exp() - T::one()) / a;

    let contract_cov = |ti: f64, tj: f64| -> T {
        let t_i = ti.max(te.val());
        let t_j = tj.max(te.val());
        let s_i = T::one() + spec.seasonal_vol_amp * (2.0 * std::f64::consts::PI * t_i).cos();
        let s_j = T::one() + spec.seasonal_vol_amp * (2.0 * std::f64::consts::PI * t_j).cos();
        let e_i = (-a * t_i).exp();
        let e_j = (-a * t_j).exp();
        s_i * s_j
            * (spec.sigma_short * spec.sigma_short * e_i * e_j * int_ss
                + spec.rho_factors * spec.sigma_short * spec.sigma_long * (e_i + e_j) * int_sl
                + spec.sigma_long * spec.sigma_long * te)
    };

    let mut m2 = T::zero();
    for i in 0..m {
        for j in 0..m {
            let cov_ij = contract_cov(spec.delivery_times[i], spec.delivery_times[j]);
            m2 += spec.forward_prices[i]
                * spec.forward_prices[j]
                * (w_norm[i] * w_norm[j])
                * cov_ij.exp();
        }
    }
    let f_sq = (f_strip * f_strip).max_val(T::from_f64(1e-16));
    let strip_var = (m2 / f_sq).max_val(T::from_f64(1.0 + 1e-12)).ln();
    let eff_vol = (strip_var / te).sqrt();
    let v_first = (contract_cov(spec.delivery_times[0], spec.delivery_times[0]) / te)
        .max_val(T::zero())
        .sqrt();
    let v_last = (contract_cov(spec.delivery_times[m - 1], spec.delivery_times[m - 1]) / te)
        .max_val(T::zero())
        .sqrt();

    let df = (-spec.rate * te).exp();
    let npv = black76_price(f_strip, spec.strike, eff_vol, te, df, spec.is_payer_call);

    ClewlowStricklandValuation {
        npv,
        forward_swap_price: f_strip,
        effective_strip_vol: eff_vol,
        first_contract_samuelson_vol: v_first,
        last_contract_samuelson_vol: v_last,
    }
}

#[cfg(test)]
mod phase15_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn gibson_parity() {
        let m = GibsonSchwartz2F {
            spot: 100.0,
            convenience_yield0: 0.03,
            rate: 0.02,
            kappa: 1.0,
            alpha_rn: 0.02,
            sigma_spot: 0.2,
            sigma_cy: 0.3,
            rho: -0.2,
        };
        let c = m.price_futures_option(1.0, 2.0, 100.0, ql_instruments::OptionType::Call);
        let p = m.price_futures_option(1.0, 2.0, 100.0, ql_instruments::OptionType::Put);
        assert!((c - p - (-0.02f64).exp() * (m.futures_price(2.0) - 100.0)).abs() < 1e-10);
        let h = 1e-5;
        let mut m_up = m;
        let mut m_dn = m;
        m_up.sigma_cy += h;
        m_dn.sigma_cy -= h;
        let fd = (m_up.futures_price(2.0) - m_dn.futures_price(2.0)) / (2.0 * h);
        let md = GibsonSchwartz2F::<Dual<1>> {
            spot: Dual::from_value(100.0),
            convenience_yield0: Dual::from_value(0.03),
            rate: Dual::from_value(0.02),
            kappa: Dual::from_value(1.0),
            alpha_rn: Dual::from_value(0.02),
            sigma_spot: Dual::from_value(0.2),
            sigma_cy: Dual::variable(0.3, 0),
            rho: Dual::from_value(-0.2),
        };
        let ad = md.futures_price(Dual::from_value(2.0));
        assert!((ad.d[0] - fd).abs() < 1e-6);
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

#[cfg(test)]
mod phase19_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn clewlow_strickland_samuelson_effect_strip_parity_and_ad() {
        let spec = ClewlowStricklandSwaptionSpec {
            option_expiry: 0.5,
            delivery_times: vec![0.75, 1.0, 1.25, 1.5],
            forward_prices: vec![75.0, 76.0, 74.0, 73.0],
            delivery_weights: vec![0.25; 4],
            strike: 74.5,
            rate: 0.03,
            sigma_short: 0.35,
            alpha_short: 1.2,
            sigma_long: 0.15,
            rho_factors: 0.2,
            seasonal_vol_amp: 0.0,
            is_payer_call: true,
        };
        let call = price_clewlow_strickland_strip_swaption(&spec);
        // Samuelson attenuation: front contract volatility > back contract volatility when alpha_short > 0
        assert!(call.first_contract_samuelson_vol > call.last_contract_samuelson_vol);

        let mut put_spec = spec.clone();
        put_spec.is_payer_call = false;
        let put = price_clewlow_strickland_strip_swaption(&put_spec);
        let df = (-0.03 * 0.5_f64).exp();
        assert!((call.npv - put.npv - df * (call.forward_swap_price - 74.5)).abs() < 1e-11);

        let h = 1e-5;
        let eval_f0 = |f0: f64| {
            let mut s = spec.clone();
            s.forward_prices[0] = f0;
            price_clewlow_strickland_strip_swaption(&s).npv
        };
        let fd = (eval_f0(75.0 + h) - eval_f0(75.0 - h)) / (2.0 * h);
        type D = Dual<1>;
        let spec_d = ClewlowStricklandSwaptionSpec {
            option_expiry: D::from_value(0.5),
            delivery_times: vec![0.75, 1.0, 1.25, 1.5],
            forward_prices: vec![
                D::variable(75.0, 0),
                D::from_value(76.0),
                D::from_value(74.0),
                D::from_value(73.0),
            ],
            delivery_weights: vec![0.25; 4],
            strike: D::from_value(74.5),
            rate: D::from_value(0.03),
            sigma_short: D::from_value(0.35),
            alpha_short: D::from_value(1.2),
            sigma_long: D::from_value(0.15),
            rho_factors: D::from_value(0.2),
            seasonal_vol_amp: D::zero(),
            is_payer_call: true,
        };
        let ad = price_clewlow_strickland_strip_swaption(&spec_d).npv;
        assert!(ad.d[0] > 0.0);
        assert!((ad.d[0] - fd).abs() < 1e-6);
    }
}
