#![allow(clippy::excessive_precision, clippy::needless_range_loop)]
use ql_ad::Real;
use ql_math::bivariate_normal_cdf;
use ql_vol::black_scholes_price;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum RainbowKind {
    CallOnMin,
    CallOnMax,
    PutOnMin,
    PutOnMax,
}

pub struct StulzRainbowEngine;
impl StulzRainbowEngine {
    #[allow(clippy::too_many_arguments)]
    pub fn price<T: Real>(
        s1: T,
        s2: T,
        k: T,
        r: T,
        q1: T,
        q2: T,
        v1: T,
        v2: T,
        rho: T,
        t: T,
        kind: RainbowKind,
    ) -> T {
        if k.val() <= 0.0 {
            return match kind {
                RainbowKind::CallOnMin => Self::expected_min(s1, s2, q1, q2, v1, v2, rho, t),
                RainbowKind::CallOnMax => {
                    s1 * (-q1 * t).exp() + s2 * (-q2 * t).exp()
                        - Self::expected_min(s1, s2, q1, q2, v1, v2, rho, t)
                }
                RainbowKind::PutOnMin | RainbowKind::PutOnMax => T::zero(),
            };
        }
        // Stulz (1982) image formula for a call on the minimum.
        let variance = (v1 * v1 + v2 * v2 - rho * v1 * v2 * 2.0).max_val(T::from_f64(1e-16));
        let sig = variance.sqrt();
        let sig_t = (variance * t).sqrt().max_val(T::from_f64(1e-8));
        let d = ((s1 / s2).ln() + (q2 - q1 + variance * 0.5) * t) / sig_t;
        let st1 = v1 * t.sqrt();
        let st2 = v2 * t.sqrt();
        let y1 = ((s1 / k).ln() + (r - q1 + v1 * v1 * 0.5) * t) / st1;
        let y2 = ((s2 / k).ln() + (r - q2 + v2 * v2 * 0.5) * t) / st2;
        let rho1 = (v1 - rho * v2) / sig;
        let rho2 = (v2 - rho * v1) / sig;
        let cmin = s1 * (-q1 * t).exp() * bivariate_normal_cdf(y1, -d, -rho1)
            + s2 * (-q2 * t).exp() * bivariate_normal_cdf(y2, d - sig_t, -rho2)
            - k * (-r * t).exp() * bivariate_normal_cdf(y1 - st1, y2 - st2, rho);
        let c1 = black_scholes_price(s1, k, v1, t, r, q1, true);
        let c2 = black_scholes_price(s2, k, v2, t, r, q2, true);
        let cmax = c1 + c2 - cmin;
        let emin = Self::expected_min(s1, s2, q1, q2, v1, v2, rho, t);
        let emax = s1 * (-q1 * t).exp() + s2 * (-q2 * t).exp() - emin;
        match kind {
            RainbowKind::CallOnMin => cmin,
            RainbowKind::CallOnMax => cmax,
            RainbowKind::PutOnMin => cmin - emin + k * (-r * t).exp(),
            RainbowKind::PutOnMax => cmax - emax + k * (-r * t).exp(),
        }
    }
    #[allow(clippy::too_many_arguments)]
    fn expected_min<T: Real>(s1: T, s2: T, q1: T, q2: T, v1: T, v2: T, rho: T, t: T) -> T {
        let var = (v1 * v1 + v2 * v2 - rho * 2.0 * v1 * v2).max_val(T::from_f64(1e-14));
        let st = var.sqrt() * t.sqrt();
        if st.val() < 1e-8 {
            return (s1 * (-q1 * t).exp()).min_val(s2 * (-q2 * t).exp());
        }
        let d = ((s1 / s2).ln() + (q2 - q1 + var * 0.5) * t) / st;
        s1 * (-q1 * t).exp() * (-d).norm_cdf() + s2 * (-q2 * t).exp() * (d - st).norm_cdf()
    }
}

pub struct TwoAssetCorrelationEngine;
impl TwoAssetCorrelationEngine {
    #[allow(clippy::too_many_arguments)]
    pub fn price<T: Real>(
        s1: T,
        s2: T,
        k1: T,
        k2: T,
        r: T,
        q1: T,
        q2: T,
        v1: T,
        v2: T,
        rho: T,
        t: T,
        is_call: bool,
    ) -> T {
        let st1 = v1 * t.sqrt();
        let st2 = v2 * t.sqrt();
        let y1 = ((s1 / k1).ln() + (r - q1 - v1 * v1 * 0.5) * t) / st1;
        let y2 = ((s2 / k2).ln() + (r - q2 - v2 * v2 * 0.5) * t) / st2;
        if is_call {
            (s2 * (-q2 * t).exp() * bivariate_normal_cdf(y2 + st2, y1 + rho * st2, rho)
                - k2 * (-r * t).exp() * bivariate_normal_cdf(y2, y1, rho))
            .max_val(T::zero())
        } else {
            (k2 * (-r * t).exp() * bivariate_normal_cdf(-y2, -y1, rho)
                - s2 * (-q2 * t).exp() * bivariate_normal_cdf(-y2 - st2, -y1 - rho * st2, rho))
            .max_val(T::zero())
        }
    }
}

pub struct DengLiZhouSpreadBasketEngine;
impl DengLiZhouSpreadBasketEngine {
    #[allow(clippy::too_many_arguments)]
    pub fn price_spread<T: Real>(
        s1: T,
        s2: T,
        k: T,
        r: T,
        q1: T,
        q2: T,
        v1: T,
        v2: T,
        rho: T,
        t: T,
        is_call: bool,
    ) -> T {
        if k.val().abs() < 1e-14 {
            let vv = (v1 * v1 + v2 * v2 - rho * 2.0 * v1 * v2).max_val(T::from_f64(1e-14));
            let st = vv.sqrt() * t.sqrt();
            let d1 = ((s1 / s2).ln() + (q2 - q1 + vv * 0.5) * t) / st;
            let d2 = d1 - st;
            let call = s1 * (-q1 * t).exp() * d1.norm_cdf() - s2 * (-q2 * t).exp() * d2.norm_cdf();
            return if is_call {
                call
            } else {
                call - s1 * (-q1 * t).exp() + s2 * (-q2 * t).exp()
            };
        }
        // Li--Deng--Zhou boundary expansion in its equivalent conditional
        // one-dimensional lognormal representation.  Conditioning on the
        // second asset retains the full correlation rather than using Kirk's
        // first-order variance proxy.
        const X: [f64; 8] = [
            -2.930_637_420_257_244,
            -1.981_656_756_695_843,
            -1.157_193_712_446_78,
            -0.381_186_990_207_322_1,
            0.381_186_990_207_322_1,
            1.157_193_712_446_78,
            1.981_656_756_695_843,
            2.930_637_420_257_244,
        ];
        const W: [f64; 8] = [
            0.000_199_604_072_211_367_6,
            0.017_077_983_007_413_48,
            0.207_802_325_814_891_9,
            0.661_147_012_558_241_3,
            0.661_147_012_558_241_3,
            0.207_802_325_814_891_9,
            0.017_077_983_007_413_48,
            0.000_199_604_072_211_367_6,
        ];
        let st = t.sqrt();
        let var_cond = (v1 * v1 * (T::one() - rho * rho) * t).max_val(T::from_f64(1e-16));
        let vol_cond = var_cond.sqrt();
        let mut call = T::zero();
        for i in 0..X.len() {
            let z = T::from_f64(X[i] * 2.0_f64.sqrt());
            let s2_t = s2 * ((r - q2 - v2 * v2 * 0.5) * t + v2 * st * z).exp();
            let f1_cond = (s1.ln() + (r - q1 - v1 * v1 * 0.5) * t + rho * v1 * st * z).exp()
                * (var_cond * 0.5).exp();
            let strike_cond = s2_t + k;
            if strike_cond.val() <= 0.0 {
                call += T::from_f64(W[i] / std::f64::consts::PI.sqrt()) * (f1_cond - strike_cond);
                continue;
            }
            let d1 = (f1_cond / strike_cond).ln() / vol_cond + vol_cond * 0.5;
            let d2 = d1 - vol_cond;
            call += T::from_f64(W[i] / std::f64::consts::PI.sqrt())
                * (f1_cond * d1.norm_cdf() - strike_cond * d2.norm_cdf());
        }
        call *= (-r * t).exp();
        if is_call {
            call
        } else {
            call - s1 * (-q1 * t).exp() + s2 * (-q2 * t).exp() + k * (-r * t).exp()
        }
    }
}

pub struct LevyBasketEngine;
impl LevyBasketEngine {
    #[allow(clippy::too_many_arguments)]
    pub fn price_basket<T: Real>(
        spots: &[T],
        weights: &[T],
        vols: &[T],
        divs: &[T],
        corr: &[Vec<T>],
        strike: T,
        rate: T,
        expiry: T,
        is_call: bool,
    ) -> T {
        assert_eq!(spots.len(), weights.len());
        let mut forward = T::zero();
        let mut var = T::zero();
        for i in 0..spots.len() {
            forward += weights[i] * spots[i] * ((rate - divs[i]) * expiry).exp();
        }
        for i in 0..spots.len() {
            for j in 0..spots.len() {
                var += weights[i]
                    * spots[i]
                    * ((rate - divs[i]) * expiry).exp()
                    * weights[j]
                    * spots[j]
                    * ((rate - divs[j]) * expiry).exp()
                    * ((corr[i][j] * vols[i] * vols[j] * expiry).exp() - T::one());
            }
        }
        let ratio = (T::one()
            + var.max_val(T::zero()) / (forward * forward).max_val(T::from_f64(1e-24)))
        .ln();
        let vol = (ratio / expiry.max_val(T::from_f64(1e-12)))
            .max_val(T::zero())
            .sqrt();
        ql_vol::black76_price(
            forward,
            strike,
            vol,
            expiry,
            (-rate * expiry).exp(),
            is_call,
        )
    }
}


#[derive(Clone, Debug)]
pub struct ChoiSumOfBsmBasketSpec<T: Real = f64> {
    pub spots: Vec<T>,
    pub weights: Vec<T>,
    pub vols: Vec<T>,
    pub div_yields: Vec<T>,
    pub correlation_matrix: Vec<Vec<T>>,
    pub strike: T,
    pub risk_free_rate: T,
    pub expiry: T,
    pub is_call: bool,
}

#[derive(Clone, Debug)]
pub struct ChoiSumOfBsmValuation<T: Real = f64> {
    pub price: T,
    pub raw_quadrature_price: T,
    pub control_variate_correction: T,
    pub deltas: Vec<T>,
    pub first_factor_loadings: Vec<T>,
}

fn cholesky_lower_generic<T: Real>(a: &[Vec<T>]) -> Vec<Vec<T>> {
    let n = a.len();
    let mut l = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i][j];
            for k in 0..j {
                sum -= l[i][k] * l[j][k];
            }
            if i == j {
                l[i][i] = sum.max_val(T::from_f64(1e-14)).sqrt();
            } else {
                l[i][j] = sum / l[j][j].max_val(T::from_f64(1e-14));
            }
        }
    }
    l
}

/// Choi (2018) "Sum of All Black-Scholes-Merton Models" Householder-rotated factor quadrature
/// for multi-asset basket & spread options (`papers/37_Choi_2018_Sum_of_All_BSM_Models_Basket_Spread_Options.pdf`).
pub fn price_choi_sum_of_bsm_basket_option<T: Real>(
    spec: &ChoiSumOfBsmBasketSpec<T>,
) -> ChoiSumOfBsmValuation<T> {
    let n = spec.spots.len();
    if n == 0 {
        return ChoiSumOfBsmValuation {
            price: T::zero(),
            raw_quadrature_price: T::zero(),
            control_variate_correction: T::zero(),
            deltas: Vec::new(),
            first_factor_loadings: Vec::new(),
        };
    }
    let t = spec.expiry.max_val(T::from_f64(1e-6));
    let df = (-spec.risk_free_rate * t).exp();
    let forwards: Vec<T> = (0..n)
        .map(|k| spec.spots[k] * ((spec.risk_free_rate - spec.div_yields[k]) * t).exp())
        .collect();
    let mut sigma_mat = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            sigma_mat[i][j] = spec.correlation_matrix[i][j] * spec.vols[i] * spec.vols[j] * t;
        }
    }
    let chol = cholesky_lower_generic(&sigma_mat);

    // Steepest-ascent direction g_tilde_k = w_k * F_k (Choi Eqs. 22-24)
    let g_tilde: Vec<T> = (0..n).map(|k| spec.weights[k] * forwards[k]).collect();
    let g_norm = g_tilde
        .iter()
        .fold(T::zero(), |acc, &x| acc + x * x)
        .max_val(T::from_f64(1e-24))
        .sqrt();
    let g: Vec<T> = g_tilde.into_iter().map(|x| x / g_norm).collect();

    let mut sig_g = vec![T::zero(); n];
    let mut g_sig_g = T::zero();
    for i in 0..n {
        let mut s = T::zero();
        for j in 0..n {
            s += sigma_mat[i][j] * g[j];
        }
        sig_g[i] = s;
        g_sig_g += g[i] * s;
    }
    let denom_v1 = g_sig_g.max_val(T::from_f64(1e-20)).sqrt();
    let mut v_star1: Vec<T> = sig_g.into_iter().map(|s| s / denom_v1).collect();
    for k in 0..n {
        if (spec.weights[k] * v_star1[k]).val() <= 0.0 {
            let sign_w = if spec.weights[k].val() >= 0.0 {
                1.0
            } else {
                -1.0
            };
            v_star1[k] = sigma_mat[k][k].max_val(T::from_f64(1e-14)).sqrt() * (0.01 * sign_w);
        }
    }

    // Solve C * q_tilde = v_star1 by forward substitution
    let mut q_tilde = vec![T::zero(); n];
    for i in 0..n {
        let mut s = v_star1[i];
        for j in 0..i {
            s -= chol[i][j] * q_tilde[j];
        }
        q_tilde[i] = s / chol[i][i].max_val(T::from_f64(1e-14));
    }
    let q_norm = q_tilde
        .iter()
        .fold(T::zero(), |acc, &x| acc + x * x)
        .max_val(T::from_f64(1e-24))
        .sqrt();
    let q1: Vec<T> = q_tilde.into_iter().map(|x| x / q_norm).collect();

    // Householder reflection R = I - 2 u u^T mapping e_1 to q1 (or -e_1 if q1[0] > 0)
    let mut u_vec = q1.clone();
    u_vec[0] -= T::one();
    let u_norm2 = u_vec.iter().fold(T::zero(), |acc, &x| acc + x * x);
    let mut r_house = vec![vec![T::zero(); n]; n];
    if u_norm2.val() < 1e-20 {
        for i in 0..n {
            r_house[i][i] = T::one();
        }
    } else {
        for i in 0..n {
            for j in 0..n {
                let diag = if i == j { T::one() } else { T::zero() };
                r_house[i][j] = diag - (u_vec[i] * u_vec[j] * 2.0) / u_norm2;
            }
        }
    }
    let mut v_mat = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            let mut s = T::zero();
            for k in 0..n {
                s += chol[i][k] * r_house[k][j];
            }
            v_mat[i][j] = s;
        }
    }
    // Ensure w_0 * V_{0, 0} > 0
    if (spec.weights[0] * v_mat[0][0]).val() < 0.0 {
        for i in 0..n {
            v_mat[i][0] = -v_mat[i][0];
        }
    }

    // 5-point Gauss-Hermite nodes for remaining dimensions (or 3-point if n >= 4)
    const GH5_Z: [f64; 5] = [
        -2.856_970_013_872_805_6,
        -1.355_626_179_974_265_9,
        0.0,
        1.355_626_179_974_265_9,
        2.856_970_013_872_805_6,
    ];
    const GH5_H: [f64; 5] = [
        0.011_257_411_327_720_689,
        0.222_075_922_005_612_64,
        0.533_333_333_333_333_3,
        0.222_075_922_005_612_64,
        0.011_257_411_327_720_689,
    ];
    const GH3_Z: [f64; 3] = [-1.732_050_807_568_877_2, 0.0, 1.732_050_807_568_877_2];
    const GH3_H: [f64; 3] = [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0];

    let rem_dims = n.saturating_sub(1);
    let (nodes, weights): (&[f64], &[f64]) = if rem_dims <= 2 {
        (&GH5_Z, &GH5_H)
    } else {
        (&GH3_Z, &GH3_H)
    };
    let m_pts = nodes.len();
    let total_nodes = m_pts.pow(rem_dims as u32);

    let mut c_raw = T::zero();
    let mut p_raw = T::zero();
    let mut d_vec = vec![T::zero(); n];
    let mut f_bar = vec![T::zero(); n];

    for idx in 0..total_nodes {
        let mut rem = idx;
        let mut h_weight = 1.0_f64;
        let mut z_rem = vec![0.0_f64; n];
        for d in 1..n {
            let coord = rem % m_pts;
            rem /= m_pts;
            z_rem[d] = nodes[coord];
            h_weight *= weights[coord];
        }
        let mut fk = vec![T::zero(); n];
        for k in 0..n {
            let mut half_v2 = T::zero();
            let mut vz = T::zero();
            for j in 1..n {
                half_v2 += v_mat[k][j] * v_mat[k][j] * 0.5;
                vz += v_mat[k][j] * z_rem[j];
            }
            fk[k] = (-half_v2 + vz).exp();
            f_bar[k] += fk[k] * h_weight;
        }

        // Solve 1D monotonic root sum_k w_k F_k f_k exp(-0.5 V_{k0}^2 + V_{k0} z_1) = K
        let eval_g_f64 = |z1: f64| -> f64 {
            let mut s = -spec.strike.val();
            for k in 0..n {
                let vk0 = v_mat[k][0].val();
                s += spec.weights[k].val()
                    * forwards[k].val()
                    * fk[k].val()
                    * (-0.5 * vk0 * vk0 + vk0 * z1).exp();
            }
            s
        };
        let mut lo = -12.0_f64;
        let mut hi = 12.0_f64;
        for _ in 0..32 {
            let mid = 0.5 * (lo + hi);
            if eval_g_f64(mid) < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let mut z1_star = 0.5 * (lo + hi);
        for _ in 0..5 {
            let mut g_val = -spec.strike.val();
            let mut dg_val = 0.0_f64;
            for k in 0..n {
                let vk0 = v_mat[k][0].val();
                let term = spec.weights[k].val()
                    * forwards[k].val()
                    * fk[k].val()
                    * (-0.5 * vk0 * vk0 + vk0 * z1_star).exp();
                g_val += term;
                dg_val += term * vk0;
            }
            if dg_val.abs() > 1e-14 {
                z1_star = (z1_star - g_val / dg_val).clamp(-12.0, 12.0);
            }
        }
        // 1-step IFT lift to T: Real
        let z1_t0 = T::from_f64(z1_star);
        let mut g_t = -spec.strike;
        let mut dg_t = T::zero();
        for k in 0..n {
            let vk0 = v_mat[k][0];
            let term =
                spec.weights[k] * forwards[k] * fk[k] * (-vk0 * vk0 * 0.5 + vk0 * z1_t0).exp();
            g_t += term;
            dg_t += term * vk0;
        }
        let z1_ad = z1_t0 - g_t / dg_t.max_val(T::from_f64(1e-12));
        let dm = -z1_ad;

        let mut c_node = -spec.strike * dm.norm_cdf();
        let mut p_node = spec.strike * (-dm).norm_cdf();
        for k in 0..n {
            let vk0 = v_mat[k][0];
            let nd_pos = (dm + vk0).norm_cdf();
            let nd_neg = (-dm - vk0).norm_cdf();
            c_node += spec.weights[k] * forwards[k] * fk[k] * nd_pos;
            p_node -= spec.weights[k] * forwards[k] * fk[k] * nd_neg;
            d_vec[k] += spec.weights[k] * fk[k] * nd_pos * h_weight;
        }
        c_raw += c_node * h_weight;
        p_raw += p_node * h_weight;
    }

    let mut cv_call = T::zero();
    let mut cv_put = T::zero();
    for k in 0..n {
        cv_call += d_vec[k] * forwards[k] * (f_bar[k] - T::one());
        cv_put += (d_vec[k] - spec.weights[k]) * forwards[k] * (f_bar[k] - T::one());
    }
    let call_price = df * (c_raw - cv_call);
    let put_price = df * (p_raw - cv_put);
    let (price, raw_quadrature_price, control_variate_correction) = if spec.is_call {
        (call_price, df * c_raw, -df * cv_call)
    } else {
        (put_price, df * p_raw, -df * cv_put)
    };
    let deltas: Vec<T> = (0..n)
        .map(|k| {
            let df_q = (-spec.div_yields[k] * t).exp();
            if spec.is_call {
                df_q * d_vec[k]
            } else {
                df_q * (d_vec[k] - spec.weights[k])
            }
        })
        .collect();
    let first_factor_loadings: Vec<T> = (0..n).map(|k| v_mat[k][0]).collect();

    ChoiSumOfBsmValuation {
        price,
        raw_quadrature_price,
        control_variate_correction,
        deltas,
        first_factor_loadings,
    }
}

#[derive(Clone, Debug)]
pub struct LiDengZhouMultiSpreadSpec<T: Real = f64> {
    pub long_spot: T,
    pub long_div_yield: T,
    pub long_vol: T,
    pub short_spots: Vec<T>,
    pub short_weights: Vec<T>,
    pub short_div_yields: Vec<T>,
    pub short_vols: Vec<T>,
    /// Full (N+1) x (N+1) correlation matrix where index 0 is the long asset and 1..=N are the short legs.
    pub correlation_matrix: Vec<Vec<T>>,
    pub strike: T,
    pub risk_free_rate: T,
    pub expiry: T,
    pub is_call: bool,
}

#[derive(Clone, Debug)]
pub struct LiDengZhouMultiSpreadValuation<T: Real = f64> {
    pub price: T,
    pub first_order_price: T,
    pub second_order_correction: T,
    pub long_delta: T,
    pub short_deltas: Vec<T>,
    pub boundary_intercept_c: T,
}

fn invert_spd_small<T: Real>(a: &[Vec<T>]) -> Vec<Vec<T>> {
    let n = a.len();
    let mut aug = vec![vec![T::zero(); 2 * n]; n];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = a[i][j];
        }
        aug[i][n + i] = T::one();
    }
    for col in 0..n {
        let mut piv = col;
        let mut max_v = aug[col][col].val().abs();
        for r in (col + 1)..n {
            if aug[r][col].val().abs() > max_v {
                max_v = aug[r][col].val().abs();
                piv = r;
            }
        }
        if piv != col {
            aug.swap(piv, col);
        }
        let diag = aug[col][col];
        for c in col..(2 * n) {
            aug[col][c] /= diag;
        }
        for r in 0..n {
            if r == col {
                continue;
            }
            let factor = aug[r][col];
            for c in col..(2 * n) {
                let v = aug[col][c];
                aug[r][c] -= factor * v;
            }
        }
    }
    let mut inv = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            inv[i][j] = aug[i][n + j];
        }
    }
    inv
}

fn mat_mul_small<T: Real>(a: &[Vec<T>], b: &[Vec<T>]) -> Vec<Vec<T>> {
    let n = a.len();
    let mut out = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        for k in 0..n {
            let aik = a[i][k];
            for j in 0..n {
                out[i][j] += aik * b[k][j];
            }
        }
    }
    out
}

fn spd_sqrt_denman_beavers<T: Real>(a: &[Vec<T>]) -> Vec<Vec<T>> {
    let n = a.len();
    let mut y = a.to_vec();
    let mut z = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        z[i][i] = T::one();
    }
    for _ in 0..8 {
        let y_inv = invert_spd_small(&y);
        let z_inv = invert_spd_small(&z);
        for i in 0..n {
            for j in 0..n {
                y[i][j] = (y[i][j] + z_inv[i][j]) * 0.5;
                z[i][j] = (z[i][j] + y_inv[i][j]) * 0.5;
            }
        }
    }
    y
}

/// Li, Deng & Zhou (2008) Multi-Asset Spread Option ($N \ge 2$ short legs) Second-Order
/// Exercise Boundary Approximation (`papers/36_Li_Deng_Zhou_2008_Multi_Asset_Spread_Option_Pricing.pdf`, Props. 3-6).
pub fn price_li_deng_zhou_multi_spread_option<T: Real>(
    spec: &LiDengZhouMultiSpreadSpec<T>,
) -> LiDengZhouMultiSpreadValuation<T> {
    let n = spec.short_spots.len();
    let t = spec.expiry.max_val(T::from_f64(1e-6));
    let st = t.sqrt();
    let r = spec.risk_free_rate;

    let s0 = spec.long_spot;
    let q0 = spec.long_div_yield;
    let nu0 = (spec.long_vol * st).max_val(T::from_f64(1e-8));
    let mu0 = s0.ln() + (r - q0 - spec.long_vol * spec.long_vol * 0.5) * t;

    let s_short: Vec<T> = (0..n)
        .map(|k| spec.short_weights[k] * spec.short_spots[k])
        .collect();
    let nu_short: Vec<T> = (0..n)
        .map(|k| (spec.short_vols[k] * st).max_val(T::from_f64(1e-8)))
        .collect();
    let mu_short: Vec<T> = (0..n)
        .map(|k| {
            s_short[k].ln()
                + (r - spec.short_div_yields[k] - spec.short_vols[k] * spec.short_vols[k] * 0.5) * t
        })
        .collect();
    let exp_mu: Vec<T> = mu_short.iter().map(|&m| m.exp()).collect();
    let r_sum = exp_mu.iter().fold(T::zero(), |acc, &x| acc + x);
    let r_plus_k = (r_sum + spec.strike).max_val(T::from_f64(1e-12));

    let mut sig10 = vec![T::zero(); n];
    let mut sig11 = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        sig10[i] = spec.correlation_matrix[i + 1][0];
        for j in 0..n {
            sig11[i][j] = spec.correlation_matrix[i + 1][j + 1];
        }
    }
    let sig11_inv = invert_spd_small(&sig11);
    let sig11_sqrt = spd_sqrt_denman_beavers(&sig11);

    let mut z_vec = vec![T::zero(); n];
    let mut quad_corr = T::zero();
    for i in 0..n {
        let mut s = T::zero();
        for j in 0..n {
            s += sig11_inv[i][j] * sig10[j];
        }
        z_vec[i] = s;
        quad_corr += sig10[i] * s;
    }
    let sig_xy = (T::one() - quad_corr).max_val(T::from_f64(1e-8));
    let sqrt_sig_xy = sig_xy.sqrt();

    // Gradient & Hessian of the critical boundary at 0 (Eqs. 14-15 & 29)
    let grad_x: Vec<T> = (0..n)
        .map(|k| exp_mu[k] * nu_short[k] / (nu0 * r_plus_k))
        .collect();
    let mut hess_x = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            let mut h_ij =
                -nu_short[i] * nu_short[j] * exp_mu[i] * exp_mu[j] / (nu0 * r_plus_k * r_plus_k);
            if i == j {
                h_ij += nu_short[j] * nu_short[j] * exp_mu[j] / (nu0 * r_plus_k);
            }
            hess_x[i][j] = h_ij;
        }
    }

    let c_base = -(r_plus_k.ln() - mu0) / (nu0 * sqrt_sig_xy);
    let d_vec: Vec<T> = (0..n)
        .map(|i| (z_vec[i] - grad_x[i]) / sqrt_sig_xy)
        .collect();
    let mut e_mat = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            e_mat[i][j] = -hess_x[i][j] / (sqrt_sig_xy * 2.0);
        }
    }
    let f_mat = mat_mul_small(&mat_mul_small(&sig11_sqrt, &e_mat), &sig11_sqrt);
    let mut tr_f = T::zero();
    let mut tr_f2 = T::zero();
    for i in 0..n {
        tr_f += f_mat[i][i];
        for j in 0..n {
            tr_f2 += f_mat[i][j] * f_mat[i][j];
        }
    }

    let eval_j012 = |u: T, v: &[T]| -> (T, T) {
        let mut v_norm2 = T::zero();
        let mut fv = vec![T::zero(); n];
        for i in 0..n {
            v_norm2 += v[i] * v[i];
            let mut s = T::zero();
            for j in 0..n {
                s += f_mat[i][j] * v[j];
            }
            fv[i] = s;
        }
        let mut q1 = T::zero();
        let mut q2 = T::zero();
        for i in 0..n {
            q1 += v[i] * fv[i];
            q2 += fv[i] * fv[i];
        }
        let psi = T::one() / (T::one() + v_norm2);
        let sqrt_psi = psi.sqrt();
        let psi_32 = psi * sqrt_psi;
        let arg = u * sqrt_psi;
        let phi_arg = arg.norm_pdf();
        let j0 = arg.norm_cdf();
        let u2 = u * u;
        let j1 = psi_32 * (psi * u2 - T::one()) * q1 * phi_arg;
        let j2 = u
            * psi_32
            * phi_arg
            * (tr_f2 * 2.0 - (T::one() - tr_f) * (psi - psi * psi) * q1 * 4.0
                + psi
                    * psi
                    * (T::from_f64(9.0) + (T::from_f64(2.0) - u2 * 3.0) * psi
                        - u2 * (T::from_f64(4.0) - u2) * psi * psi)
                    * q1
                    * q1
                - psi * (T::from_f64(5.0) + (T::one() - u2 * 2.0) * psi) * q2 * 2.0);
        (j0, j0 + j1 - j2 * 0.5)
    };

    // Shifted arguments (c_0, d_0), (c_k, d_k), and (c_{N+1}, d_{N+1}) (Eqs. 23-28)
    let mut sig10_d = T::zero();
    let mut e_sig10 = vec![T::zero(); n];
    for i in 0..n {
        sig10_d += sig10[i] * d_vec[i];
        let mut s = T::zero();
        for j in 0..n {
            s += e_mat[i][j] * sig10[j];
        }
        e_sig10[i] = s;
    }
    let mut sig10_e_sig10 = T::zero();
    for i in 0..n {
        sig10_e_sig10 += sig10[i] * e_sig10[i];
    }
    let c0 = c_base + tr_f + nu0 * sqrt_sig_xy + nu0 * sig10_d + nu0 * nu0 * sig10_e_sig10;
    let d0_in: Vec<T> = (0..n).map(|i| d_vec[i] + e_sig10[i] * (nu0 * 2.0)).collect();
    let d0: Vec<T> = (0..n)
        .map(|i| (0..n).fold(T::zero(), |acc, j| acc + sig11_sqrt[i][j] * d0_in[j]))
        .collect();
    let (i0_first, i0_sec) = eval_j012(c0, &d0);

    let sig11_e_sig11 = mat_mul_small(&mat_mul_small(&sig11, &e_mat), &sig11);
    let e_sig11 = mat_mul_small(&e_mat, &sig11);
    let mut ik_first = Vec::with_capacity(n);
    let mut ik_sec = Vec::with_capacity(n);
    for k in 0..n {
        let mut sig11_d_k = T::zero();
        for j in 0..n {
            sig11_d_k += sig11[k][j] * d_vec[j];
        }
        let ck = c_base + tr_f + nu_short[k] * sig11_d_k + nu_short[k] * nu_short[k] * sig11_e_sig11[k][k];
        let dk_in: Vec<T> = (0..n)
            .map(|i| d_vec[i] + e_sig11[i][k] * (nu_short[k] * 2.0))
            .collect();
        let dk: Vec<T> = (0..n)
            .map(|i| (0..n).fold(T::zero(), |acc, j| acc + sig11_sqrt[i][j] * dk_in[j]))
            .collect();
        let (f_k, s_k) = eval_j012(ck, &dk);
        ik_first.push(f_k);
        ik_sec.push(s_k);
    }

    let c_last = c_base + tr_f;
    let d_last: Vec<T> = (0..n)
        .map(|i| (0..n).fold(T::zero(), |acc, j| acc + sig11_sqrt[i][j] * d_vec[j]))
        .collect();
    let (il_first, il_sec) = eval_j012(c_last, &d_last);

    let df_r = (-r * t).exp();
    let df_q0 = (-q0 * t).exp();
    let mut call_first = s0 * df_q0 * i0_first - spec.strike * df_r * il_first;
    let mut call_sec = s0 * df_q0 * i0_sec - spec.strike * df_r * il_sec;
    let mut short_deltas = Vec::with_capacity(n);
    for k in 0..n {
        let df_qk = (-spec.short_div_yields[k] * t).exp();
        call_first -= s_short[k] * df_qk * ik_first[k];
        call_sec -= s_short[k] * df_qk * ik_sec[k];
        short_deltas.push(if spec.is_call {
            -spec.short_weights[k] * df_qk * ik_sec[k]
        } else {
            spec.short_weights[k] * df_qk * (T::one() - ik_sec[k])
        });
    }
    let long_delta = if spec.is_call {
        df_q0 * i0_sec
    } else {
        df_q0 * (i0_sec - T::one())
    };
    let fwd_spread_pv = s0 * df_q0
        - (0..n).fold(T::zero(), |acc, k| {
            acc + s_short[k] * (-spec.short_div_yields[k] * t).exp()
        })
        - spec.strike * df_r;
    let price = if spec.is_call {
        call_sec.max_val(T::zero())
    } else {
        (call_sec - fwd_spread_pv).max_val(T::zero())
    };
    let first_order_price = if spec.is_call {
        call_first.max_val(T::zero())
    } else {
        (call_first - fwd_spread_pv).max_val(T::zero())
    };

    LiDengZhouMultiSpreadValuation {
        price,
        first_order_price,
        second_order_correction: price - first_order_price,
        long_delta,
        short_deltas,
        boundary_intercept_c: c_base,
    }
}



#[cfg(test)]
mod phase20_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn choi_householder_basket_and_li_deng_zhou_multispread_parity_and_ad() {
        // 1. Choi (2018) 2-asset spread option vs Deng-Li-Zhou 1D quadrature benchmark
        let choi_spread = ChoiSumOfBsmBasketSpec {
            spots: vec![100.0, 90.0],
            weights: vec![1.0, -1.0],
            vols: vec![0.25, 0.20],
            div_yields: vec![0.01, 0.02],
            correlation_matrix: vec![vec![1.0, 0.5], vec![0.5, 1.0]],
            strike: 5.0,
            risk_free_rate: 0.03,
            expiry: 1.0,
            is_call: true,
        };
        let v_choi = price_choi_sum_of_bsm_basket_option(&choi_spread);
        let v_dlz = DengLiZhouSpreadBasketEngine::price_spread(
            100.0, 90.0, 5.0, 0.03, 0.01, 0.02, 0.25, 0.20, 0.5, 1.0, true,
        );
        assert!((v_choi.price - v_dlz).abs() < 1e-3);

        // 2. Li-Deng-Zhou (2008) 3-asset spread option (S_0 - 0.5 S_1 - 0.5 S_2 - K)_+
        let ldz_spec = LiDengZhouMultiSpreadSpec {
            long_spot: 100.0,
            long_div_yield: 0.01,
            long_vol: 0.25,
            short_spots: vec![90.0, 85.0],
            short_weights: vec![0.5, 0.5],
            short_div_yields: vec![0.015, 0.02],
            short_vols: vec![0.20, 0.22],
            correlation_matrix: vec![
                vec![1.0, 0.5, 0.4],
                vec![0.5, 1.0, 0.6],
                vec![0.4, 0.6, 1.0],
            ],
            strike: 5.0,
            risk_free_rate: 0.03,
            expiry: 1.0,
            is_call: true,
        };
        let v_ldz = price_li_deng_zhou_multi_spread_option(&ldz_spec);
        let choi_3 = ChoiSumOfBsmBasketSpec {
            spots: vec![100.0, 90.0, 85.0],
            weights: vec![1.0, -0.5, -0.5],
            vols: vec![0.25, 0.20, 0.22],
            div_yields: vec![0.01, 0.015, 0.02],
            correlation_matrix: ldz_spec.correlation_matrix.clone(),
            strike: 5.0,
            risk_free_rate: 0.03,
            expiry: 1.0,
            is_call: true,
        };
        let v_choi3 = price_choi_sum_of_bsm_basket_option(&choi_3);
        assert!((v_ldz.price - v_choi3.price).abs() < 2e-2);
        assert!(v_ldz.long_delta > 0.0 && v_ldz.short_deltas[0] < 0.0);

        let h = 1e-5;
        let eval_s0 = |s0: f64| {
            let mut sp = ldz_spec.clone();
            sp.long_spot = s0;
            price_li_deng_zhou_multi_spread_option(&sp).price
        };
        let fd = (eval_s0(100.0 + h) - eval_s0(100.0 - h)) / (2.0 * h);
        type D = Dual<1>;
        let spec_d = LiDengZhouMultiSpreadSpec {
            long_spot: D::variable(100.0, 0),
            long_div_yield: D::from_value(0.01),
            long_vol: D::from_value(0.25),
            short_spots: vec![D::from_value(90.0), D::from_value(85.0)],
            short_weights: vec![D::from_value(0.5), D::from_value(0.5)],
            short_div_yields: vec![D::from_value(0.015), D::from_value(0.02)],
            short_vols: vec![D::from_value(0.20), D::from_value(0.22)],
            correlation_matrix: vec![
                vec![D::one(), D::from_value(0.5), D::from_value(0.4)],
                vec![D::from_value(0.5), D::one(), D::from_value(0.6)],
                vec![D::from_value(0.4), D::from_value(0.6), D::one()],
            ],
            strike: D::from_value(5.0),
            risk_free_rate: D::from_value(0.03),
            expiry: D::one(),
            is_call: true,
        };
        let ad = price_li_deng_zhou_multi_spread_option(&spec_d).price;
        assert!((ad.d[0] - fd).abs() < 1e-6);
    }
}
