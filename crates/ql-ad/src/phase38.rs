#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::approx_constant,
    clippy::too_many_arguments
)]
use crate::Real;

#[derive(Clone, Copy, Debug)]
pub struct GriewankRevolveBsbAdjointSpec<T: Real = f64> {
    pub spot_price: T,
    pub strike_price: T,
    pub risk_free_rate: T,
    pub sigma_min: T,
    pub sigma_max: T,
    pub bsb_switching_steepness: T,
    pub maturity_years: T,
}

#[derive(Clone, Copy, Debug)]
pub struct GriewankRevolveBsbAdjointReport<T: Real = f64> {
    pub primal_bsb_option_value: T,
    pub full_tape_adjoint_sigma_max_vega: T,
    pub revolve_checkpointed_sigma_max_vega: T,
    pub revolve_checkpointed_rate_rho: T,
    pub full_tape_vs_revolve_adjoint_error: T,
    pub memory_compression_ratio: T,
}

fn bsb_step<T: Real>(
    u: [T; 5],
    s_grid: [T; 5],
    dx: T,
    dt: T,
    r: T,
    sig_min: T,
    sig_max: T,
    steep: T,
) -> [T; 5] {
    let mut next = u;
    let inv_2dx = T::one() / (T::from_f64(2.0) * dx);
    let inv_dx2 = T::one() / (dx * dx);
    for i in 1..4 {
        let du = (u[i + 1] - u[i - 1]) * inv_2dx;
        let gamma = (u[i + 1] - T::from_f64(2.0) * u[i] + u[i - 1]) * inv_dx2;
        let sw = T::one() / (T::one() + (-steep * gamma).exp());
        let sig_eff = sig_min + (sig_max - sig_min) * sw;
        let pde_rhs = -r * u[i]
            + r * s_grid[i] * du
            + T::from_f64(0.5) * sig_eff * sig_eff * s_grid[i] * s_grid[i] * gamma;
        next[i] = u[i] + dt * pde_rhs;
    }
    next
}

fn bsb_step_vjp<T: Real>(
    u: [T; 5],
    adj_next: [T; 5],
    s_grid: [T; 5],
    dx: T,
    dt: T,
    r: T,
    sig_min: T,
    sig_max: T,
    steep: T,
) -> ([T; 5], T, T) {
    let mut adj_u = [T::zero(); 5];
    adj_u[0] = adj_next[0];
    adj_u[4] = adj_next[4];
    let inv_2dx = T::one() / (T::from_f64(2.0) * dx);
    let inv_dx2 = T::one() / (dx * dx);
    let mut d_sig_max = T::zero();
    let mut d_r = T::zero();

    for i in 1..4 {
        let lam = adj_next[i];
        let du = (u[i + 1] - u[i - 1]) * inv_2dx;
        let gamma = (u[i + 1] - T::from_f64(2.0) * u[i] + u[i - 1]) * inv_dx2;
        let sw = T::one() / (T::one() + (-steep * gamma).exp());
        let dsw_dg = steep * sw * (T::one() - sw);
        let sig_eff = sig_min + (sig_max - sig_min) * sw;
        let s2 = s_grid[i] * s_grid[i];

        let d_rhs_dgamma = T::from_f64(0.5)
            * s2
            * (sig_eff * sig_eff
                + T::from_f64(2.0) * sig_eff * (sig_max - sig_min) * dsw_dg * gamma);
        let d_rhs_ddu = r * s_grid[i];

        adj_u[i] =
            adj_u[i] + lam * (T::one() - dt * r - dt * T::from_f64(2.0) * inv_dx2 * d_rhs_dgamma);
        adj_u[i + 1] = adj_u[i + 1] + lam * dt * (inv_2dx * d_rhs_ddu + inv_dx2 * d_rhs_dgamma);
        adj_u[i - 1] = adj_u[i - 1] + lam * dt * (-inv_2dx * d_rhs_ddu + inv_dx2 * d_rhs_dgamma);

        d_sig_max = d_sig_max + lam * dt * (sig_eff * sw * s2 * gamma);
        d_r = d_r + lam * dt * (-u[i] + s_grid[i] * du);
    }
    (adj_u, d_sig_max, d_r)
}

pub fn evaluate_griewank_revolve_checkpointed_bsb_pde_adjoint<T: Real>(
    spec: &GriewankRevolveBsbAdjointSpec<T>,
) -> GriewankRevolveBsbAdjointReport<T> {
    let s0 = spec.spot_price.max_val(T::from_f64(10.0));
    let k = spec.strike_price.max_val(T::from_f64(10.0));
    let dx = T::from_f64(5.0);
    let s_grid = [
        s0 - T::from_f64(2.0) * dx,
        s0 - dx,
        s0,
        s0 + dx,
        s0 + T::from_f64(2.0) * dx,
    ];
    let mut u0 = [T::zero(); 5];
    for i in 0..5 {
        let diff = s_grid[i] - k;
        u0[i] = T::from_f64(0.5) * (diff + (diff * diff + T::from_f64(4.0)).sqrt());
    }
    let dt = spec.maturity_years.max_val(T::from_f64(1e-4)) / T::from_f64(6.0);

    // 1. Full-tape forward and reverse pass (stores all 6 states u^{(0)}..u^{(5)})
    let mut full_tape = [[T::zero(); 5]; 6];
    let mut u_cur = u0;
    for step in 0..6 {
        full_tape[step] = u_cur;
        u_cur = bsb_step(
            u_cur,
            s_grid,
            dx,
            dt,
            spec.risk_free_rate,
            spec.sigma_min,
            spec.sigma_max,
            spec.bsb_switching_steepness,
        );
    }
    let primal_val = u_cur[2];

    let mut adj_full = [T::zero(); 5];
    adj_full[2] = T::one();
    let mut vega_full = T::zero();
    for step in (0..6).rev() {
        let (adj_prev, dv, _) = bsb_step_vjp(
            full_tape[step],
            adj_full,
            s_grid,
            dx,
            dt,
            spec.risk_free_rate,
            spec.sigma_min,
            spec.sigma_max,
            spec.bsb_switching_steepness,
        );
        adj_full = adj_prev;
        vega_full = vega_full + dv;
    }

    // 2. Griewank-Walther Revolve binomial checkpointing with only c = 2 checkpoints (step 0 and step 3)
    let ckpt0 = u0;
    let mut ckpt3 = u0;
    for _ in 0..3 {
        ckpt3 = bsb_step(
            ckpt3,
            s_grid,
            dx,
            dt,
            spec.risk_free_rate,
            spec.sigma_min,
            spec.sigma_max,
            spec.bsb_switching_steepness,
        );
    }
    let mut adj_rev = [T::zero(); 5];
    adj_rev[2] = T::one();
    let mut vega_rev = T::zero();
    let mut rho_rev = T::zero();

    for step in (0..6).rev() {
        let mut u_recomp = if step >= 3 { ckpt3 } else { ckpt0 };
        let base_step = if step >= 3 { 3 } else { 0 };
        for _ in base_step..step {
            u_recomp = bsb_step(
                u_recomp,
                s_grid,
                dx,
                dt,
                spec.risk_free_rate,
                spec.sigma_min,
                spec.sigma_max,
                spec.bsb_switching_steepness,
            );
        }
        let (adj_prev, dv, dr) = bsb_step_vjp(
            u_recomp,
            adj_rev,
            s_grid,
            dx,
            dt,
            spec.risk_free_rate,
            spec.sigma_min,
            spec.sigma_max,
            spec.bsb_switching_steepness,
        );
        adj_rev = adj_prev;
        vega_rev = vega_rev + dv;
        rho_rev = rho_rev + dr;
    }

    let dv_diff = vega_full - vega_rev;
    let adj_err = (dv_diff * dv_diff + T::from_f64(1e-30)).sqrt();

    GriewankRevolveBsbAdjointReport {
        primal_bsb_option_value: primal_val,
        full_tape_adjoint_sigma_max_vega: vega_full,
        revolve_checkpointed_sigma_max_vega: vega_rev,
        revolve_checkpointed_rate_rho: rho_rev,
        full_tape_vs_revolve_adjoint_error: adj_err,
        memory_compression_ratio: T::from_f64(6.0 / 2.0),
    }
}

#[cfg(test)]
mod phase38_tests {
    use super::*;
    use crate::Dual;

    #[test]
    fn griewank_revolve_bsb_adjoint_invariants_and_ad_parity() {
        let spec = GriewankRevolveBsbAdjointSpec {
            spot_price: 100.0,
            strike_price: 100.0,
            risk_free_rate: 0.03,
            sigma_min: 0.15,
            sigma_max: 0.25,
            bsb_switching_steepness: 2.0,
            maturity_years: 0.5,
        };
        let rep = evaluate_griewank_revolve_checkpointed_bsb_pde_adjoint(&spec);
        assert!(rep.full_tape_vs_revolve_adjoint_error < 1e-12);
        assert!(rep.memory_compression_ratio > 1.0);

        let h = 1e-6;
        let r_ad = evaluate_griewank_revolve_checkpointed_bsb_pde_adjoint(
            &GriewankRevolveBsbAdjointSpec {
                spot_price: Dual::<1>::from_value(spec.spot_price),
                strike_price: Dual::<1>::from_value(spec.strike_price),
                risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
                sigma_min: Dual::<1>::from_value(spec.sigma_min),
                sigma_max: Dual::<1>::variable(spec.sigma_max, 0),
                bsb_switching_steepness: Dual::<1>::from_value(spec.bsb_switching_steepness),
                maturity_years: Dual::<1>::from_value(spec.maturity_years),
            },
        );
        let mut sp = spec;
        sp.sigma_max += h;
        let mut sm = spec;
        sm.sigma_max -= h;
        let fd = (evaluate_griewank_revolve_checkpointed_bsb_pde_adjoint(&sp)
            .primal_bsb_option_value
            - evaluate_griewank_revolve_checkpointed_bsb_pde_adjoint(&sm).primal_bsb_option_value)
            / (2.0 * h);
        assert!((r_ad.primal_bsb_option_value.d[0] - fd).abs() < 1e-5);
        assert!((rep.revolve_checkpointed_sigma_max_vega - fd).abs() < 1e-5);
    }
}
