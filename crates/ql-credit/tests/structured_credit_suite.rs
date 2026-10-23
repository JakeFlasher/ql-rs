use ql_ad::Dual;
use ql_credit::*;

#[test]
fn asb_conservation_and_dual_sensitivity() {
    let probs = vec![0.02; 20];
    let units = vec![1usize; 20];
    let tranches = [(0.0, 0.03), (0.03, 0.07), (0.07, 0.15), (0.15, 1.0)];
    let total: f64 = tranches
        .iter()
        .map(|&(a, d)| asb_recursive_expected_tranche_loss(&probs, &units, 1.0 / 20.0, 0.2, a, d))
        .sum();
    assert!((total - 0.02).abs() < 1e-12, "{total}");
    let total_high_corr: f64 = tranches
        .iter()
        .map(|&(a, d)| asb_recursive_expected_tranche_loss(&probs, &units, 1.0 / 20.0, 0.95, a, d))
        .sum();
    assert!((total_high_corr - 0.02).abs() < 1e-12);

    let p: Vec<Dual<5>> = probs
        .iter()
        .enumerate()
        .map(|(i, &q)| {
            if i < 4 {
                Dual::variable(q, i)
            } else {
                Dual::from_value(q)
            }
        })
        .collect();
    let d = asb_recursive_expected_tranche_loss(
        &p,
        &units,
        Dual::from_value(1.0 / 20.0),
        Dual::variable(0.2, 4),
        Dual::from_value(0.03),
        Dual::from_value(0.07),
    );
    let h = 1e-5;
    for j in 0..4 {
        let mut up = probs.clone();
        let mut dn = probs.clone();
        up[j] += h;
        dn[j] -= h;
        let fd = (asb_recursive_expected_tranche_loss(&up, &units, 1.0 / 20.0, 0.2, 0.03, 0.07)
            - asb_recursive_expected_tranche_loss(&dn, &units, 1.0 / 20.0, 0.2, 0.03, 0.07))
            / (2.0 * h);
        assert!((d.d[j] - fd).abs() < 1e-6);
    }
}

#[test]
fn lhp_and_cds_option_front_end_protection() {
    let l = gaussian_lhp_expected_tranche_loss(0.02, 0.4, 0.2, 0.03, 0.07);
    assert!((0.0..=0.04).contains(&l));
    let no_fep = black_cds_option_price(
        0.02,
        0.015,
        4.0,
        0.3,
        1.0,
        0.02,
        0.97,
        0.4,
        1_000_000.0,
        true,
        true,
    );
    let with_fep = black_cds_option_price(
        0.02,
        0.015,
        4.0,
        0.3,
        1.0,
        0.02,
        0.97,
        0.4,
        1_000_000.0,
        true,
        false,
    );
    assert!(with_fep > no_fep);
    let receiver = black_cds_option_price(
        0.02,
        0.015,
        4.0,
        0.3,
        1.0,
        0.02,
        0.97,
        0.4,
        1_000_000.0,
        false,
        true,
    );
    assert!(receiver >= 0.0);
}

#[test]
fn cdo_tranche_spread_seniority_and_lhp_convergence() {
    let probs = vec![0.02; 20];
    let units = vec![1usize; 20];
    let make = |attachment, detachment| SyntheticCdoTranche {
        attachment,
        detachment,
        running_spread: 0.01,
        payment_times: vec![1.0],
        accruals: vec![1.0],
    };
    let fair = |a, d| price_cdo_tranche(&make(a, d), &probs, &units, 0.006, 0.2, &[1.0]).2;
    let fair_eq = fair(0.0, 0.03);
    let fair_mezz = fair(0.03, 0.07);
    let fair_sen = fair(0.07, 0.15);
    assert!(fair_eq > fair_mezz && fair_mezz > fair_sen && fair_sen > 0.0);

    let homogeneous_probs = vec![0.02; 100];
    let homogeneous_units = vec![1usize; 100];
    let asb = asb_recursive_expected_tranche_loss(
        &homogeneous_probs,
        &homogeneous_units,
        0.006,
        0.2,
        0.03,
        0.07,
    );
    let lhp = gaussian_lhp_expected_tranche_loss(0.02, 0.4, 0.2, 0.03, 0.07);
    assert!((asb - lhp).abs() < 2e-3);
}
