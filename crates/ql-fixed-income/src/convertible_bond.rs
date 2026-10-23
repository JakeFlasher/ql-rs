use ql_ad::{Dual, Real};
#[derive(Clone, Debug)]
pub struct ConvertibleBondSpec<T: Real> {
    pub spot: T,
    pub conversion_ratio: T,
    pub face_value: T,
    pub coupon_rate: T,
    pub coupon_times: Vec<T>,
    pub maturity: T,
    pub rate: T,
    pub div_yield: T,
    pub credit_spread: T,
    pub vol: T,
    pub call_price: Option<T>,
    pub put_price: Option<T>,
}
pub fn price_convertible_bond_tf<T: Real>(s: &ConvertibleBondSpec<T>, steps: usize) -> (T, T) {
    let n = steps.max(1);
    let dt = s.maturity / n as f64;
    let u = (s.vol * dt.sqrt()).exp();
    let d = T::one() / u;
    let p = (((s.rate - s.div_yield) * dt).exp() - d) / (u - d);
    let mut v = vec![T::zero(); n + 1];
    let mut b = v.clone();
    let coupon_at = |time: T| -> T {
        let mut previous = T::zero();
        let mut amount = T::zero();
        for &ct in &s.coupon_times {
            if (ct - time).val().abs() <= dt.val().abs() * 0.5 + 1e-12 {
                amount = s.face_value * s.coupon_rate * (ct - previous);
                break;
            }
            previous = ct;
        }
        amount
    };
    let maturity_coupon = coupon_at(s.maturity);
    for i in 0..=n {
        let spot = s.spot * u.powf((n - i) as f64) * d.powf(i as f64);
        let conv = s.conversion_ratio * spot;
        let cash = s.face_value + maturity_coupon;
        v[i] = conv.max_val(cash);
        b[i] = if conv.val() >= cash.val() {
            T::zero()
        } else {
            cash
        };
    }
    for k in (0..n).rev() {
        for i in 0..=k {
            let ev = p * v[i] + (T::one() - p) * v[i + 1];
            let eb = p * b[i] + (T::one() - p) * b[i + 1];
            let mut val =
                (-s.rate * dt).exp() * (ev - eb) + (-(s.rate + s.credit_spread) * dt).exp() * eb;
            let mut bb = (-(s.rate + s.credit_spread) * dt).exp() * eb;
            let coupon = coupon_at(dt * k as f64);
            val += coupon;
            bb += coupon;
            if let Some(c) = s.call_price
                && val.val() > c.val()
            {
                val = c;
                bb = c;
            }
            if let Some(q) = s.put_price
                && val.val() < q.val()
            {
                val = q;
                bb = q;
            }
            let spot = s.spot * u.powf((k - i) as f64) * d.powf(i as f64);
            let conv = s.conversion_ratio * spot;
            if conv.val() > val.val() {
                val = conv;
                bb = T::zero();
            }
            v[i] = val;
            b[i] = bb;
        }
    }
    (v[0], b[0])
}
pub fn convertible_bond_greeks_dual4(
    s: &ConvertibleBondSpec<f64>,
    steps: usize,
) -> (f64, [f64; 4]) {
    let d: ConvertibleBondSpec<Dual<4>> = ConvertibleBondSpec {
        spot: Dual::<4>::variable(s.spot, 0),
        conversion_ratio: Dual::from_value(s.conversion_ratio),
        face_value: Dual::from_value(s.face_value),
        coupon_rate: Dual::from_value(s.coupon_rate),
        coupon_times: s
            .coupon_times
            .iter()
            .map(|&x| Dual::from_value(x))
            .collect(),
        maturity: Dual::from_value(s.maturity),
        rate: Dual::variable(s.rate, 2),
        div_yield: Dual::from_value(s.div_yield),
        credit_spread: Dual::variable(s.credit_spread, 3),
        vol: Dual::variable(s.vol, 1),
        call_price: s.call_price.map(Dual::from_value),
        put_price: s.put_price.map(Dual::from_value),
    };
    let (v, _) = price_convertible_bond_tf(&d, steps);
    (v.val, [v.d[0], v.d[1], v.d[2], v.d[3]])
}
