#[derive(Clone, Debug)]
pub struct SinhGrid1D {
    pub nodes: Vec<f64>,
    pub first_minus: Vec<f64>,
    pub first_zero: Vec<f64>,
    pub first_plus: Vec<f64>,
    pub second_minus: Vec<f64>,
    pub second_zero: Vec<f64>,
    pub second_plus: Vec<f64>,
}
impl SinhGrid1D {
    pub fn new(n: usize, min: f64, max: f64, focus: f64, c: f64) -> Self {
        let n = n.max(3);
        let lo = ((min - focus) / c).asinh();
        let hi = ((max - focus) / c).asinh();
        let nodes = (0..n)
            .map(|i| focus + c * (lo + (hi - lo) * i as f64 / (n - 1) as f64).sinh())
            .collect::<Vec<_>>();
        let mut x = Self {
            nodes,
            first_minus: vec![0.; n],
            first_zero: vec![0.; n],
            first_plus: vec![0.; n],
            second_minus: vec![0.; n],
            second_zero: vec![0.; n],
            second_plus: vec![0.; n],
        };
        for i in 1..n - 1 {
            let hm = x.nodes[i] - x.nodes[i - 1];
            let hp = x.nodes[i + 1] - x.nodes[i];
            x.first_minus[i] = -hp / (hm * (hm + hp));
            x.first_zero[i] = (hp - hm) / (hm * hp);
            x.first_plus[i] = hm / (hp * (hm + hp));
            x.second_minus[i] = 2. / (hm * (hm + hp));
            x.second_zero[i] = -2. / (hm * hp);
            x.second_plus[i] = 2. / (hp * (hm + hp));
        }
        x
    }
}
#[derive(Clone, Debug)]
pub struct HestonGrid2D {
    pub x: SinhGrid1D,
    pub v: SinhGrid1D,
}
impl HestonGrid2D {
    pub fn new(nx: usize, nv: usize, spot: f64, v0: f64) -> Self {
        Self::new_centered(nx, nv, spot.ln(), v0)
    }
    pub fn new_centered(nx: usize, nv: usize, center_x: f64, v0: f64) -> Self {
        let v_focus = if v0 > 0.1 { v0 } else { 0.04 };
        Self {
            x: SinhGrid1D::new(nx, center_x - 3.0, center_x + 3.0, center_x, 0.4),
            v: SinhGrid1D::new(nv, 0.0, (v_focus * 8.0).max(0.6), v_focus, 0.05),
        }
    }
}
