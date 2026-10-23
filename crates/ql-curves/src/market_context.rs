use crate::curves::DiscountCurve;
use ql_ad::Real;
use std::collections::HashMap;
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DiscountCurveId {
    UsdSofr,
    EurEstr,
    GbpSonia,
    Custom(u16),
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ForwardIndexId {
    UsdSofr,
    UsdIbor3m,
    EurEuribor3m,
    EurEuribor6m,
    GbpSonia,
    Custom(u16),
}
#[derive(Debug, Clone)]
pub struct MarketContext<T: Real = f64> {
    pub discount_curves: HashMap<DiscountCurveId, DiscountCurve<T>>,
    pub forward_curves: HashMap<ForwardIndexId, DiscountCurve<T>>,
}
impl<T: Real> Default for MarketContext<T> {
    fn default() -> Self {
        Self {
            discount_curves: HashMap::new(),
            forward_curves: HashMap::new(),
        }
    }
}
impl<T: Real> MarketContext<T> {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn with_discount_curve(mut self, id: DiscountCurveId, c: DiscountCurve<T>) -> Self {
        self.discount_curves.insert(id, c);
        self
    }
    pub fn with_forward_curve(mut self, id: ForwardIndexId, c: DiscountCurve<T>) -> Self {
        self.forward_curves.insert(id, c);
        self
    }
    pub fn insert_discount_curve(&mut self, id: DiscountCurveId, c: DiscountCurve<T>) {
        self.discount_curves.insert(id, c);
    }
    pub fn insert_forward_curve(&mut self, id: ForwardIndexId, c: DiscountCurve<T>) {
        self.forward_curves.insert(id, c);
    }
    pub fn discount_curve(&self, id: DiscountCurveId) -> Option<&DiscountCurve<T>> {
        self.discount_curves.get(&id)
    }
    pub fn forward_curve(&self, id: ForwardIndexId) -> Option<&DiscountCurve<T>> {
        self.forward_curves.get(&id)
    }
}

impl DiscountCurveId {
    pub const USD_SOFR: Self = Self::UsdSofr;
    pub const EUR_ESTR: Self = Self::EurEstr;
    pub const GBP_SONIA: Self = Self::GbpSonia;
}
impl ForwardIndexId {
    pub const USD_SOFR: Self = Self::UsdSofr;
    pub const USD_IBOR_3M: Self = Self::UsdIbor3m;
    pub const EUR_EURIBOR_3M: Self = Self::EurEuribor3m;
    pub const EUR_EURIBOR_6M: Self = Self::EurEuribor6m;
    pub const GBP_SONIA: Self = Self::GbpSonia;
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CurveId {
    Discount(DiscountCurveId),
    Forward(ForwardIndexId),
}
