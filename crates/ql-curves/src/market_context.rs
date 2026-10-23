use crate::curves::DiscountCurve;
use ql_ad::Real;
use ql_time::Date;
use std::collections::HashMap;
use std::sync::Arc;
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
/// Immutable market snapshot used by higher-level pricing workflows.  The
/// typed identifiers are retained for backwards compatibility; string keys
/// provide an escape hatch for user-defined indices and are the canonical
/// representation at application boundaries.
#[derive(Clone)]
pub struct MarketContext<T: Real = f64> {
    pub reference_date: Date,
    pub discount_curves: HashMap<DiscountCurveId, DiscountCurve<T>>,
    pub forward_curves: HashMap<ForwardIndexId, DiscountCurve<T>>,
    pub named_discount_curves: HashMap<String, Arc<dyn crate::curves::YieldTermStructure<T>>>,
    pub named_forward_curves: HashMap<String, Arc<dyn crate::curves::YieldTermStructure<T>>>,
    pub spots: HashMap<String, T>,
    pub vols: HashMap<String, T>,
}
impl<T: Real> std::fmt::Debug for MarketContext<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MarketContext")
            .field("reference_date", &self.reference_date)
            .field(
                "discount_curves",
                &self.discount_curves.keys().collect::<Vec<_>>(),
            )
            .field(
                "forward_curves",
                &self.forward_curves.keys().collect::<Vec<_>>(),
            )
            .field(
                "named_discount_curves",
                &self.named_discount_curves.keys().collect::<Vec<_>>(),
            )
            .field(
                "named_forward_curves",
                &self.named_forward_curves.keys().collect::<Vec<_>>(),
            )
            .field("spots", &self.spots)
            .field("vols", &self.vols)
            .finish()
    }
}
impl<T: Real> Default for MarketContext<T> {
    fn default() -> Self {
        Self {
            reference_date: Date(0),
            discount_curves: HashMap::new(),
            forward_curves: HashMap::new(),
            named_discount_curves: HashMap::new(),
            named_forward_curves: HashMap::new(),
            spots: HashMap::new(),
            vols: HashMap::new(),
        }
    }
}
impl<T: Real> MarketContext<T> {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn with_reference_date(mut self, reference_date: Date) -> Self {
        self.reference_date = reference_date;
        self
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
    pub fn with_discount(
        mut self,
        name: impl Into<String>,
        curve: Arc<dyn crate::curves::YieldTermStructure<T>>,
    ) -> Self {
        self.named_discount_curves.insert(name.into(), curve);
        self
    }
    pub fn with_forward(
        mut self,
        name: impl Into<String>,
        curve: Arc<dyn crate::curves::YieldTermStructure<T>>,
    ) -> Self {
        self.named_forward_curves.insert(name.into(), curve);
        self
    }
    pub fn with_spot(mut self, name: impl Into<String>, value: T) -> Self {
        self.spots.insert(name.into(), value);
        self
    }
    pub fn with_vol(mut self, name: impl Into<String>, value: T) -> Self {
        self.vols.insert(name.into(), value);
        self
    }
    pub fn discount(&self, name: &str) -> Option<&dyn crate::curves::YieldTermStructure<T>> {
        self.named_discount_curves.get(name).map(Arc::as_ref)
    }
    pub fn forward(&self, name: &str) -> Option<&dyn crate::curves::YieldTermStructure<T>> {
        self.named_forward_curves.get(name).map(Arc::as_ref)
    }
    pub fn spot(&self, name: &str) -> Option<T> {
        self.spots.get(name).copied()
    }
    pub fn vol(&self, name: &str) -> Option<T> {
        self.vols.get(name).copied()
    }
}

pub type MarketEnvironment<T = f64> = MarketContext<T>;

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
