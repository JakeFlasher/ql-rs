use crate::curves::DiscountCurve;
use ql_ad::Real;
use ql_ad::TapeVar;
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
impl<T: Real + 'static> MarketContext<T> {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn with_reference_date(mut self, reference_date: Date) -> Self {
        self.reference_date = reference_date;
        self
    }
    pub fn with_discount_curve(mut self, id: DiscountCurveId, c: DiscountCurve<T>) -> Self {
        self.insert_discount_curve(id, c);
        self
    }
    pub fn with_forward_curve(mut self, id: ForwardIndexId, c: DiscountCurve<T>) -> Self {
        self.insert_forward_curve(id, c);
        self
    }
    pub fn insert_discount_curve(&mut self, id: DiscountCurveId, c: DiscountCurve<T>) {
        let name = format!("{id:?}");
        self.named_discount_curves.insert(name, Arc::new(c.clone()));
        self.discount_curves.insert(id, c);
    }
    pub fn insert_forward_curve(&mut self, id: ForwardIndexId, c: DiscountCurve<T>) {
        let name = format!("{id:?}");
        self.named_forward_curves.insert(name, Arc::new(c.clone()));
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
        if let Some(c) = self.named_discount_curves.get(name) {
            return Some(c.as_ref());
        }
        let id = match name.to_ascii_uppercase().as_str() {
            "USD_SOFR" | "USD-SOFR" => DiscountCurveId::UsdSofr,
            "EUR_ESTR" | "EUR-ESTR" => DiscountCurveId::EurEstr,
            "GBP_SONIA" | "GBP-SONIA" => DiscountCurveId::GbpSonia,
            _ => return None,
        };
        self.discount_curves
            .get(&id)
            .map(|c| c as &dyn crate::curves::YieldTermStructure<T>)
    }
    pub fn forward(&self, name: &str) -> Option<&dyn crate::curves::YieldTermStructure<T>> {
        if let Some(c) = self.named_forward_curves.get(name) {
            return Some(c.as_ref());
        }
        let id = match name.to_ascii_uppercase().as_str() {
            "USD_SOFR" | "USD-SOFR" => ForwardIndexId::UsdSofr,
            "USD_IBOR_3M" | "USD-IBOR-3M" => ForwardIndexId::UsdIbor3m,
            "EUR_EURIBOR_3M" | "EUR-EURIBOR-3M" => ForwardIndexId::EurEuribor3m,
            "EUR_EURIBOR_6M" | "EUR-EURIBOR-6M" => ForwardIndexId::EurEuribor6m,
            "GBP_SONIA" | "GBP-SONIA" => ForwardIndexId::GbpSonia,
            _ => return None,
        };
        self.forward_curves
            .get(&id)
            .map(|c| c as &dyn crate::curves::YieldTermStructure<T>)
    }
    pub fn spot(&self, name: &str) -> Option<T> {
        self.spots.get(name).copied()
    }
    pub fn vol(&self, name: &str) -> Option<T> {
        self.vols.get(name).copied()
    }
}

impl MarketContext<f64> {
    /// Lift all typed curve nodes and scalar market observables onto a fresh
    /// reverse-mode scalar type.  The structure (pillar times and identifiers)
    /// is retained, so callers can record one coherent tape pass over the
    /// context rather than rebuilding ad-hoc maps.
    pub fn to_ad(&self) -> MarketContext<TapeVar<f64>> {
        let mut out =
            MarketContext::<TapeVar<f64>>::default().with_reference_date(self.reference_date);
        for (id, curve) in &self.discount_curves {
            let vars = curve
                .discounts
                .iter()
                .copied()
                .map(TapeVar::from_f64)
                .collect::<Vec<_>>();
            out.insert_discount_curve(*id, curve.to_ad_curve(&vars));
        }
        for (id, curve) in &self.forward_curves {
            let vars = curve
                .discounts
                .iter()
                .copied()
                .map(TapeVar::from_f64)
                .collect::<Vec<_>>();
            out.insert_forward_curve(*id, curve.to_ad_curve(&vars));
        }
        out.spots = self
            .spots
            .iter()
            .map(|(k, v)| (k.clone(), TapeVar::from_f64(*v)))
            .collect();
        out.vols = self
            .vols
            .iter()
            .map(|(k, v)| (k.clone(), TapeVar::from_f64(*v)))
            .collect();
        out
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
