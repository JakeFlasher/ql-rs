use ql_ad::Real;
use ql_vol::{bachelier_price, black76_price};
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OptionType {
    Call,
    Put,
}
impl OptionType {
    pub fn is_call(self) -> bool {
        matches!(self, Self::Call)
    }
}
#[derive(Debug, Clone, Copy)]
pub struct EuropeanOption {
    pub forward: f64,
    pub strike: f64,
    pub expiry: f64,
    pub volatility: f64,
    pub discount: f64,
    pub option_type: OptionType,
}
impl EuropeanOption {
    pub fn price<T: Real>(&self) -> T {
        black76_price(
            T::from_f64(self.forward),
            T::from_f64(self.strike),
            T::from_f64(self.volatility),
            T::from_f64(self.expiry),
            T::from_f64(self.discount),
            self.option_type.is_call(),
        )
    }
    pub fn black76_price<T: Real>(&self) -> T {
        self.price()
    }
    pub fn bachelier_price<T: Real>(&self) -> T {
        bachelier_price(
            T::from_f64(self.forward),
            T::from_f64(self.strike),
            T::from_f64(self.volatility),
            T::from_f64(self.expiry),
            T::from_f64(self.discount),
            self.option_type.is_call(),
        )
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BarrierType {
    UpAndOut,
    UpAndIn,
    DownAndOut,
    DownAndIn,
}
// Historical ql-exotics names are retained as associated aliases while the
// workspace uses this single canonical barrier type.
impl BarrierType {
    #[allow(non_upper_case_globals)]
    pub const UpOut: Self = Self::UpAndOut;
    #[allow(non_upper_case_globals)]
    pub const UpIn: Self = Self::UpAndIn;
    #[allow(non_upper_case_globals)]
    pub const DownOut: Self = Self::DownAndOut;
    #[allow(non_upper_case_globals)]
    pub const DownIn: Self = Self::DownAndIn;
}
#[derive(Debug, Clone, Copy)]
pub struct BarrierOption {
    pub underlying: f64,
    pub strike: f64,
    pub barrier: f64,
    pub expiry: f64,
    pub volatility: f64,
    pub rate: f64,
    pub option_type: OptionType,
    pub barrier_type: BarrierType,
    pub rebate: f64,
}
impl BarrierOption {
    pub fn price<T: Real>(&self) -> T {
        crate::pde1d::Fd1dBlackScholesEngine::default().price_barrier_option(self)
    }
}
#[derive(Debug, Clone, Copy)]
pub struct AmericanOption {
    pub spot: f64,
    pub strike: f64,
    pub expiry: f64,
    pub volatility: f64,
    pub rate: f64,
    pub option_type: OptionType,
}
impl AmericanOption {
    pub fn price<T: Real>(&self) -> T {
        crate::pde1d::Fd1dBlackScholesEngine::default().price_american(self)
    }
    pub fn new(
        spot: f64,
        strike: f64,
        expiry: f64,
        volatility: f64,
        rate: f64,
        option_type: OptionType,
    ) -> Self {
        Self {
            spot,
            strike,
            expiry,
            volatility,
            rate,
            option_type,
        }
    }
}
impl EuropeanOption {
    pub fn new(
        forward: f64,
        strike: f64,
        expiry: f64,
        volatility: f64,
        discount: f64,
        option_type: OptionType,
    ) -> Self {
        Self {
            forward,
            strike,
            expiry,
            volatility,
            discount,
            option_type,
        }
    }
    pub fn greeks<T: Real>(&self) -> (T, T, T, T, T, T, T, T) {
        ql_vol::black76_greeks(
            T::from_f64(self.forward),
            T::from_f64(self.strike),
            T::from_f64(self.volatility),
            T::from_f64(self.expiry),
            T::from_f64(self.discount),
            self.option_type.is_call(),
        )
    }
}
impl BarrierOption {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        underlying: f64,
        strike: f64,
        barrier: f64,
        expiry: f64,
        volatility: f64,
        rate: f64,
        option_type: OptionType,
        barrier_type: BarrierType,
    ) -> Self {
        Self {
            underlying,
            strike,
            barrier,
            expiry,
            volatility,
            rate,
            option_type,
            barrier_type,
            rebate: 0.0,
        }
    }
}
