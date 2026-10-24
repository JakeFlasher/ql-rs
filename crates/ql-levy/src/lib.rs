use ql_ad::Real;
pub trait LevyModel<T: Real> {
    fn martingale_compensator(&self) -> T;
    fn char_exponent(&self, u: T) -> (T, T);
    fn cumulants(&self, rate: T, div_yield: T, expiry: T) -> (T, T, T);
}
pub mod cgmy;
pub mod cos;
pub mod fft;
pub mod kou;
pub mod merton;
pub mod vg;
pub use cgmy::*;
pub use cos::*;
pub use fft::*;
pub use kou::*;
pub use merton::*;
pub use vg::*;
