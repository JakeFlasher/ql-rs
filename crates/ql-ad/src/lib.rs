#![feature(portable_simd)]
#![cfg_attr(feature = "enzyme", feature(autodiff))]
pub mod dual;
pub mod ift;
pub mod scalar;
pub mod tape;
pub use dual::*;
pub use ift::*;
pub use scalar::*;
pub use tape::*;

mod phase30;
#[allow(unused_imports)]
pub use phase30::*;
mod phase38;
#[allow(unused_imports)]
pub use phase38::*;
