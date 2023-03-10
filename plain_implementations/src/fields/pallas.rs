use ark_ff::fields::{Fp256, MontBackend, MontConfig};
use std::convert::TryInto;

#[derive(MontConfig)]
#[modulus = "28948022309329048855892746252171976963363056481941560715954676764349967630337"]
#[generator = "5"]
pub struct FqConfig;
pub type FpPallas = Fp256<MontBackend<FqConfig, 4>>;