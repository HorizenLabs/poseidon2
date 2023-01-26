use lazy_static::lazy_static;
use std::sync::Arc;

use crate::{
    fields::{bls12::FpBLS12, bn256::FpBN256, pallas::FpPallas, vesta::FpVesta, goldilocks::FpGoldiLocks, babybear::FpBabyBear},
    neptune::neptune_params::NeptuneParams,
};

lazy_static! {
    // Number of partial rounds:
    // ceil(1.125 * ceil(log_d(2) * (min(kappa, log_2(p)) - 6) + 3 + t + log_d(t)))
    // BN256
    pub static ref NEPTUNE_BN_PARAMS: Arc<NeptuneParams<FpBN256>> = Arc::new(NeptuneParams::new(4, 5, 6, 68));
    // BLS12
    pub static ref NEPTUNE_BLS_PARAMS: Arc<NeptuneParams<FpBLS12>> = Arc::new(NeptuneParams::new(4, 5, 6, 68));
    // Goldilocks
    pub static ref NEPTUNE_GOLDILOCKS_PARAMS: Arc<NeptuneParams<FpGoldiLocks>> = Arc::new(NeptuneParams::new(12, 7, 6, 42));
    // BabyBear
    pub static ref NEPTUNE_BABYBEAR_PARAMS: Arc<NeptuneParams<FpBabyBear>> = Arc::new(NeptuneParams::new(24, 7, 6, 43));
    // Pallas
    pub static ref NEPTUNE_PALLAS_PARAMS: Arc<NeptuneParams<FpPallas>> = Arc::new(NeptuneParams::new(4, 5, 6, 68));
    // Vesta
    pub static ref NEPTUNE_VESTA_PARAMS: Arc<NeptuneParams<FpVesta>> = Arc::new(NeptuneParams::new(4, 5, 6, 68));
}