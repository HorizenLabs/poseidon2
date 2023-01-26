use super::gmimc_params::GmimcParams;
use crate::fields::goldilocks::FpGoldiLocks;

use lazy_static::lazy_static;
use std::sync::Arc;

type Scalar = FpGoldiLocks;

lazy_static! {
    // Number of rounds:
    // max(2 + 2 * (t + t^2), ceil(2 * log_d(p)) + 2 * t)
    pub static ref GMIMC_GOLDILOCKS_PARAMS: Arc<GmimcParams<Scalar>> =
        Arc::new(GmimcParams::new(12, 7, 314));
}
