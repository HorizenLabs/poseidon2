use super::poseidon_params::PoseidonParams;
use crate::merkle_tree::merkle_tree_fp::MerkleTreeHash;
use ark_ff::PrimeField;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct Poseidon<S: PrimeField> {
    pub(crate) params: Arc<PoseidonParams<S>>,
}

impl<S: PrimeField> Poseidon<S> {
    pub fn new(params: &Arc<PoseidonParams<S>>) -> Self {
        Poseidon {
            params: Arc::clone(params),
        }
    }

    pub fn get_t(&self) -> usize {
        self.params.t
    }

    pub fn permutation(&self, input: &[S]) -> Vec<S> {
        let t = self.params.t;
        assert_eq!(input.len(), t);

        let mut current_state = input.to_owned();
        for r in 0..self.params.rounds_f_beginning {
            current_state = self.add_rc(&current_state, &self.params.round_constants[r]);
            current_state = self.sbox(&current_state);
            current_state = self.matmul(&current_state, &self.params.mds);
        }
        let p_end = self.params.rounds_f_beginning + self.params.rounds_p;
        current_state = self.add_rc(&current_state, &self.params.opt_round_constants[0]);
        current_state = self.matmul(&current_state, &self.params.m_i);

        for r in self.params.rounds_f_beginning..p_end {
            current_state[0] = self.sbox_p(&current_state[0]);
            if r < p_end - 1 {
                current_state[0].add_assign(
                    &self.params.opt_round_constants[r + 1 - self.params.rounds_f_beginning][0],
                );
            }
            current_state = self.cheap_matmul(&current_state, p_end - r - 1);
        }
        for r in p_end..self.params.rounds {
            current_state = self.add_rc(&current_state, &self.params.round_constants[r]);
            current_state = self.sbox(&current_state);
            current_state = self.matmul(&current_state, &self.params.mds);
        }
        current_state
    }

    pub fn permutation_not_opt(&self, input: &[S]) -> Vec<S> {
        let t = self.params.t;
        assert_eq!(input.len(), t);

        let mut current_state = input.to_owned();

        for r in 0..self.params.rounds_f_beginning {
            current_state = self.add_rc(&current_state, &self.params.round_constants[r]);
            current_state = self.sbox(&current_state);
            current_state = self.matmul(&current_state, &self.params.mds);
        }
        let p_end = self.params.rounds_f_beginning + self.params.rounds_p;
        for r in self.params.rounds_f_beginning..p_end {
            current_state = self.add_rc(&current_state, &self.params.round_constants[r]);
            current_state[0] = self.sbox_p(&current_state[0]);
            current_state = self.matmul(&current_state, &self.params.mds);
        }
        for r in p_end..self.params.rounds {
            current_state = self.add_rc(&current_state, &self.params.round_constants[r]);
            current_state = self.sbox(&current_state);
            current_state = self.matmul(&current_state, &self.params.mds);
        }
        current_state
    }

    fn sbox(&self, input: &[S]) -> Vec<S> {
        input.iter().map(|el| self.sbox_p(el)).collect()
    }

    fn sbox_p(&self, input: &S) -> S {
        let mut input2 = *input;
        input2.square_in_place();

        match self.params.d {
            3 => {
                let mut out = input2;
                out.mul_assign(input);
                out
            }
            5 => {
                let mut out = input2;
                out.square_in_place();
                out.mul_assign(input);
                out
            }
            7 => {
                let mut out = input2;
                out.square_in_place();
                out.mul_assign(&input2);
                out.mul_assign(input);
                out
            }
            _ => {
                panic!()
            }
        }
    }

    fn cheap_matmul(&self, input: &[S], r: usize) -> Vec<S> {
        let v = &self.params.v[r];
        let w_hat = &self.params.w_hat[r];
        let t = self.params.t;

        let mut new_state = vec![S::zero(); t];
        new_state[0] = self.params.mds[0][0];
        new_state[0].mul_assign(&input[0]);
        for i in 1..t {
            let mut tmp = w_hat[i - 1];
            tmp.mul_assign(&input[i]);
            new_state[0].add_assign(&tmp);
        }
        for i in 1..t {
            new_state[i] = input[0];
            new_state[i].mul_assign(&v[i - 1]);
            new_state[i].add_assign(&input[i]);
        }

        new_state
    }

    fn matmul(&self, input: &[S], mat: &[Vec<S>]) -> Vec<S> {
        let t = mat.len();
        debug_assert!(t == input.len());
        let mut out = vec![S::zero(); t];
        for row in 0..t {
            for (col, inp) in input.iter().enumerate().take(t) {
                let mut tmp = mat[row][col];
                tmp.mul_assign(inp);
                out[row].add_assign(&tmp);
            }
        }
        out
    }

    fn add_rc(&self, input: &[S], rc: &[S]) -> Vec<S> {
        input
            .iter()
            .zip(rc.iter())
            .map(|(a, b)| {
                let mut r = *a;
                r.add_assign(b);
                r
            })
            .collect()
    }
}

impl<F: PrimeField> MerkleTreeHash<F> for Poseidon<F> {
    fn compress(&self, input: &[&F; 2]) -> F {
        self.permutation(&[input[0].to_owned(), input[1].to_owned(), F::zero()])[0]
    }
}

#[cfg(test)]
mod poseidon_tests_bls12 {
    use super::*;
    use crate::fields::{bls12::FpBLS12, utils::from_hex, utils::random_scalar};
    use crate::poseidon::poseidon_instance_bls12::{
        POSEIDON_BLS_3_PARAMS,
        POSEIDON_BLS_4_PARAMS,
        POSEIDON_BLS_8_PARAMS,
    };

    type Scalar = FpBLS12;

    static TESTRUNS: usize = 5;

    #[test]
    fn consistent_perm() {
        let instances = vec![
            Poseidon::new(&POSEIDON_BLS_3_PARAMS),
            Poseidon::new(&POSEIDON_BLS_4_PARAMS),
            Poseidon::new(&POSEIDON_BLS_8_PARAMS)
        ];
        for instance in instances {
            let t = instance.params.t;
            for _ in 0..TESTRUNS {
                let input1: Vec<Scalar> = (0..t).map(|_| random_scalar()).collect();

                let mut input2: Vec<Scalar>;
                loop {
                    input2 = (0..t).map(|_| random_scalar()).collect();
                    if input1 != input2 {
                        break;
                    }
                }

                let perm1 = instance.permutation(&input1);
                let perm2 = instance.permutation(&input1);
                let perm3 = instance.permutation(&input2);
                assert_eq!(perm1, perm2);
                assert_ne!(perm1, perm3);
            }
        }
    }

    #[test]
    fn kats() {
        let poseidon = Poseidon::new(&POSEIDON_BLS_3_PARAMS);
        let input: Vec<Scalar> = vec![Scalar::from(0), Scalar::from(1), Scalar::from(2)];
        let perm = poseidon.permutation(&input);
        assert_eq!(
            perm[0],
            from_hex("0x22a1f1595d99e4a04fc0a5b16be51a844a7cb5b5d69627ebbd1ee8142e7532ce")
        );
        assert_eq!(
            perm[1],
            from_hex("0x5f7ae6d6c380c90510de9c045ee75163eae24054ba8cd88d254cd1c343f43176")
        );
        assert_eq!(
            perm[2],
            from_hex("0x1b7e4da7d1ac6accb2e0470a83ba87d7bb585f4ba8c9a34f936faf3b3dfc695b")
        );
    }
    #[test]
    fn opt_equals_not_opt() {
        let instances = vec![
            Poseidon::new(&POSEIDON_BLS_3_PARAMS),
            Poseidon::new(&POSEIDON_BLS_4_PARAMS),
            Poseidon::new(&POSEIDON_BLS_8_PARAMS)
        ];
        for instance in instances {
            let t = instance.params.t;
            for _ in 0..TESTRUNS {
                let input: Vec<Scalar> = (0..t).map(|_| random_scalar()).collect();

                let perm1 = instance.permutation(&input);
                let perm2 = instance.permutation_not_opt(&input);
                assert_eq!(perm1, perm2);
            }
        }
    }
}

#[cfg(test)]
mod poseidon_tests_bn256 {
    use super::*;
    use crate::fields::{bn256::FpBN256, utils::from_hex, utils::random_scalar};
    use crate::poseidon::poseidon_instance_bn256::POSEIDON_BN_PARAMS;

    type Scalar = FpBN256;

    static TESTRUNS: usize = 5;

    #[test]
    fn consistent_perm() {
        let poseidon = Poseidon::new(&POSEIDON_BN_PARAMS);
        let t = poseidon.params.t;
        for _ in 0..TESTRUNS {
            let input1: Vec<Scalar> = (0..t).map(|_| random_scalar()).collect();

            let mut input2: Vec<Scalar>;
            loop {
                input2 = (0..t).map(|_| random_scalar()).collect();
                if input1 != input2 {
                    break;
                }
            }

            let perm1 = poseidon.permutation(&input1);
            let perm2 = poseidon.permutation(&input1);
            let perm3 = poseidon.permutation(&input2);
            assert_eq!(perm1, perm2);
            assert_ne!(perm1, perm3);
        }
    }

    #[test]
    fn kats() {
        let poseidon = Poseidon::new(&POSEIDON_BN_PARAMS);
        let input: Vec<Scalar> = vec![Scalar::from(0), Scalar::from(1), Scalar::from(2)];
        let perm = poseidon.permutation(&input);
        assert_eq!(
            perm[0],
            from_hex("0x2e72c60509a284872f62830b58ed8524a58c362dd3ddb98b2767f36b566596bd")
        );
        assert_eq!(
            perm[1],
            from_hex("0x180a812301272545f79ae1012b0425162a1833ac39101e070732f4d8a8bc4718")
        );
        assert_eq!(
            perm[2],
            from_hex("0x1828343d70eed99aae404e3ea58209f45743f3d54983fe250ce1526a9d8cf88e")
        );
    }
    #[test]
    fn opt_equals_not_opt() {
        let poseidon = Poseidon::new(&POSEIDON_BN_PARAMS);
        let t = poseidon.params.t;
        for _ in 0..TESTRUNS {
            let input: Vec<Scalar> = (0..t).map(|_| random_scalar()).collect();

            let perm1 = poseidon.permutation(&input);
            let perm2 = poseidon.permutation_not_opt(&input);
            assert_eq!(perm1, perm2);
        }
    }
}

#[cfg(test)]
#[allow(unused_imports)]
mod poseidon_tests_goldilocks {
    use super::*;
    use crate::fields::{goldilocks::FpGoldiLocks, utils::from_hex, utils::random_scalar};
    use crate::poseidon::poseidon_instance_goldilocks::{
        POSEIDON_GOLDILOCKS_8_PARAMS,
        POSEIDON_GOLDILOCKS_12_PARAMS,
        POSEIDON_GOLDILOCKS_16_PARAMS,
        POSEIDON_GOLDILOCKS_20_PARAMS,
    };
    use std::convert::TryFrom;

    type Scalar = FpGoldiLocks;
    use ark_ff::UniformRand;

    static TESTRUNS: usize = 5;

    #[test]
    fn consistent_perm() {
        let instances = vec![
            Poseidon::new(&POSEIDON_GOLDILOCKS_8_PARAMS),
            Poseidon::new(&POSEIDON_GOLDILOCKS_12_PARAMS),
            Poseidon::new(&POSEIDON_GOLDILOCKS_16_PARAMS),
            Poseidon::new(&POSEIDON_GOLDILOCKS_20_PARAMS),
        ];
        for instance in instances {
            let t = instance.params.t;
            for _ in 0..TESTRUNS {
                let input1: Vec<Scalar> = (0..t).map(|_| random_scalar()).collect();

                let mut input2: Vec<Scalar>;
                loop {
                    input2 = (0..t).map(|_| random_scalar()).collect();
                    if input1 != input2 {
                        break;
                    }
                }

                let perm1 = instance.permutation(&input1);
                let perm2 = instance.permutation(&input1);
                let perm3 = instance.permutation(&input2);
                assert_eq!(perm1, perm2);
                assert_ne!(perm1, perm3);
            }
        }
    }

    #[test]
    fn kats() {
        let poseidon = Poseidon::new(&POSEIDON_GOLDILOCKS_12_PARAMS);
        // let input: Vec<Scalar> = vec![Scalar::zero(), Scalar::one(), utils::from_u64::<Scalar>(2)];
        let mut input: Vec<Scalar> = vec![];
        for i in 0..poseidon.params.t {
            input.push(Scalar::from(i as u64));
        }
        let perm = poseidon.permutation(&input);
        assert_eq!(
            perm[0],
            from_hex("0xe9ad770762f48ef5")
        );
        assert_eq!(
            perm[1],
            from_hex("0xc12796961ddc7859")
        );
        assert_eq!(
            perm[2],
            from_hex("0xa61b71de9595e016")
        );
        assert_eq!(
            perm[3],
            from_hex("0xead9e6aa583aafa3")
        );
        assert_eq!(
            perm[4],
            from_hex("0x93e297beff76e95b")
        );
        assert_eq!(
            perm[5],
            from_hex("0x53abd3c5c2a0e924")
        );
        assert_eq!(
            perm[6],
            from_hex("0xf3bc50e655c74f51")
        );
        assert_eq!(
            perm[7],
            from_hex("0x246cac41b9a45d84")
        );
        assert_eq!(
            perm[8],
            from_hex("0xcc7f9314b2341f4f")
        );
        assert_eq!(
            perm[9],
            from_hex("0xf5f071587c83415c")
        );
        assert_eq!(
            perm[10],
            from_hex("0x09486cf35116fba3")
        );
        assert_eq!(
            perm[11],
            from_hex("0x9d82aaf136b5c38a")
        );
    }

    #[test]
    fn opt_equals_not_opt() {
        let instances = vec![
            Poseidon::new(&POSEIDON_GOLDILOCKS_8_PARAMS),
            Poseidon::new(&POSEIDON_GOLDILOCKS_12_PARAMS),
            Poseidon::new(&POSEIDON_GOLDILOCKS_16_PARAMS),
            Poseidon::new(&POSEIDON_GOLDILOCKS_20_PARAMS),
        ];
        for instance in instances {
            let t = instance.params.t;
            for _ in 0..TESTRUNS {
                let input: Vec<Scalar> = (0..t).map(|_| random_scalar()).collect();

                let perm1 = instance.permutation(&input);
                let perm2 = instance.permutation_not_opt(&input);
                assert_eq!(perm1, perm2);
            }
        }
    }
}

#[cfg(test)]
#[allow(unused_imports)]
mod poseidon_tests_babybear {
    use super::*;
    use crate::fields::{babybear::FpBabyBear, utils::from_hex, utils::random_scalar};
    use crate::poseidon::poseidon_instance_babybear::POSEIDON_BABYBEAR_PARAMS;

    type Scalar = FpBabyBear;
    use ark_ff::UniformRand;

    static TESTRUNS: usize = 5;

    #[test]
    fn consistent_perm() {
        let poseidon = Poseidon::new(&POSEIDON_BABYBEAR_PARAMS);
        let t = poseidon.params.t;
        for _ in 0..TESTRUNS {
            let input1: Vec<Scalar> = (0..t).map(|_| random_scalar()).collect();

            let mut input2: Vec<Scalar>;
            loop {
                input2 = (0..t).map(|_| random_scalar()).collect();
                if input1 != input2 {
                    break;
                }
            }

            let perm1 = poseidon.permutation(&input1);
            let perm2 = poseidon.permutation(&input1);
            let perm3 = poseidon.permutation(&input2);
            assert_eq!(perm1, perm2);
            assert_ne!(perm1, perm3);
        }
    }

    #[test]
    fn opt_equals_not_opt() {
        let poseidon = Poseidon::new(&POSEIDON_BABYBEAR_PARAMS);
        let t = poseidon.params.t;
        for _ in 0..TESTRUNS {
            let input: Vec<Scalar> = (0..t).map(|_| random_scalar()).collect();

            let perm1 = poseidon.permutation(&input);
            let perm2 = poseidon.permutation_not_opt(&input);
            assert_eq!(perm1, perm2);
        }
    }
}

#[cfg(test)]
#[allow(unused_imports)]
mod poseidon_tests_pallas {
    use super::*;
    use crate::fields::{pallas::FpPallas, utils::from_hex, utils::random_scalar};
    use crate::poseidon::poseidon_instance_pallas::{
        POSEIDON_PALLAS_3_PARAMS,
        POSEIDON_PALLAS_4_PARAMS,
        POSEIDON_PALLAS_8_PARAMS,
    };


    type Scalar = FpPallas;
    use ark_ff::UniformRand;

    static TESTRUNS: usize = 5;

    #[test]
    fn consistent_perm() {
        let instances = vec![
            Poseidon::new(&POSEIDON_PALLAS_3_PARAMS),
            Poseidon::new(&POSEIDON_PALLAS_4_PARAMS),
            Poseidon::new(&POSEIDON_PALLAS_8_PARAMS)
        ];
        for instance in instances {
            let t = instance.params.t;
            for _ in 0..TESTRUNS {
                let input1: Vec<Scalar> = (0..t).map(|_| random_scalar()).collect();

                let mut input2: Vec<Scalar>;
                loop {
                    input2 = (0..t).map(|_| random_scalar()).collect();
                    if input1 != input2 {
                        break;
                    }
                }

                let perm1 = instance.permutation(&input1);
                let perm2 = instance.permutation(&input1);
                let perm3 = instance.permutation(&input2);
                assert_eq!(perm1, perm2);
                assert_ne!(perm1, perm3);
            }
        }
    }

    #[test]
    fn kats() {
        let poseidon = Poseidon::new(&POSEIDON_PALLAS_3_PARAMS);
        let input: Vec<Scalar> = vec![Scalar::from(0), Scalar::from(1), Scalar::from(2)];
        let perm = poseidon.permutation(&input);
        assert_eq!(
            perm[0],
            from_hex("0x05d05891ddea8307b31c867776064d69da60aad0d14e99eee943c94282c30526")
        );
        assert_eq!(
            perm[1],
            from_hex("0x077a24130de1934e6f6cdbae982cc71615ddafce084fd06500336376d5d35796")
        );
        assert_eq!(
            perm[2],
            from_hex("0x1f6c3401c4f3fb94f1821e9ddaad8e173ed6c11ed4fb8bcf792f97ed8107e4ec")
        );
    }

    #[test]
    fn opt_equals_not_opt() {
        let instances = vec![
            Poseidon::new(&POSEIDON_PALLAS_3_PARAMS),
            Poseidon::new(&POSEIDON_PALLAS_4_PARAMS),
            Poseidon::new(&POSEIDON_PALLAS_8_PARAMS)
        ];
        for instance in instances {
            let t = instance.params.t;
            for _ in 0..TESTRUNS {
                let input: Vec<Scalar> = (0..t).map(|_| random_scalar()).collect();

                let perm1 = instance.permutation(&input);
                let perm2 = instance.permutation_not_opt(&input);
                assert_eq!(perm1, perm2);
            }
        }
    }
}

#[cfg(test)]
#[allow(unused_imports)]
mod poseidon_tests_vesta {
    use super::*;
    use crate::fields::{vesta::FpVesta, utils::from_hex, utils::random_scalar};
    use crate::poseidon::poseidon_instance_vesta::POSEIDON_VESTA_PARAMS;

    type Scalar = FpVesta;
    use ark_ff::UniformRand;

    static TESTRUNS: usize = 5;

    #[test]
    fn consistent_perm() {
        let poseidon = Poseidon::new(&POSEIDON_VESTA_PARAMS);
        let t = poseidon.params.t;
        for _ in 0..TESTRUNS {
            let input1: Vec<Scalar> = (0..t).map(|_| random_scalar()).collect();

            let mut input2: Vec<Scalar>;
            loop {
                input2 = (0..t).map(|_| random_scalar()).collect();
                if input1 != input2 {
                    break;
                }
            }

            let perm1 = poseidon.permutation(&input1);
            let perm2 = poseidon.permutation(&input1);
            let perm3 = poseidon.permutation(&input2);
            assert_eq!(perm1, perm2);
            assert_ne!(perm1, perm3);
        }
    }

    #[test]
    fn kats() {
        let poseidon = Poseidon::new(&POSEIDON_VESTA_PARAMS);
        let input: Vec<Scalar> = vec![Scalar::from(0), Scalar::from(1), Scalar::from(2)];
        let perm = poseidon.permutation(&input);
        assert_eq!(
            perm[0],
            from_hex("0x27a8d19d5d875f59e2439822e5f4ee80f2d113a522cab5d589e3cca52f66f868")
        );
        assert_eq!(
            perm[1],
            from_hex("0x33186d9b217a7ff8e125439b8ac8b3188870343d06db1691dec0215cc188e1de")
        );
        assert_eq!(
            perm[2],
            from_hex("0x1a4a058ffbfef2006edbd773a1f731d200be5ab8d631dfdc158f603b845caf91")
        );
    }

    #[test]
    fn opt_equals_not_opt() {
        let poseidon = Poseidon::new(&POSEIDON_VESTA_PARAMS);
        let t = poseidon.params.t;
        for _ in 0..TESTRUNS {
            let input: Vec<Scalar> = (0..t).map(|_| random_scalar()).collect();

            let perm1 = poseidon.permutation(&input);
            let perm2 = poseidon.permutation_not_opt(&input);
            assert_eq!(perm1, perm2);
        }
    }
}