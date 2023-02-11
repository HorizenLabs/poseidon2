use super::poseidon2_params::Poseidon2Params;
use crate::merkle_tree::merkle_tree_fp::MerkleTreeHash;
use ark_ff::PrimeField;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct Poseidon2<F: PrimeField> {
    pub(crate) params: Arc<Poseidon2Params<F>>,
}

impl<F: PrimeField> Poseidon2<F> {
    pub fn new(params: &Arc<Poseidon2Params<F>>) -> Self {
        Poseidon2 {
            params: Arc::clone(params),
        }
    }

    pub fn get_t(&self) -> usize {
        self.params.t
    }

    pub fn permutation(&self, input: &[F]) -> Vec<F> {
        let t = self.params.t;
        assert_eq!(input.len(), t);

        let mut current_state = input.to_owned();

        // Linear layer at beginning
        self.matmul_external(&mut current_state);

        for r in 0..self.params.rounds_f_beginning {
            current_state = self.add_rc(&current_state, &self.params.round_constants[r]);
            current_state = self.sbox(&current_state);
            self.matmul_external(&mut current_state);
        }

        let p_end = self.params.rounds_f_beginning + self.params.rounds_p;
        for r in self.params.rounds_f_beginning..p_end {
            current_state[0].add_assign(&self.params.round_constants[r][0]);
            current_state[0] = self.sbox_p(&current_state[0]);
            self.matmul_internal(&mut current_state, &self.params.mat_internal_diag_m_1);
        }
        
        for r in p_end..self.params.rounds {
            current_state = self.add_rc(&current_state, &self.params.round_constants[r]);
            current_state = self.sbox(&current_state);
            self.matmul_external(&mut current_state);
        }
        current_state
    }

    fn sbox(&self, input: &[F]) -> Vec<F> {
        input.iter().map(|el| self.sbox_p(el)).collect()
    }

    fn sbox_p(&self, input: &F) -> F {
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

    fn matmul_external(&self, input: &mut[F]) {
        let t = self.params.t;
        // let divisible_four = if t % 4 == 0 {
        //     true
        // } else {
        //     false
        // };

        match t {
            3 => {
                // Matrix circ(2, 1, 1)
                let mut sum = input[0];
                sum.add_assign(&input[1]);
                sum.add_assign(&input[2]);
                input[0].add_assign(&sum);
                input[1].add_assign(&sum);
                input[2].add_assign(&sum);
            }
            4 | 8 | 12 | 16 | 20 | 24 => {
                // Applying cheap 4x4 MDS matrix to each 4-element part of the state
                let t4 = t / 4;
                for i in 0..t4 {
                    let start_index = i * 4;
                    let mut t_0 = input[start_index];
                    t_0.add_assign(&input[start_index + 1]);
                    let mut t_1 = input[start_index + 2];
                    t_1.add_assign(&input[start_index + 3]);
                    let mut t_2 = input[start_index + 1];
                    t_2.double_in_place();
                    t_2.add_assign(&t_1);
                    let mut t_3 = input[start_index + 3];
                    t_3.double_in_place();
                    t_3.add_assign(&t_0);
                    let mut t_4 = t_1;
                    t_4.double_in_place();
                    t_4.double_in_place();
                    t_4.add_assign(&t_3);
                    let mut t_5 = t_0;
                    t_5.double_in_place();
                    t_5.double_in_place();
                    t_5.add_assign(&t_2);
                    let mut t_6 = t_3;
                    t_6.add_assign(&t_5);
                    let mut t_7 = t_2;
                    t_7.add_assign(&t_4);
                    input[start_index] = t_6;
                    input[start_index + 1] = t_5;
                    input[start_index + 2] = t_7;
                    input[start_index + 3] = t_4;
                }

                // Applying second cheap matrix
                let mut stored = [F::zero(); 4];
                for l in 0..4 {
                    stored[l] = input[l];
                    for j in 1..t4 {
                        stored[l].add_assign(&input[4 * j + l]);
                    }
                }
                for i in 0..input.len() {
                    input[i].add_assign(&stored[i % 4]);
                }
            }
            _ => {
                panic!()
            }
        }

    }

    fn matmul_internal(&self, input: &mut[F], mat_internal_diag_m_1: &[F]) {
        let t = self.params.t;

        match t {
            3 => {
                // [2, 1, 1]
                // [1, 2, 1]
                // [1, 1, 3]
                let mut sum = input[0];
                sum.add_assign(&input[1]);
                sum.add_assign(&input[2]);
                input[0].add_assign(&sum);
                input[1].add_assign(&sum);
                input[2].double_in_place();
                input[2].add_assign(&sum);
            }
            4 | 8 | 12 | 16 | 20 | 24 => {
                // Compute input sum
                let mut sum = input[0];
                input
                    .iter()
                    .skip(1)
                    .take(t-1)
                    .for_each(|el| sum.add_assign(el));
                // Add sum + diag entry * element to each element
                for i in 0..input.len() {
                    input[i].mul_assign(&mat_internal_diag_m_1[i]);
                    input[i].add_assign(&sum);
                }
            }
            _ => {
                panic!()
            }
        }
    }

    fn add_rc(&self, input: &[F], rc: &[F]) -> Vec<F> {
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

impl<F: PrimeField> MerkleTreeHash<F> for Poseidon2<F> {
    fn compress(&self, input: &[&F]) -> F {
        self.permutation(&[input[0].to_owned(), input[1].to_owned(), F::zero()])[0]
    }
}

#[allow(unused_imports)]
#[cfg(test)]
mod poseidon2_tests_goldilocks {
    use super::*;
    use crate::{fields::{goldilocks::FpGoldiLocks, utils::from_hex, utils::random_scalar}};
    use crate::poseidon2::poseidon2_instance_goldilocks::{
        POSEIDON2_GOLDILOCKS_8_PARAMS,
        POSEIDON2_GOLDILOCKS_12_PARAMS,
        POSEIDON2_GOLDILOCKS_16_PARAMS,
        POSEIDON2_GOLDILOCKS_20_PARAMS,
    };
    use std::convert::TryFrom;

    type Scalar = FpGoldiLocks;

    static TESTRUNS: usize = 5;

    #[test]
    fn consistent_perm() {
        let instances = vec![
            Poseidon2::new(&POSEIDON2_GOLDILOCKS_8_PARAMS),
            Poseidon2::new(&POSEIDON2_GOLDILOCKS_12_PARAMS),
            Poseidon2::new(&POSEIDON2_GOLDILOCKS_16_PARAMS),
            Poseidon2::new(&POSEIDON2_GOLDILOCKS_20_PARAMS),
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
        let poseidon2 = Poseidon2::new(&POSEIDON2_GOLDILOCKS_12_PARAMS);
        let mut input: Vec<Scalar> = vec![];
        for i in 0..poseidon2.params.t {
            input.push(Scalar::from(i as u64));
        }
        let perm = poseidon2.permutation(&input);
        assert_eq!(perm[0], from_hex("0xed3dbcc4ff1e8d33"));
        assert_eq!(perm[1], from_hex("0xfb85eac6ac91a150"));
        assert_eq!(perm[2], from_hex("0xd41e1e237ed3e2ef"));
        assert_eq!(perm[3], from_hex("0x5e289bf0a4c11897"));
        assert_eq!(perm[4], from_hex("0x4398b20f93e3ba6b"));
        assert_eq!(perm[5], from_hex("0x5659a48ffaf2901d"));
        assert_eq!(perm[6], from_hex("0xe44d81e89a88f8ae"));
        assert_eq!(perm[7], from_hex("0x08efdb285f8c3dbc"));
        assert_eq!(perm[8], from_hex("0x294ab7503297850e"));
        assert_eq!(perm[9], from_hex("0xa11c61f4870b9904"));
        assert_eq!(perm[10], from_hex("0xa6855c112cc08968"));
        assert_eq!(perm[11], from_hex("0x17c6d53d2fb3e8c1"));

    }
}

#[allow(unused_imports)]
#[cfg(test)]
mod poseidon2_tests_babybear {
    use super::*;
    use crate::{fields::{babybear::FpBabyBear, utils::from_hex, utils::random_scalar}};
    use crate::poseidon2::poseidon2_instance_babybear::{
        POSEIDON2_BABYBEAR_16_PARAMS,
        POSEIDON2_BABYBEAR_24_PARAMS,
    };
    use std::convert::TryFrom;

    type Scalar = FpBabyBear;

    static TESTRUNS: usize = 5;

    #[test]
    fn consistent_perm() {
        let instances = vec![
            Poseidon2::new(&POSEIDON2_BABYBEAR_16_PARAMS),
            Poseidon2::new(&POSEIDON2_BABYBEAR_24_PARAMS)
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
        let poseidon2 = Poseidon2::new(&POSEIDON2_BABYBEAR_24_PARAMS);
        let mut input: Vec<Scalar> = vec![];
        for i in 0..poseidon2.params.t {
            input.push(Scalar::from(i as u64));
        }
        let perm = poseidon2.permutation(&input);
        assert_eq!(perm[0], from_hex("0x7264c7b4"));
        assert_eq!(perm[1], from_hex("0x57f54b89"));
        assert_eq!(perm[2], from_hex("0x08e94bc9"));
        assert_eq!(perm[3], from_hex("0x4b05f544"));
        assert_eq!(perm[4], from_hex("0x68e1ccac"));
        assert_eq!(perm[5], from_hex("0x5aff4a53"));
        assert_eq!(perm[6], from_hex("0x685c3665"));
        assert_eq!(perm[7], from_hex("0x58b53bb8"));
        assert_eq!(perm[8], from_hex("0x0bdc0a03"));
        assert_eq!(perm[9], from_hex("0x0376e981"));
        assert_eq!(perm[10], from_hex("0x498926e6"));
        assert_eq!(perm[11], from_hex("0x3185cef3"));
        assert_eq!(perm[12], from_hex("0x17778434"));
        assert_eq!(perm[13], from_hex("0x551bbf78"));
        assert_eq!(perm[14], from_hex("0x46e0df61"));
        assert_eq!(perm[15], from_hex("0x04222917"));
        assert_eq!(perm[16], from_hex("0x09b8d147"));
        assert_eq!(perm[17], from_hex("0x6e0192ce"));
        assert_eq!(perm[18], from_hex("0x25215927"));
        assert_eq!(perm[19], from_hex("0x06b43026"));
        assert_eq!(perm[20], from_hex("0x093b34f3"));
        assert_eq!(perm[21], from_hex("0x6ed57a4c"));
        assert_eq!(perm[22], from_hex("0x73a7c7a3"));
        assert_eq!(perm[23], from_hex("0x511e49f5"));

    }
}

#[allow(unused_imports)]
#[cfg(test)]
mod poseidon2_tests_bls12 {
    use super::*;
    use crate::{fields::{bls12::FpBLS12, utils::from_hex, utils::random_scalar}};
    use crate::poseidon2::poseidon2_instance_bls12::{
        POSEIDON2_BLS_3_PARAMS,
        POSEIDON2_BLS_4_PARAMS,
        POSEIDON2_BLS_8_PARAMS,
    };
    use std::convert::TryFrom;

    type Scalar = FpBLS12;

    static TESTRUNS: usize = 5;

    #[test]
    fn consistent_perm() {
        let instances = vec![
            Poseidon2::new(&POSEIDON2_BLS_3_PARAMS),
            Poseidon2::new(&POSEIDON2_BLS_4_PARAMS),
            Poseidon2::new(&POSEIDON2_BLS_8_PARAMS)
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
        let poseidon2 = Poseidon2::new(&POSEIDON2_BLS_3_PARAMS);
        let mut input: Vec<Scalar> = vec![];
        for i in 0..poseidon2.params.t {
            input.push(Scalar::from(i as u64));
        }
        let perm = poseidon2.permutation(&input);
        assert_eq!(perm[0], from_hex("0x562af4b3710cdba6cea53e1f73325b21bb97ac810943b74d863d87163ee8042e"));
        assert_eq!(perm[1], from_hex("0x4674eba4cef166510c0d7a9ddf08cf813637bc2081e2c40c5047dce7ecdf2b95"));
        assert_eq!(perm[2], from_hex("0x0cf55ec35287dca6195eb6dd43e9ac1aba8857b4d3e4501be8bd8e9946a8dc54"));
    }
}

#[allow(unused_imports)]
#[cfg(test)]
mod poseidon2_tests_bn256 {
    use super::*;
    use crate::{fields::{bn256::FpBN256, utils::from_hex, utils::random_scalar}, poseidon2::poseidon2_instance_bn256::POSEIDON2_BN256_PARAMS};
    use std::convert::TryFrom;

    type Scalar = FpBN256;

    static TESTRUNS: usize = 5;

    #[test]
    fn consistent_perm() {
        let poseidon2 = Poseidon2::new(&POSEIDON2_BN256_PARAMS);
        let t = poseidon2.params.t;
        for _ in 0..TESTRUNS {
            let input1: Vec<Scalar> = (0..t).map(|_| random_scalar()).collect();

            let mut input2: Vec<Scalar>;
            loop {
                input2 = (0..t).map(|_| random_scalar()).collect();
                if input1 != input2 {
                    break;
                }
            }

            let perm1 = poseidon2.permutation(&input1);
            let perm2 = poseidon2.permutation(&input1);
            let perm3 = poseidon2.permutation(&input2);
            assert_eq!(perm1, perm2);
            assert_ne!(perm1, perm3);
        }
    }

    #[test]
    fn kats() {
        let poseidon2 = Poseidon2::new(&POSEIDON2_BN256_PARAMS);
        let mut input: Vec<Scalar> = vec![];
        for i in 0..poseidon2.params.t {
            input.push(Scalar::from(i as u64));
        }
        let perm = poseidon2.permutation(&input);
        assert_eq!(perm[0], from_hex("0x30610a447b7dec194697fb50786aa7421494bd64c221ba4d3b1af25fb07bd103"));
        assert_eq!(perm[1], from_hex("0x13f731d6ffbad391be22d2ac364151849e19fa38eced4e761bcd21dbdc600288"));
        assert_eq!(perm[2], from_hex("0x1433e2c8f68382c447c5c14b8b3df7cbfd9273dd655fe52f1357c27150da786f"));

    }
}

#[allow(unused_imports)]
#[cfg(test)]
mod poseidon2_tests_pallas {
    use super::*;
    use crate::{fields::{pallas::FpPallas, utils::from_hex, utils::random_scalar}};
    use crate::poseidon2::poseidon2_instance_pallas::{
        POSEIDON2_PALLAS_3_PARAMS,
        POSEIDON2_PALLAS_4_PARAMS,
        POSEIDON2_PALLAS_8_PARAMS,
    };
    use std::convert::TryFrom;

    type Scalar = FpPallas;

    static TESTRUNS: usize = 5;

    #[test]
    fn consistent_perm() {
        let instances = vec![
            Poseidon2::new(&POSEIDON2_PALLAS_3_PARAMS),
            Poseidon2::new(&POSEIDON2_PALLAS_4_PARAMS),
            Poseidon2::new(&POSEIDON2_PALLAS_8_PARAMS)
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
        let poseidon2 = Poseidon2::new(&POSEIDON2_PALLAS_3_PARAMS);
        let mut input: Vec<Scalar> = vec![];
        for i in 0..poseidon2.params.t {
            input.push(Scalar::from(i as u64));
        }
        let perm = poseidon2.permutation(&input);
        assert_eq!(perm[0], from_hex("0x1d6a0158834521fd53c786debb6456171b462f2ee3beb243bba0fa87e2fe7024"));
        assert_eq!(perm[1], from_hex("0x1d5caec6db09d78608971fd010944e5d81c036012031db1ae6f79b6886eff724"));
        assert_eq!(perm[2], from_hex("0x1622fe4b0dccd5e5393a31ba40e15c3e5a93a83dbe416f5e3d44f86998ce4878"));

    }
}

#[allow(unused_imports)]
#[cfg(test)]
mod poseidon2_tests_vesta {
    use super::*;
    use crate::{fields::{vesta::FpVesta, utils::from_hex, utils::random_scalar}, poseidon2::poseidon2_instance_vesta::POSEIDON2_VESTA_PARAMS};
    use std::convert::TryFrom;

    type Scalar = FpVesta;

    static TESTRUNS: usize = 5;

    #[test]
    fn consistent_perm() {
        let poseidon2 = Poseidon2::new(&POSEIDON2_VESTA_PARAMS);
        let t = poseidon2.params.t;
        for _ in 0..TESTRUNS {
            let input1: Vec<Scalar> = (0..t).map(|_| random_scalar()).collect();

            let mut input2: Vec<Scalar>;
            loop {
                input2 = (0..t).map(|_| random_scalar()).collect();
                if input1 != input2 {
                    break;
                }
            }

            let perm1 = poseidon2.permutation(&input1);
            let perm2 = poseidon2.permutation(&input1);
            let perm3 = poseidon2.permutation(&input2);
            assert_eq!(perm1, perm2);
            assert_ne!(perm1, perm3);
        }
    }

    #[test]
    fn kats() {
        let poseidon2 = Poseidon2::new(&POSEIDON2_VESTA_PARAMS);
        let mut input: Vec<Scalar> = vec![];
        for i in 0..poseidon2.params.t {
            input.push(Scalar::from(i as u64));
        }
        let perm = poseidon2.permutation(&input);
        assert_eq!(perm[0], from_hex("0x2dbc40552a2b2a785f63489278361134609229297aafcde2bc8958d45546966f"));
        assert_eq!(perm[1], from_hex("0x373efbf666c89d0653612fadf73ab1db47df46c38e86a959d642801a7c4b0201"));
        assert_eq!(perm[2], from_hex("0x09cde65715aa5d7a6a5bc5feedd549489f805de1c7951e2135b8135ef1112f59"));

    }
}