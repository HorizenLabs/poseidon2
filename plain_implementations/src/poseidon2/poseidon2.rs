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
        current_state = self.add_rc(&current_state, &self.params.opt_round_constants[0]);
        for r in 1..(self.params.rounds_p + 1) {
            current_state[0] = self.sbox_p(&current_state[0]);
            current_state[0].add_assign(&self.params.opt_round_constants[r][0]);
            self.matmul_internal(&mut current_state, &self.params.mat_internal_diag_m_1);
        }
        
        for r in p_end..self.params.rounds {
            current_state = self.add_rc(&current_state, &self.params.round_constants[r]);
            current_state = self.sbox(&current_state);
            self.matmul_external(&mut current_state);
        }
        current_state
    }

    pub fn permutation_not_opt(&self, input: &[F]) -> Vec<F> {
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
            current_state = self.add_rc(&current_state, &self.params.round_constants[r]);
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
        assert_eq!(perm[0], from_hex("0x12a0e34a20879ec7"));
        assert_eq!(perm[1], from_hex("0xdff531b65bcd2770"));
        assert_eq!(perm[2], from_hex("0xacfde64ace0bb075"));
        assert_eq!(perm[3], from_hex("0x0d9044b55d697e2c"));
        assert_eq!(perm[4], from_hex("0x3afffd14aeaab03d"));
        assert_eq!(perm[5], from_hex("0x49f51fc35517972a"));
        assert_eq!(perm[6], from_hex("0x74da7d0508892cc7"));
        assert_eq!(perm[7], from_hex("0x04cdca75b5646fe2"));
        assert_eq!(perm[8], from_hex("0xb110bcea2fe730c1"));
        assert_eq!(perm[9], from_hex("0x9558a7f7af657d1e"));
        assert_eq!(perm[10], from_hex("0x9751e0cb1ee16c17"));
        assert_eq!(perm[11], from_hex("0xd2574bee1659fee7"));

    }

    #[test]
    fn opt_equals_not_opt() {
        let instances = vec![
            Poseidon2::new(&POSEIDON2_GOLDILOCKS_8_PARAMS),
            Poseidon2::new(&POSEIDON2_GOLDILOCKS_12_PARAMS),
            Poseidon2::new(&POSEIDON2_GOLDILOCKS_16_PARAMS),
            Poseidon2::new(&POSEIDON2_GOLDILOCKS_20_PARAMS),
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
        assert_eq!(perm[0], from_hex("0x2ec08c0d"));
        assert_eq!(perm[1], from_hex("0x5fe607c5"));
        assert_eq!(perm[2], from_hex("0x10ab945e"));
        assert_eq!(perm[3], from_hex("0x1b8973cc"));
        assert_eq!(perm[4], from_hex("0x15870998"));
        assert_eq!(perm[5], from_hex("0x66fce31d"));
        assert_eq!(perm[6], from_hex("0x3f38ea43"));
        assert_eq!(perm[7], from_hex("0x6b9c3705"));
        assert_eq!(perm[8], from_hex("0x02f07a0b"));
        assert_eq!(perm[9], from_hex("0x4b052f69"));
        assert_eq!(perm[10], from_hex("0x65bfca0d"));
        assert_eq!(perm[11], from_hex("0x3a4baba8"));
        assert_eq!(perm[12], from_hex("0x71d9b602"));
        assert_eq!(perm[13], from_hex("0x46335095"));
        assert_eq!(perm[14], from_hex("0x0c68c3c2"));
        assert_eq!(perm[15], from_hex("0x4133626e"));
        assert_eq!(perm[16], from_hex("0x109e0b39"));
        assert_eq!(perm[17], from_hex("0x452ae1e0"));
        assert_eq!(perm[18], from_hex("0x6d63b8f5"));
        assert_eq!(perm[19], from_hex("0x1f2fc257"));
        assert_eq!(perm[20], from_hex("0x12894b7d"));
        assert_eq!(perm[21], from_hex("0x4a03d9b8"));
        assert_eq!(perm[22], from_hex("0x653f5994"));
        assert_eq!(perm[23], from_hex("0x1ba2f443"));

    }

    #[test]
    fn opt_equals_not_opt() {
        let instances = vec![
            Poseidon2::new(&POSEIDON2_BABYBEAR_16_PARAMS),
            Poseidon2::new(&POSEIDON2_BABYBEAR_24_PARAMS)
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
        assert_eq!(perm[0], from_hex("0x0c07cd926a6b1cdc4671b5063f1d958fd9ed157ce0afc5a3a07296bd44a03f71"));
        assert_eq!(perm[1], from_hex("0x719a8dfa8df976c1cd0d02e30eb7ba48eabbef2d94c0ea3a864238466812a314"));
        assert_eq!(perm[2], from_hex("0x052937bf500adbb26dc1a5ae50b37cbd3bbb7c9863be8976506199dd5179d6da"));
    }

    #[test]
    fn opt_equals_not_opt() {
        let instances = vec![
            Poseidon2::new(&POSEIDON2_BLS_3_PARAMS),
            Poseidon2::new(&POSEIDON2_BLS_4_PARAMS),
            Poseidon2::new(&POSEIDON2_BLS_8_PARAMS)
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
        assert_eq!(perm[0], from_hex("0x12df045f8784b551534ecbfb3d039c6840577a508f67fcdf2624c77f189f828b"));
        assert_eq!(perm[1], from_hex("0x209b2c3f875fff41ed81576b77089a2cfa30d3c840ac66307cb7cd9292a5018f"));
        assert_eq!(perm[2], from_hex("0x16414ee2d17ef3598e807b2a5b732d25f2b1fc821126b21fce447b2ee0fc0d0c"));

    }

    #[test]
    fn opt_equals_not_opt() {
        let poseidon2 = Poseidon2::new(&POSEIDON2_BN256_PARAMS);
        let t = poseidon2.params.t;
        for _ in 0..TESTRUNS {
            let input: Vec<Scalar> = (0..t).map(|_| random_scalar()).collect();

            let perm1 = poseidon2.permutation(&input);
            let perm2 = poseidon2.permutation_not_opt(&input);
            assert_eq!(perm1, perm2);
        }
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
        assert_eq!(perm[0], from_hex("0x0db342985726e25dd266eed3e5acb59ac009d9489ee2b3dff31425378ff8c207"));
        assert_eq!(perm[1], from_hex("0x2faba4e66cc58fc3ac1a35d1575b71188eaf53b297ae96edbc88257f73250690"));
        assert_eq!(perm[2], from_hex("0x1a002ed5e9e105dc7614a5220c6b98e33a83e494aeeae51c89c68222523f8772"));

    }

    #[test]
    fn opt_equals_not_opt() {
        let instances = vec![
            Poseidon2::new(&POSEIDON2_PALLAS_3_PARAMS),
            Poseidon2::new(&POSEIDON2_PALLAS_4_PARAMS),
            Poseidon2::new(&POSEIDON2_PALLAS_8_PARAMS)
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
        assert_eq!(perm[0], from_hex("0x1d440fc641aac892570b84e71acf7c5318dd644fee0e422da954035fdf23a41d"));
        assert_eq!(perm[1], from_hex("0x2fa30f04971aafde79157ab0e08a817407d420a991c43a375945f0530db32ef5"));
        assert_eq!(perm[2], from_hex("0x2ec7abde6fc2af15997f510e0995d0c39c9136099509b392bc3de44052349cc4"));

    }

    #[test]
    fn opt_equals_not_opt() {
        let poseidon2 = Poseidon2::new(&POSEIDON2_VESTA_PARAMS);
        let t = poseidon2.params.t;
        for _ in 0..TESTRUNS {
            let input: Vec<Scalar> = (0..t).map(|_| random_scalar()).collect();

            let perm1 = poseidon2.permutation(&input);
            let perm2 = poseidon2.permutation_not_opt(&input);
            assert_eq!(perm1, perm2);
        }
    }
}