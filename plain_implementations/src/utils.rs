// use std::cmp::min;

use ark_ff::PrimeField;

// pub fn from_u64<F: PrimeField>(val: u64) -> F {
//     F::from_repr(F::Repr::from(val)).unwrap()
// }

// guassian elimination
pub fn mat_inverse<F: PrimeField>(mat: &[Vec<F>]) -> Vec<Vec<F>> {
    let n = mat.len();
    assert!(mat[0].len() == n);

    let mut m = mat.to_owned();
    let mut inv = vec![vec![F::zero(); n]; n];
    for (i, invi) in inv.iter_mut().enumerate() {
        invi[i] = F::one();
    }

    // upper triangle
    for row in 0..n {
        for j in 0..row {
            // subtract from these rows
            let el = m[row][j];
            for col in 0..n {
                // do subtraction for each col
                if col < j {
                    m[row][col] = F::zero();
                } else {
                    let mut tmp = m[j][col];
                    tmp.mul_assign(&el);
                    m[row][col].sub_assign(&tmp);
                }
                if col > row {
                    inv[row][col] = F::zero();
                } else {
                    let mut tmp = inv[j][col];
                    tmp.mul_assign(&el);
                    inv[row][col].sub_assign(&tmp);
                }
            }
        }
        // make 1 in diag
        let el_inv = m[row][row].inverse().unwrap();
        for col in 0..n {
            match col.cmp(&row) {
                std::cmp::Ordering::Less => inv[row][col].mul_assign(&el_inv),
                std::cmp::Ordering::Equal => {
                    m[row][col] = F::one();
                    inv[row][col].mul_assign(&el_inv)
                }
                std::cmp::Ordering::Greater => m[row][col].mul_assign(&el_inv),
            }
        }
    }

    // upper triangle
    for row in (0..n).rev() {
        for j in (row + 1..n).rev() {
            // subtract from these rows
            let el = m[row][j];
            for col in 0..n {
                // do subtraction for each col

                #[cfg(debug_assertions)]
                {
                    if col >= j {
                        m[row][col] = F::zero();
                    }
                }
                let mut tmp = inv[j][col];
                tmp.mul_assign(&el);
                inv[row][col].sub_assign(&tmp);
            }
        }
    }

    #[cfg(debug_assertions)]
    {
        for (row, mrow) in m.iter().enumerate() {
            for (col, v) in mrow.iter().enumerate() {
                if row == col {
                    debug_assert!(*v == F::one());
                } else {
                    debug_assert!(*v == F::zero());
                }
            }
        }
    }

    inv
}

pub fn mat_transpose<F: PrimeField>(mat: &[Vec<F>]) -> Vec<Vec<F>> {
    let rows = mat.len();
    let cols = mat[0].len();
    let mut transpose = vec![vec![F::zero(); rows]; cols];

    for (row, matrow) in mat.iter().enumerate() {
        for col in 0..cols {
            transpose[col][row] = matrow[col];
        }
    }
    transpose
}

// pub fn random_scalar_rng<F: PrimeField, R: Rng>(allow_zero: bool, rng: &mut R) -> F {
//     loop {
//         let s = F::rand(rng);
//         if allow_zero || s != F::zero() {
//             return s;
//         }
//     }
// }

// pub fn random_scalar<F: PrimeField>(allow_zero: bool) -> F {
//     loop {
//         let s = F::rand(&mut thread_rng());
//         if allow_zero || s != F::zero() {
//             return s;
//         }
//     }
// }

// pub fn into_limbs<F: PrimeField>(val: &F) -> Vec<u64> {
//     val.into_repr().as_ref().to_owned()
// }

// pub fn from_limbs<F: PrimeField>(repr: &[u64]) -> F {
//     let mut tmp = F::Repr::default();
//     tmp.as_mut().copy_from_slice(repr);
//     F::from_repr(tmp).unwrap()
// }

//-----------------------------------------------------------------------------
// mod_inverse
//-----------------------------------------------------------------------------
// pub fn mod_inverse<F: PrimeField>(val: u16, modulus: &F::Repr) -> F::Repr {
//     if val == 0 {
//         panic!("0 has no inverse!");
//     }

//     let mut m = val;
//     let mut tmp_v = modulus.to_owned();

//     let (q, tmp) = divide_long::<F>(&tmp_v, m);
//     let mut v = m;
//     m = tmp;
//     let mut a = q;
//     let mut a_neg = true;
//     let mut prev_a = F::Repr::from(1);
//     let mut prev_a_neg = false;

//     while m != 0 {
//         let q = v / m;
//         let tmp = v % m;
//         v = m;
//         m = tmp;

//         let tmp_a = a;
//         let tmp_a_neg = a_neg;

//         let (qa, _) = mul_by_single_word::<F>(&a, q as u64);
//         if a_neg != prev_a_neg {
//             a = prev_a;
//             a.add_nocarry(&qa);
//             a_neg = prev_a_neg;
//         } else if prev_a > qa {
//             a = prev_a;
//             a.sub_noborrow(&qa);
//             a_neg = prev_a_neg;
//         } else {
//             a = qa;
//             a.sub_noborrow(&prev_a);
//             a_neg = !a_neg;
//         }

//         prev_a = tmp_a;
//         prev_a_neg = tmp_a_neg;
//     }

//     if v != 1 {
//         panic!("{} has no inverse!", val);
//     }

//     if prev_a_neg {
//         tmp_v.sub_noborrow(&prev_a);
//         tmp_v
//     } else {
//         prev_a
//     }
// }

// #[cfg(test)]
// fn add<F: PrimeField>(lhs: &F::Repr, rhs: &F::Repr) -> (F::Repr, u64) {
//     let mut c = 0;
//     let mut res = F::Repr::default();
//     res.as_mut()
//         .iter_mut()
//         .zip(lhs.as_ref().iter())
//         .zip(rhs.as_ref().iter())
//         .for_each(|((res, l), r)| {
//             let tmp = l.overflowing_add(*r);
//             let c_ = tmp.1 as u64;
//             let tmp = tmp.0.overflowing_add(c);
//             *res = tmp.0;
//             c = c_ + tmp.1 as u64;
//         });

//     (res, c)
// }

// #[cfg(test)]
// fn subtract<F: PrimeField>(c: u64, lhs: &F::Repr, rhs: &F::Repr) -> (F::Repr, u64) {
//     if lhs >= rhs {
//         let mut res = lhs.to_owned();
//         res.sub_noborrow(rhs);
//         (res, c)
//     } else {
//         let mut res = rhs.to_owned();
//         res.sub_noborrow(lhs);
//         res.as_mut().iter_mut().for_each(|r| *r = !*r);
//         res = add_single_word::<F>(&res, 1).0;
//         (res, c - 1)
//     }
// }

// -----------------------------------------------------------------------------
// Credit shamatar
//-----------------------------------------------------------------------------

// pub fn mul_by_single_word<F: PrimeField>(u: &F::Repr, w: u64) -> (F::Repr, u64) {
//     let mut res = F::Repr::default();

//     let u_ref = u.as_ref();
//     let res_mut = res.as_mut();

//     let w_ = w as u128;

//     let mut tmp = (u_ref[0] as u128) * w_;
//     res_mut[0] = tmp as u64;
//     res_mut
//         .iter_mut()
//         .zip(u_ref.iter())
//         .skip(1)
//         .for_each(|(r, u_)| {
//             tmp = (*u_ as u128) * w_ + (tmp >> 64);
//             *r = tmp as u64;
//         });
//     (res, (tmp >> 64) as u64)
// }

// pub fn add_single_word<F: PrimeField>(u: &F::Repr, w: u64) -> (F::Repr, u64) {
//     let mut res = F::Repr::default();

//     let mut of = w;
//     res.as_mut()
//         .iter_mut()
//         .zip(u.as_ref().iter())
//         .for_each(|(r, u)| {
//             let (tmp, o) = u.overflowing_add(of);
//             *r = tmp;
//             of = o as u64;
//         });

//     (res, of)
// }

// #[inline(always)]
// fn div_mod_word_by_short(hi: u64, lo: u64, y: u16) -> (u64, u64) {
//     let t = ((hi as u128) << 64) + lo as u128;
//     let q = (t / (y as u128)) as u64;
//     let r = (t % (y as u128)) as u64;

//     (q, r)
// }

// pub fn divide_long<F: PrimeField>(a: &F::Repr, divisor: u16) -> (F::Repr, u16) {
//     let mut result = F::Repr::default();

//     let a_ref = a.as_ref();
//     let result_mut = result.as_mut();

//     let len = a_ref.len();

//     result_mut[len - 1] = a_ref[len - 1] / (divisor as u64);
//     let mut r = a.as_ref()[len - 1] % (divisor as u64);

//     result_mut
//         .iter_mut()
//         .zip(a_ref.iter())
//         .rev()
//         .skip(1)
//         .for_each(|(res, a_)| {
//             let (q, m) = div_mod_word_by_short(r, *a_, divisor);
//             *res = q;
//             r = m;
//         });

//     (result, r as u16)
// }

// #[cfg(test)]
// mod utils_test_bn {
//     use super::*;
//     use bellman_ce::pairing::{
//         bn256,
//         ff::{Field, PrimeField},
//         to_hex,
//     };

//     static TESTRUNS: usize = 5;

//     type Scalar = bn256::Fr;

//     #[test]
//     fn random() {
//         let rands: Vec<Scalar> = (0..TESTRUNS).map(|_| random_scalar(true)).collect();
//         for i in 0..TESTRUNS {
//             for j in i + 1..TESTRUNS {
//                 assert_ne!(rands[i], rands[j]);
//             }
//         }
//     }

//     #[test]
//     fn from_u64() {
//         let ten = super::from_u64::<Scalar>(10);
//         assert_eq!(
//             to_hex(&ten),
//             "000000000000000000000000000000000000000000000000000000000000000a"
//         )
//     }

//     #[test]
//     fn limbs() {
//         let ten = super::from_u64::<Scalar>(10);
//         let ten_limbs = [10, 0, 0, 0];
//         assert_eq!(into_limbs::<Scalar>(&ten), ten_limbs);
//         assert_eq!(ten, from_limbs::<Scalar>(&ten_limbs));
//         let input: Scalar = random_scalar(true);

//         for _ in 0..TESTRUNS {
//             assert_eq!(input, from_limbs::<Scalar>(&into_limbs::<Scalar>(&input)));
//         }
//     }

//     #[test]
//     fn div_mod_multiply_add() {
//         let mut rng = thread_rng();

//         // KAT
//         let ten = super::from_u64::<Scalar>(10);
//         let ten_repr = ten.into_repr();
//         let div = 3;

//         let (res, m) = divide_long::<Scalar>(&ten_repr, div);

//         assert_eq!(m, 1);
//         assert_eq!(
//             to_hex(&Scalar::from_repr(res).unwrap()),
//             "0000000000000000000000000000000000000000000000000000000000000003"
//         );

//         let (tmp, _) = mul_by_single_word::<Scalar>(&res, div as u64);
//         let (tmp, _) = add_single_word::<Scalar>(&tmp, m as u64);
//         assert_eq!(Scalar::from_repr(tmp).unwrap(), ten);

//         let mut res = Scalar::from_repr(res).unwrap();
//         let div = super::from_u64::<Scalar>(div as u64);
//         let m = super::from_u64::<Scalar>(m as u64);

//         res.mul_assign(&div);
//         res.add_assign(&m);

//         assert_eq!(ten, res);

//         // rand tests
//         for _ in 0..TESTRUNS {
//             let input: Scalar = random_scalar_rng(true, &mut rng);
//             let mut div = rng.gen::<u16>();
//             if div == 0 {
//                 div = 1;
//             }
//             let (res, m) = divide_long::<Scalar>(&input.into_repr(), div);

//             let (tmp, _) = mul_by_single_word::<Scalar>(&res, div as u64);
//             let (tmp, _) = add_single_word::<Scalar>(&tmp, m as u64);
//             assert_eq!(Scalar::from_repr(tmp).unwrap(), input);
//         }
//     }

//     #[test]
//     fn add_sub() {
//         for _ in 0..TESTRUNS {
//             let input1: Scalar = random_scalar(true);
//             let input2: Scalar = random_scalar(true);

//             let repr1 = input1.into_repr();
//             let repr2 = input2.into_repr();

//             let tmp = add::<Scalar>(&repr1, &repr2);

//             let res1 = subtract::<Scalar>(tmp.1, &tmp.0, &repr1);
//             let res2 = subtract::<Scalar>(tmp.1, &tmp.0, &repr2);

//             assert_eq!(res1.1, 0);
//             assert_eq!(res2.1, 0);
//             assert_eq!(Scalar::from_repr(res1.0).unwrap(), input2);
//             assert_eq!(Scalar::from_repr(res2.0).unwrap(), input1);
//         }
//     }

//     #[test]
//     fn mod_inverse_test() {
//         let d = 5;
//         let mut p_1 = Scalar::one();
//         p_1.negate();
//         let p_1 = p_1.into_repr();
//         let d_inv = mod_inverse::<Scalar>(d as u16, &p_1);

//         let (r, c) = mul_by_single_word::<Scalar>(&d_inv, d);
//         let mut r = r;
//         if c != 0 {
//             let mut c = c;
//             while c != 0 {
//                 let tmp = subtract::<Scalar>(c, &r, &p_1);
//                 c = tmp.1;
//                 r = tmp.0;
//             }
//         }
//         while r >= p_1 {
//             r.sub_noborrow(&p_1);
//         }
//         assert_eq!(Scalar::from_repr(r).unwrap(), Scalar::one());
//     }
// }

// #[cfg(test)]
// mod utils_test_bls {
//     use super::*;
//     use bellman_ce::pairing::{
//         bls12_381,
//         ff::{Field, PrimeField},
//         to_hex,
//     };

//     static TESTRUNS: usize = 5;

//     type Scalar = bls12_381::Fr;

//     #[test]
//     fn random() {
//         let rands: Vec<Scalar> = (0..TESTRUNS).map(|_| random_scalar(true)).collect();
//         for i in 0..TESTRUNS {
//             for j in i + 1..TESTRUNS {
//                 assert_ne!(rands[i], rands[j]);
//             }
//         }
//     }

//     #[test]
//     fn from_u64() {
//         let ten = super::from_u64::<Scalar>(10);
//         assert_eq!(
//             to_hex(&ten),
//             "000000000000000000000000000000000000000000000000000000000000000a"
//         )
//     }

//     #[test]
//     fn limbs() {
//         let ten = super::from_u64::<Scalar>(10);
//         let ten_limbs = [10, 0, 0, 0];
//         assert_eq!(into_limbs::<Scalar>(&ten), ten_limbs);
//         assert_eq!(ten, from_limbs::<Scalar>(&ten_limbs));
//         let input: Scalar = random_scalar(true);

//         for _ in 0..TESTRUNS {
//             assert_eq!(input, from_limbs::<Scalar>(&into_limbs::<Scalar>(&input)));
//         }
//     }

//     #[test]
//     fn div_mod_multiply_add() {
//         let mut rng = thread_rng();

//         // KAT
//         let ten = super::from_u64::<Scalar>(10);
//         let ten_repr = ten.into_repr();
//         let div = 3;

//         let (res, m) = divide_long::<Scalar>(&ten_repr, div);

//         assert_eq!(m, 1);
//         assert_eq!(
//             to_hex(&Scalar::from_repr(res).unwrap()),
//             "0000000000000000000000000000000000000000000000000000000000000003"
//         );

//         let (tmp, _) = mul_by_single_word::<Scalar>(&res, div as u64);
//         let (tmp, _) = add_single_word::<Scalar>(&tmp, m as u64);
//         assert_eq!(Scalar::from_repr(tmp).unwrap(), ten);

//         let mut res = Scalar::from_repr(res).unwrap();
//         let div = super::from_u64::<Scalar>(div as u64);
//         let m = super::from_u64::<Scalar>(m as u64);

//         res.mul_assign(&div);
//         res.add_assign(&m);

//         assert_eq!(ten, res);

//         // rand tests
//         for _ in 0..TESTRUNS {
//             let input: Scalar = random_scalar_rng(true, &mut rng);
//             let mut div = rng.gen::<u16>();
//             if div == 0 {
//                 div = 1;
//             }
//             let (res, m) = divide_long::<Scalar>(&input.into_repr(), div);

//             let (tmp, _) = mul_by_single_word::<Scalar>(&res, div as u64);
//             let (tmp, _) = add_single_word::<Scalar>(&tmp, m as u64);
//             assert_eq!(Scalar::from_repr(tmp).unwrap(), input);
//         }
//     }

//     #[test]
//     fn add_sub() {
//         for _ in 0..TESTRUNS {
//             let input1: Scalar = random_scalar(true);
//             let input2: Scalar = random_scalar(true);

//             let repr1 = input1.into_repr();
//             let repr2 = input2.into_repr();

//             let tmp = add::<Scalar>(&repr1, &repr2);

//             let res1 = subtract::<Scalar>(tmp.1, &tmp.0, &repr1);
//             let res2 = subtract::<Scalar>(tmp.1, &tmp.0, &repr2);

//             assert_eq!(res1.1, 0);
//             assert_eq!(res2.1, 0);
//             assert_eq!(Scalar::from_repr(res1.0).unwrap(), input2);
//             assert_eq!(Scalar::from_repr(res2.0).unwrap(), input1);
//         }
//     }

//     #[test]
//     fn mod_inverse_test() {
//         let d = 5;
//         let mut p_1 = Scalar::one();
//         p_1.negate();
//         let p_1 = p_1.into_repr();
//         let d_inv = mod_inverse::<Scalar>(d as u16, &p_1);

//         let (r, c) = mul_by_single_word::<Scalar>(&d_inv, d);
//         let mut r = r;
//         if c != 0 {
//             let mut c = c;
//             while c != 0 {
//                 let tmp = subtract::<Scalar>(c, &r, &p_1);
//                 c = tmp.1;
//                 r = tmp.0;
//             }
//         }

//         while r >= p_1 {
//             r.sub_noborrow(&p_1);
//         }
//         assert_eq!(Scalar::from_repr(r).unwrap(), Scalar::one());
//     }
// }
