use p3_air::{Air, AirBuilder, BaseAir};
use p3_field::{AbstractField, Field};
use p3_matrix:: Matrix;
use p3_matrix::dense::RowMajorMatrix;
use std::ops::{Add, Mul};
// use ark_ff::fields::models::fp::{Fp64, MontBackend, MontConfig};
// use ark_poly::{polynomial::univariate::DensePolynomial, DenseUVPolynomial};
// use ark_poly::Polynomial;
// use ark_ff::PrimeField;
use crate::params::N;

// #[derive(MontConfig)]
// #[modulus = "1085276161"]
// #[generator = "11"]
// pub struct FqConfig;
// pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

// Define AIR constraint
pub struct PolyMulAir {
	pub a: Vec<u32>,
	pub b: Vec<u32>,
    pub modulus: u32
}

/*
Polynomial Multiplication Air
Input:
- a = a[0] + a[1] * X + ... + a[N-1] * X^{N-1}
- b = b[0] + b[1] * X + ... + b[N-1] * X^{N-1}
Output:
- out = out[0] + out[1] * X + ... + out[2N-2] * X^{2N-2}

Note:
- PolyMulAir does not have a state transition. Values required for constraints are all stored in one row.
- While output polynomial `out` is calculated manually by generate_polymul_trace(),
we prove that this multiplication was done correctly, by enforcing a constraint such that a(x)*b(x) === out(x)  at x = [0..2N-1) based on Lagrange polynomial interpolation.
*/
impl<F: Field> BaseAir<F> for PolyMulAir {
    // Air Table looks like this
    // row:[     a: N     ][     b: N     ][               out(x): 2N-1               ]
    //     ^----------inputs------------- ^^---calculated by generate_polymul_trace---^
    //     [0........................................................................0]
    //     [0........................................................................0]
    //     [0........................................................................0]
    fn width(&self) -> usize {
         4*N-1
    }
}

fn mod_exp(mut base: u64, mut exp: u64, modulus: u64) -> u64 {
    if modulus == 1 {
        return 0;
    }
    let mut result = 1;
    base %= modulus; // Initial reduction of base

    // perform exponentiation by iterating exponents in binary representation from the LSB to MSB
    while exp > 0 {
        // when the bit is 1: base * result
        if exp % 2 == 1 {
            result = base * result % modulus;
        }
        // right shift exponent to the right by 1
        exp >>= 1;
        base = base * base % modulus;
    }
    result
}

// Define constraints
impl<AB: AirBuilder> Air<AB> for PolyMulAir {
    fn eval(&self, builder: &mut AB) {

        let main = builder.main();
        let row = main.row_slice(0);

        // Enforce self.a and self.b as 2 input polynomials
		for i in 0..N {
            builder.when_first_row().assert_eq(row[i], AB::Expr::from_canonical_u32(self.a[i]));
			builder.when_first_row().assert_eq(row[i+N], AB::Expr::from_canonical_u32(self.b[i]));
		}

        let mut a_eval: Vec<<AB as AirBuilder>::Expr> = Vec::with_capacity(2*N-1);
        let mut b_eval: Vec<<AB as AirBuilder>::Expr> = Vec::with_capacity(2*N-1);
        let mut out_eval: Vec<<AB as AirBuilder>::Expr> = Vec::with_capacity(2*N-1);

        // Evaluate 2 input polynomial a(x) and b(x) at x = [0..2N-1)
        // a = a[0] + a[1] * X + ... + a[N-1] * X^{N-1}
        // when x = 0, a_eval[0] = a[0] + a[1]*0 + a[2]*0^2 + ... + a[N-1] * 0^{N-1}
        // when x = 1, a_eval[1] = a[0] + a[1]*1 + a[2]*1^2 + ... + a[N-1] * 1^{N-1}
        // ...
        // when x = 2N-1, a_eval[2N-1] = a[0] + a[1]*(2N-1) + ... + a[N-1] * (2N-1)^{N-1}
        for i in 0..2*N-1 {
            a_eval.push(AB::Expr::zero());
            b_eval.push(AB::Expr::zero());
            for j in 0..N {
                let power = AB::Expr::from_canonical_u64(mod_exp(i as u64, j as u64, self.modulus as u64));

                let _ = a_eval[i].clone().add(row[j].mul(power.clone()));
                let _ = b_eval[i].clone().add(row[j+N].mul(power));
            }
        }

        // Evaluate output polynomial out(x) at x = [0..2N-1)
        for i in 0..2*N-1 {
            out_eval.push(AB::Expr::zero());
            for j in 0..2*N-1 {
                let power = AB::Expr::from_canonical_u64(mod_exp(i as u64, j as u64, self.modulus as u64));

                let _ = out_eval[i].clone().add(row[j+2*N].mul(power));
            }
        }

       // Enforce a[x] * b[x] === out[x] at x = [0..2N-1)
       // TODO: add non-native modular reduction
       // currently this is under-constrained
        for i in 0..2*N-1 {
            builder.assert_eq(a_eval[i].clone().mul(b_eval[i].clone()), out_eval[i].clone());
        }
    }

}

// Define a function to generate execution trace
pub fn generate_polymul_trace<F: Field>(a:Vec<u32>, b:Vec<u32>, modulus: u32) -> RowMajorMatrix<F> {
    let mut values: Vec<F>= Vec::with_capacity(4 * (4*N-1)); // 4 is the minimum number of rows required

	// Assign input polynomials to values vector
	for i in 0..N {
		values.push(F::from_canonical_u32(a[i]));
	}
	for i in 0..N {
		values.push(F::from_canonical_u32(b[i]));
	}

    let mut out:Vec<u128> = Vec::with_capacity(2*N-1);

	// Multiply the 2 polynomials manually and assign coefficients to values vector
    // Temporarily using u128 for intermediate values to avoid overflow
	for i in 0..2*N-1 {
        if i < N {
            // a's index increases from 0 to i, b's index decreases from i to 0
            // ex. N = 3 where N is the number of coefficients
            // when i = 0, a[0] * b[0]
            // when i = 1, a[0] * b[1] + a[1] * b[0]
            // when i = 2, a[0] * b[2] + a[1] * b[1] + a[2] * b[0]
            out.push(0);
            for a_idx in 0..i+1 {
                let b_idx = i - a_idx;
                out[i] += a[a_idx] as u128 * b[b_idx] as u128 % modulus as u128;
            }

        } else {
            // a's index increases from i-N+1 to N-1, which is the highest degree of input polynomial, b's index decreases from N-1 to i-(N-1)
            // ex. N = 3 where N is the number of coefficients
            // when i = 3, a[1] * b[2] + a[2] * b[1]
            // when i = 4, a[2] * b[2]
            out.push(0);
            for a_idx in i-N+1..N {
                let b_idx = i - a_idx;
                out[i] += a[a_idx] as u128 * b[b_idx] as u128 % modulus as u128;
            }
        }

        out[i] %= modulus as u128;
        values.push(F::from_canonical_u32(out[i] as u32));

	}

    // check a(x) * b(x) == out(x)
    // this check is done outside the circuit/constraints just for test purposes

    // let mut a_eval: Vec<Fq> = Vec::new();
    // let mut b_eval: Vec<Fq> = Vec::new();
    // let mut in_eval: Vec<u64> = Vec::new();
    // let mut out_eval: Vec<u64> = Vec::new();

    // let a_fq_coeffs: Vec<Fq> = a.iter().map(|&x: &u32| Fq::from(x)).collect();
    // let a_poly = DensePolynomial::from_coefficients_slice(&a_fq_coeffs);

    // let b_fq_coeffs: Vec<Fq> = b.iter().map(|&x| Fq::from(x)).collect();
    // let b_poly = DensePolynomial::from_coefficients_slice(&b_fq_coeffs);

    // let out_fq_coeffs: Vec<Fq> = out.iter().map(|&x| Fq::from(x)).collect();
    // let out_poly = DensePolynomial::from_coefficients_slice(&out_fq_coeffs);

    // for i in 0..2*N-1 {

    //     let i_as_fq = Fq::from(i as u64);

    //     // Evaluate a(x) and b(x) at x = [0..2N-1)
    //     a_eval.push(a_poly.evaluate(&i_as_fq));
    //     b_eval.push(b_poly.evaluate(&i_as_fq));

    //     // Calculate a(x) * b(x) at x = [0..2N-1)
    //     in_eval.push((a_eval[i] * b_eval[i]).into_bigint().as_ref()[0]);

    //     // Evaluate out(x) at x = [0..2N-1)
    //     out_eval.push(out_poly.evaluate(&i_as_fq).into_bigint().as_ref()[0]);

    //     println!("in_eval[{}]: {}", i, in_eval[i]);
    //     println!("out_eval[{}]: {}", i, out_eval[i]);

    // }

    // Fill in the last 3 rows with 0
    for _i in 0..3*(4*N-1) {
        values.push(F::zero());
    }
    RowMajorMatrix::new(values, 4*N-1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fmt::Debug;
    use p3_mersenne_31::Mersenne31;
    use p3_keccak::Keccak256Hash;
    use rand::{thread_rng, Rng};
    use p3_challenger::{HashChallenger, SerializingChallenger32};
    use p3_uni_stark::{prove, verify};
    use crate::gadgets::config::{initialize_config, ZkConfig, Challenger, Val};
    use crate::params::P1;

    #[test]
    fn test_poly_mul() -> Result<(), impl Debug> {

        let ZkConfig { config, byte_hash } = initialize_config();

        // generate 2 random input polynomials with N coefficients in the range of [0, N]
        let mut rng = thread_rng();
        let random_poly1: Vec<u32> = (0..N).map(|_| {
            rng.gen_range(0..P1)
        }).collect();

        let random_poly2: Vec<u32> = (0..N).map(|_| {
            rng.gen_range(0..P1)
        }).collect();

        let air = PolyMulAir { a:random_poly1.clone(), b:random_poly2.clone(), modulus:P1};

        let trace = generate_polymul_trace::<Val>(random_poly1, random_poly2, P1);

        let mut challenger: SerializingChallenger32<Mersenne31, HashChallenger<u8, Keccak256Hash, 32>> = Challenger::from_hasher(vec![], byte_hash);

        let proof = prove(&config, &air, &mut challenger, trace, &vec![]);

        let mut challenger = Challenger::from_hasher(vec![], byte_hash);
        verify(&config, &air, &mut challenger, &proof, &vec![])
    }
}
