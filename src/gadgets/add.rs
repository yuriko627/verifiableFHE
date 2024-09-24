use p3_air::{Air, AirBuilder, BaseAir};
use p3_field::{AbstractField, Field};
use p3_matrix::Matrix;
use p3_matrix::dense::RowMajorMatrix;
use crate:: params::N;
use std::ops::{Add, Sub};

// Define AIR constraint inputs
pub struct PolyAddAir {
	pub a: Vec<u32>,
	pub b: Vec<u32>,
	pub modulus: u32
}

/*
Polynomial Addition Air
Input:
- a = a[0] + a[1] * X + ... + a[N-1] * X^{N-1}
- b = b[0] + b[1] * X + ... + b[N-1] * X^{N-1}
- mod: FHE ciphertext modulus
Output:
- out = out[0] + out[1] * X + ... + out[2N-2] * X^{2N-2}

Note:
- PolyAddAir does not have a state transition. Values required for constraints are all stored in one row.
- While output polynomial `out` is calculated manually by generate_polyadd_trace(), we prove that this addition was done correctly, by enforcing a constraint such that a(x)+b(x) === out(x)  at x = [0..2N-1) based on Lagrange polynomial interpolation.
*/
impl<F: Field> BaseAir<F> for PolyAddAir {
    // Air Table looks like this
    // row:[      a: N      ][      b: N      ][mod:1][      out(x): N      ]
    //     ^------------------inputs-----------------^^-calculated by generate_polyadd_trace
    //     [0..............................................................0]
    //     [0..............................................................0]
    //     [0..............................................................0]
    fn width(&self) -> usize {
        3*N+1
    }
}

// Define constraints
impl<AB: AirBuilder> Air<AB> for PolyAddAir {
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let row = main.row_slice(0);

        // Enforce self.a and self.b as 2 input polynomials
		for i in 0..N {
			builder.when_first_row().assert_eq(row[i], AB::Expr::from_canonical_u32(self.a[i]));
			builder.when_first_row().assert_eq(row[i+N], AB::Expr::from_canonical_u32(self.b[i]));
		}

        // Enforce self.modulus as mod
        builder.when_first_row().assert_eq(row[2*N], AB::Expr::from_canonical_u32(self.modulus));

        /*
        We want to ensure a[i] + b[i]) === out[i] mod p
        where p = non-native 31-bits modulus for FHE, which is smaller than n = native modulus for ZK (Mersenne31).
        -> We can enforce a[i] + b[i] === q[i] * p[i] + out[i] for N coefficients ...(1)
        However, a[i] + b[i] is at most p-1 + p-1 = 2*p-2, which overflows n.
        So, we will *virtually* expand the field size to 2^t*n > 2*p,
        and break it down into two constraints by CRT:
        1) a[i] + b[i] ===  q * p + out[i] (mod 2^t)
        2) a[i] + b[i] ===  q * p + out[i] (mod n)
        Note:
        - we can apply CRT because 2^t and n are co-prime to each other.
        - qotient q_1 for 1) and q_2 for 2) are both pre-computed outside the circuit.

        Toy example:
        Suppose p = 5, n = 7, a = 3, b = 4, out = 2
        3 + 4 = 7 = 2 === 2 mod 5 (inside FHE)
        3 + 4 = 7 = 0 =/= 2 mod 7 (inside ZKP)
        -> LHS evaluates to inconsistent values

        5*2 = 10 < 2^1 * 7 = 14 -> let's expand the field to mod 14
        Now, break the constraint down by CRT.
        1) a + b ===  q_1 * p + out (mod 2) where q_1 = 1 is precomputed
        2) a + b ===  q_2 * p + out (mod 7) where q_2 = 1 is precomputed
        1) (3 + 4) % 2 evaluates to 1 === (1 * 5 + 2) % 2 evaluates to 1 (mod 2)
        2) (3 + 4) % 7 evaluates to 0 === (1 * 5 + 2) % 7 evaluates to 0 (mod 7)
        */

    }
}

// Define a function to generate execution trace
pub fn generate_polyadd_trace<F: Field>(a:Vec<u32>, b:Vec<u32>, modulus: u32) -> RowMajorMatrix<F> {
    let mut values: Vec<F>= Vec::with_capacity(4*(3*N+1)); // 4 is the minimum number of rows required

	// Add input polynomials to values vector
	for i in 0..N {
        println!("a[{}]: {}", i, a[i]);
		values.push(F::from_canonical_u32(a[i]));
	}
	for i in 0..N {
        println!("b[{}]: {}", i, b[i]);
		values.push(F::from_canonical_u32(b[i]));
	}
    // Add modulus to values vector
    values.push(F::from_canonical_u32(modulus));

	// Add the 2 polynomials and push it to values vector
	for i in 0..N {
		values.push(F::from_canonical_u32((a[i] + b[i]) % modulus));
        println!("out[{}]: {}", i, (a[i] + b[i]) % modulus);
	}

	// Fill in the rest of the slots (last 3 rows) with 0
	for _ in 0..3*(3*N+1) {
		values.push(F::zero());
	}
    RowMajorMatrix::new(values, 3*N+1)

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
    fn test_poly_add() -> Result<(), impl Debug> {

        let ZkConfig { config, byte_hash } = initialize_config();

        // generate 2 random input polynomials with N coefficients in the range of [0, N]
        let mut rng = thread_rng();
        let random_poly1: Vec<u32> = (0..N).map(|_| {
            rng.gen_range(0..P1)
        }).collect();

        let random_poly2: Vec<u32> = (0..N).map(|_| {
            rng.gen_range(0..P1)
        }).collect();

        let air = PolyAddAir { a:random_poly1.clone(), b:random_poly2.clone(), modulus:P1 };

        let trace = generate_polyadd_trace::<Val>(random_poly1, random_poly2, P1);

        let mut challenger: SerializingChallenger32<Mersenne31, HashChallenger<u8, Keccak256Hash, 32>> = Challenger::from_hasher(vec![], byte_hash);
        let proof = prove(&config, &air, &mut challenger, trace, &vec![]);

        let mut challenger = Challenger::from_hasher(vec![], byte_hash);
        verify(&config, &air, &mut challenger, &proof, &vec![])
    }
}

