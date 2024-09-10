use p3_air::{Air, AirBuilder, BaseAir};
use p3_field::{AbstractField, Field};
use p3_matrix::Matrix;
use p3_matrix::dense::RowMajorMatrix;


// Define AIR constraint inputs
pub struct PolyAddAir {
	pub poly1: Vec<u32>,
	pub poly2: Vec<u32>,
	pub added_poly: Vec<u32>
}

impl<F: Field> BaseAir<F> for PolyAddAir {
    fn width(&self) -> usize {
        10500 // poly1, pol2, added_poly's coeffs (3500 for each) in one row
    }
}

// Define your constraints
impl<AB: AirBuilder> Air<AB> for PolyAddAir {
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let current = main.row_slice(0);

        // Enforce input polynomial values and the added polynomial values
		for i in 0..3500 {
			builder.when_first_row().assert_eq(current[i], AB::Expr::from_canonical_u32(self.poly1[i]));
			builder.when_first_row().assert_eq(current[i+3500], AB::Expr::from_canonical_u32(self.poly2[i]));
            builder.when_first_row().assert_eq(current[i+7000], AB::Expr::from_canonical_u32(self.added_poly[i]));
		}

    }
}

// Define a function to generate execution trace
pub fn generate_polyadd_trace<F: Field>(poly1:Vec<u32>, poly2:Vec<u32>) -> RowMajorMatrix<F> {
    let mut values: Vec<F>= Vec::with_capacity(4 * 10500); // 4 is the minimum number of rows required

	// fill in the states in each iteration in the `values` vector
	for i in 0..3500 {
		values.push(F::from_canonical_u32(poly1[i]));
	}
	for i in 0..3500 {
		values.push(F::from_canonical_u32(poly2[i]));
	}

	// Add the 2 polynomials and push it to values vector
	for i in 0..3500 {
		values.push(F::from_canonical_u32((poly1[i] + poly2[i]) % 536870939));
	}

	// Fill in the rest of the slots (last 3 rows) with 0
	for _ in 0..(3 * 10500) {
		values.push(F::zero());
	}
    RowMajorMatrix::new(values, 10500)

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

    #[test]
    fn test_poly_add() -> Result<(), impl Debug> {

        let ZkConfig { config, byte_hash } = initialize_config();

        // generate 2 random polynomial with 3500 terms each
        let mut rng = thread_rng();
        let random_poly1: Vec<u32> = (0..3500).map(|_| {
            rng.gen_range(0..536870939) // chose an arbitrary 30-bits prime number
        }).collect();

        let random_poly2: Vec<u32> = (0..3500).map(|_| {
            rng.gen_range(0..536870939)
        }).collect();

        let mut added_poly: Vec<u32> = Vec::with_capacity(3500);

        // Add the 2 polynomials
        for i in 0..3500 {
            added_poly.push((random_poly1[i] + random_poly2[i]) % 536870939);
        }

        let air = PolyAddAir { poly1: random_poly1.clone(), poly2: random_poly2.clone(), added_poly };

        let trace = generate_polyadd_trace::<Val>(random_poly1, random_poly2);

        let mut challenger: SerializingChallenger32<Mersenne31, HashChallenger<u8, Keccak256Hash, 32>> = Challenger::from_hasher(vec![], byte_hash);
        let proof = prove(&config, &air, &mut challenger, trace, &vec![]);

        let mut challenger = Challenger::from_hasher(vec![], byte_hash);
        verify(&config, &air, &mut challenger, &proof, &vec![])
    }
}

