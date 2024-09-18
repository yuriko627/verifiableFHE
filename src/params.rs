// N: number of ciphertext polynomial coefficients/terms
pub const N: usize = 3500;

// 3 ciphertext modulus in the RNS-decomposed fields
pub const P1: u32 = 1085276161; // 31-bits, generator: 11
pub const P2: u32 = 1092616193; // 31-bits, generator: 3
pub const P3: u32 = 1095761921; // 31-bits, generator: 3

// P: ciphertext modulus in the original ring
// 1299343865123888653488095233: 91-bits
pub const P: u128 = P1 as u128 * P2 as u128 * P3 as u128;
