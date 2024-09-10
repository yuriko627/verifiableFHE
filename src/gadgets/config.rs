use std::marker::PhantomData;
use p3_mersenne_31::Mersenne31;
use p3_challenger::{HashChallenger, SerializingChallenger32};
use p3_circle::CirclePcs;
use p3_commit::ExtensionMmcs;
use p3_field::extension::BinomialExtensionField;
use p3_fri::FriConfig;
use p3_keccak::Keccak256Hash;
use p3_merkle_tree::FieldMerkleTreeMmcs;
use p3_symmetric::{CompressionFunctionFromHasher, SerializingHasher32};
use p3_uni_stark::StarkConfig;
use tracing_forest::util::LevelFilter;
use tracing_forest::ForestLayer;
use tracing_subscriber::layer::SubscriberExt;
use tracing_subscriber::util::SubscriberInitExt;
use tracing_subscriber::{EnvFilter, Registry};

// Define a struct to hold all configuration types
pub struct ZkConfig {
    pub config: StarkConfig<Pcs, BinomialExtensionField<Mersenne31, 3>, SerializingChallenger32<Mersenne31, HashChallenger<u8, Keccak256Hash, 32>>>,
    pub byte_hash: Keccak256Hash,
}

// Type aliases for the ZK system configuration
pub type Val = Mersenne31;
pub type Challenge = BinomialExtensionField<Val, 3>;
pub type ByteHash = Keccak256Hash;
pub type FieldHash = SerializingHasher32<ByteHash>;
pub type MyCompress = CompressionFunctionFromHasher<u8, ByteHash, 2, 32>;
pub type ValMmcs = FieldMerkleTreeMmcs<Val, u8, FieldHash, MyCompress, 32>;
pub type ChallengeMmcs = ExtensionMmcs<Val, Challenge, ValMmcs>;
pub type Challenger = SerializingChallenger32<Val, HashChallenger<u8, ByteHash, 32>>;
pub type Pcs = CirclePcs<Val, ValMmcs, ChallengeMmcs>;

pub fn initialize_config() -> ZkConfig {

    let env_filter = EnvFilter::builder()
        .with_default_directive(LevelFilter::INFO.into())
        .from_env_lossy();

    Registry::default()
        .with(env_filter)
        .with(ForestLayer::default())
        .try_init() // Use try_init() to prevent conflicts
        .ok(); // Ignore errors if already initialized

    // Initialize zk system configuration
    let byte_hash = ByteHash {};
    let field_hash = FieldHash::new(Keccak256Hash {});
    let compress = MyCompress::new(byte_hash);

    let val_mmcs = ValMmcs::new(field_hash, compress);
    let challenge_mmcs = ChallengeMmcs::new(val_mmcs.clone());

    let fri_config = FriConfig {
        log_blowup: 1,
        num_queries: 100,
        proof_of_work_bits: 16,
        mmcs: challenge_mmcs,
    };

    let pcs = Pcs {
        mmcs: val_mmcs,
        fri_config,
        _phantom: PhantomData,
    };

    let config = StarkConfig::new(pcs);

    ZkConfig {
        config,
        byte_hash,
    }
}
