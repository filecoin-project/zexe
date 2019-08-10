use algebra::bytes::ToBytes;
use rand::Rng;
use std::hash::Hash;

pub mod injective_map;
pub mod pedersen;

use crate::Error;

pub trait FixedLengthCRH {
    const INPUT_SIZE_BITS: usize;
    type Output: ToBytes + Clone + Eq + Hash + Default;
    type Parameters: Clone + Default;

    fn setup<R: Rng>(r: &mut R) -> Result<Self::Parameters, Error>;
    fn evaluate(parameters: &Self::Parameters, input: &[u8]) -> Result<Self::Output, Error>;
}

#[cfg(test)]
mod test {
    use super::{
        pedersen::{PedersenCRH, PedersenWindow},
        FixedLengthCRH,
    };
    use crate::Error;
    use algebra::{
        bytes::ToBytes,
        curves::jubjub::JubJubAffine as JubJub,
        fields::{jubjub::fr::Fr, Field},
        to_bytes,
    };
    use rand::thread_rng;

    fn crh_test<C: FixedLengthCRH>(input: &[u8]) -> Result<C::Output, Error> {
        let rng = &mut thread_rng();
        let parameters = C::setup::<_>(rng).unwrap();
        C::evaluate(&parameters, input)
    }

    #[ignore]
    #[test]
    fn pedersen_crh_test() {
        #[derive(Clone)]
        pub(super) struct Window;

        impl PedersenWindow for Window {
            const WINDOW_SIZE: usize = 4;
            const NUM_WINDOWS: usize = 64;
        }

        let input = to_bytes!(Fr::zero()).unwrap();

        println!("{}", input.len());

        crh_test::<PedersenCRH<JubJub, Window>>(&input).unwrap();
    }

    #[test]
    fn test_pedersen_crh_fixed() {
        use algebra::{
            biginteger::BigInteger,
            curves::{jubjub::JubJubProjective as JubJub, ProjectiveCurve},
            fields::PrimeField,
        };

        let input = b"some bytes";

        #[derive(Clone)]
        pub(super) struct Window;

        impl PedersenWindow for Window {
            const WINDOW_SIZE: usize = 4;
            const NUM_WINDOWS: usize = 64;
        }

        let hash = crh_test::<PedersenCRH<JubJub, Window>>(&input[..]).unwrap();
        let mut hash_bytes = Vec::new();
        hash.into_affine()
            .x
            .into_repr()
            .write_le(&mut hash_bytes)
            .expect("failed to write result hash");
        let expected = vec![
            237, 70, 41, 231, 39, 180, 131, 120, 36, 36, 119, 199, 200, 225, 153, 242, 106, 116,
            70, 9, 12, 249, 169, 84, 105, 38, 225, 115, 165, 188, 98, 25,
        ];

        assert_eq!(hash_bytes, expected);
    }

    #[test]
    fn pedersen_crh_benchmark() {
        use std::time::Instant;

        let rng = &mut thread_rng();

        const SAMPLES: u64 = 2;

        println!("Running {} times and averaging the time spent.", SAMPLES);

        #[derive(Clone)]
        pub(super) struct Window;

        impl PedersenWindow for Window {
            const WINDOW_SIZE: usize = 4;
            const NUM_WINDOWS: usize = 64;
        }

        let input = to_bytes!(Fr::one()).unwrap();

        let setup_start = Instant::now();
        for _ in 0..SAMPLES {
            let _parameters = PedersenCRH::<JubJub, Window>::setup::<_>(rng).unwrap();
        }
        let mut elapsed = setup_start.elapsed();
        let mut in_ms = elapsed.as_secs() * 1000 + elapsed.subsec_nanos() as u64 / 1_000_000;
        println!("Setup time: {:?}ms", in_ms / SAMPLES);

        let parameters = PedersenCRH::<JubJub, Window>::setup::<_>(rng).unwrap();

        let evaluate_start = Instant::now();
        for _ in 0..SAMPLES {
            let _hash =
                PedersenCRH::<JubJub, Window>::evaluate(&parameters, input.as_slice()).unwrap();
        }
        elapsed = evaluate_start.elapsed();
        in_ms = elapsed.as_secs() * 1000 + elapsed.subsec_nanos() as u64 / 1_000_000;
        println!("Average evaluation time: {:?}ms", in_ms / SAMPLES);
    }
}
