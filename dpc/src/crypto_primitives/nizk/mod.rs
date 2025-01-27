use algebra::bytes::ToBytes;
use rand::Rng;

pub mod gm17;
pub use self::gm17::Gm17;
pub mod groth16;
pub use self::groth16::Groth16;

use crate::Error;

pub trait NIZK {
    type Circuit;
    type AssignedCircuit;
    type VerifierInput: ?Sized;
    type ProvingParameters: Clone;
    type VerificationParameters: Clone + Default + From<Self::PreparedVerificationParameters>;
    type PreparedVerificationParameters: Clone + Default + From<Self::VerificationParameters>;
    type Proof: ToBytes + Clone + Default;

    fn setup<R: Rng>(
        circuit: Self::Circuit,
        rng: &mut R,
    ) -> Result<
        (
            Self::ProvingParameters,
            Self::PreparedVerificationParameters,
        ),
        Error,
    >;

    fn prove<R: Rng>(
        parameter: &Self::ProvingParameters,
        input_and_witness: Self::AssignedCircuit,
        rng: &mut R,
    ) -> Result<Self::Proof, Error>;

    fn verify(
        verifier_key: &Self::PreparedVerificationParameters,
        input: &Self::VerifierInput,
        proof: &Self::Proof,
    ) -> Result<bool, Error>;
}

#[cfg(test)]
mod test {
    use rand::thread_rng;
    use std::ops::AddAssign;

    #[test]
    fn test_gm17() {
        use crate::crypto_primitives::nizk::{gm17::Gm17, NIZK};
        use algebra::{curves::bls12_381::Bls12_381, fields::bls12_381::Fr, Field};
        use snark::{Circuit, ConstraintSystem, SynthesisError};

        #[derive(Copy, Clone)]
        struct R1CSCircuit {
            x:   Option<Fr>,
            sum: Option<Fr>,
            w:   Option<Fr>,
        }

        impl R1CSCircuit {
            pub(super) fn new(x: Fr, sum: Fr, w: Fr) -> Self {
                Self {
                    x:   Some(x),
                    sum: Some(sum),
                    w:   Some(w),
                }
            }
        }

        impl Circuit<Bls12_381> for R1CSCircuit {
            fn synthesize<CS: ConstraintSystem<Bls12_381>>(
                self,
                cs: &mut CS,
            ) -> Result<(), SynthesisError> {
                let input = cs.alloc_input(|| "x", || Ok(self.x.unwrap()))?;
                let sum = cs.alloc_input(|| "sum", || Ok(self.sum.unwrap()))?;
                let witness = cs.alloc(|| "w", || Ok(self.w.unwrap()))?;

                cs.enforce(
                    || "check_one",
                    |lc| lc + sum,
                    |lc| lc + CS::one(),
                    |lc| lc + input + witness,
                );
                Ok(())
            }
        }

        let mut sum = Fr::one();
        sum.add_assign(&Fr::one());
        let circuit = R1CSCircuit::new(Fr::one(), sum, Fr::one());

        let rng = &mut thread_rng();

        let parameters = Gm17::<Bls12_381, R1CSCircuit, [Fr]>::setup(circuit, rng).unwrap();

        let proof =
            Gm17::<Bls12_381, R1CSCircuit, [Fr]>::prove(&parameters.0, circuit, rng).unwrap();

        let result =
            Gm17::<Bls12_381, R1CSCircuit, [Fr]>::verify(&parameters.1, &[Fr::one(), sum], &proof)
                .unwrap();
        assert!(result);
    }
}
