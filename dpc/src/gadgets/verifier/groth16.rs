use crate::{crypto_primitives::nizk::groth16::Groth16, gadgets::verifier::NIZKVerifierGadget};
use algebra::{utils::ToEngineFr, AffineCurve, PairingEngine};
use snark::{
    groth16::{Proof, VerifyingKey},
    Circuit, ConstraintSystem, SynthesisError,
};
use snark_gadgets::{
    fields::FieldGadget,
    groups::GroupGadget,
    pairing::PairingGadget,
    uint8::UInt8,
    utils::{AllocGadget, EqGadget, ToBitsGadget, ToBytesGadget},
};
use std::{borrow::Borrow, marker::PhantomData};

#[derive(Derivative)]
#[derivative(Clone(bound = "P::G1Gadget: Clone, P::G2Gadget: Clone"))]
pub struct ProofGadget<
    PairingE: PairingEngine,
    ConstraintE: PairingEngine,
    P: PairingGadget<PairingE, ConstraintE>,
> {
    pub a: P::G1Gadget,
    pub b: P::G2Gadget,
    pub c: P::G1Gadget,
}

#[derive(Derivative)]
#[derivative(Clone(
    bound = "P::G1Gadget: Clone, P::GTGadget: Clone, P::G1PreparedGadget: Clone, \
             P::G2PreparedGadget: Clone, "
))]
pub struct VerifyingKeyGadget<
    PairingE: PairingEngine,
    ConstraintE: PairingEngine,
    P: PairingGadget<PairingE, ConstraintE>,
> {
    pub alpha_g1: P::G1Gadget,
    pub beta_g1:  P::G1Gadget,
    pub beta_g2:  P::G2Gadget,
    pub gamma_g2: P::G2Gadget,
    pub delta_g1: P::G1Gadget,
    pub delta_g2: P::G2Gadget,
    pub ic:       Vec<P::G1Gadget>,
}

impl<
        PairingE: PairingEngine,
        ConstraintE: PairingEngine,
        P: PairingGadget<PairingE, ConstraintE>,
    > VerifyingKeyGadget<PairingE, ConstraintE, P>
{
    pub fn prepare<CS: ConstraintSystem<ConstraintE>>(
        &self,
        mut cs: CS,
    ) -> Result<PreparedVerifyingKeyGadget<PairingE, ConstraintE, P>, SynthesisError> {
        let mut cs = cs.ns(|| "Preparing verifying key");

        let alpha_g1_pc = P::prepare_g1(&mut cs.ns(|| "Prepare alpha_g1_pc"), &self.alpha_g1)?;
        let beta_g2_pc = P::prepare_g2(&mut cs.ns(|| "Prepare beta_g2_pc"), &self.beta_g2)?;
        let gamma_g2_pc = P::prepare_g2(&mut cs.ns(|| "Prepare gamma_g2_pc"), &self.gamma_g2)?;
        let delta_g2_pc = P::prepare_g2(&mut cs.ns(|| "Prepare delta_g2_pc"), &self.delta_g2)?;

        Ok(PreparedVerifyingKeyGadget {
            alpha_g1: self.alpha_g1.clone(),
            beta_g1: self.beta_g1.clone(),
            beta_g2: self.beta_g2.clone(),
            gamma_g2: self.gamma_g2.clone(),
            delta_g1: self.delta_g1.clone(),
            delta_g2: self.delta_g2.clone(),
            alpha_g1_pc,
            beta_g2_pc,
            gamma_g2_pc,
            delta_g2_pc,
            ic: self.ic.clone(),
        })
    }
}

#[derive(Derivative)]
#[derivative(Clone(
    bound = "P::G1Gadget: Clone, P::GTGadget: Clone, P::G1PreparedGadget: Clone, \
             P::G2PreparedGadget: Clone, "
))]
pub struct PreparedVerifyingKeyGadget<
    PairingE: PairingEngine,
    ConstraintE: PairingEngine,
    P: PairingGadget<PairingE, ConstraintE>,
> {
    pub alpha_g1:    P::G1Gadget,
    pub beta_g1:     P::G1Gadget,
    pub beta_g2:     P::G2Gadget,
    pub gamma_g2:    P::G2Gadget,
    pub delta_g1:    P::G1Gadget,
    pub delta_g2:    P::G2Gadget,
    pub alpha_g1_pc: P::G1PreparedGadget,
    pub beta_g2_pc:  P::G2PreparedGadget,
    pub gamma_g2_pc: P::G2PreparedGadget, // could be negated during preparation
    pub delta_g2_pc: P::G2PreparedGadget, // could be negated during preparation
    pub ic:          Vec<P::G1Gadget>,
}

pub struct Groth16VerifierGadget<PairingE, ConstraintE, P>
where
    PairingE: PairingEngine,
    ConstraintE: PairingEngine,
    P: PairingGadget<PairingE, ConstraintE>,
{
    _pairing_engine: PhantomData<PairingE>,
    _engine:         PhantomData<ConstraintE>,
    _pairing_gadget: PhantomData<P>,
}

impl<PairingE, ConstraintE, P, C, V> NIZKVerifierGadget<Groth16<PairingE, C, V>, ConstraintE>
    for Groth16VerifierGadget<PairingE, ConstraintE, P>
where
    PairingE: PairingEngine,
    ConstraintE: PairingEngine,
    C: Circuit<PairingE>,
    V: ToEngineFr<PairingE>,
    P: PairingGadget<PairingE, ConstraintE>,
{
    type VerificationKeyGadget = VerifyingKeyGadget<PairingE, ConstraintE, P>;
    type ProofGadget = ProofGadget<PairingE, ConstraintE, P>;

    fn check_verify<'a, CS, I, T>(
        mut cs: CS,
        vk: &Self::VerificationKeyGadget,
        mut public_inputs: I,
        proof: &Self::ProofGadget,
    ) -> Result<(), SynthesisError>
    where
        CS: ConstraintSystem<ConstraintE>,
        I: Iterator<Item = &'a T>,
        T: 'a + ToBitsGadget<ConstraintE> + ?Sized,
    {
        // create PreparedVerifyingKeyGadget
        let pvk = vk.prepare(&mut cs.ns(|| "Prepare vk"))?;

        // A * B + inputs * (-gamma) + C * (-delta) = alpha * beta
        // where input  = \sum_{i=0}^l input_i pvk.ic[i]

        // TODO: rename to inputs
        let g_psi = {
            let mut cs = cs.ns(|| "Process input");
            let mut g_psi = pvk.ic[0].clone();
            let mut input_len = 1;
            for (i, (input, b)) in public_inputs
                .by_ref()
                .zip(pvk.ic.iter().skip(1))
                .enumerate()
            {
                let input_bits = input.to_bits(cs.ns(|| format!("Input {}", i)))?;
                g_psi = b.mul_bits(cs.ns(|| format!("Mul {}", i)), &g_psi, input_bits.iter())?;
                input_len += 1;
            }
            // Check that the input and the query in the verification are of the
            // same length.
            assert!(input_len == pvk.ic.len() && public_inputs.next().is_none());
            g_psi
        };

        let test1_exp = {
            let a_prep = P::prepare_g1(cs.ns(|| "A prep"), &proof.a)?;
            let b_prep = P::prepare_g2(cs.ns(|| "B prep"), &proof.b)?;
            let neg_c = proof.c.clone().negate(cs.ns(|| "neg C"))?;
            let neg_c_prep = P::prepare_g1(cs.ns(|| "C prep"), &neg_c)?;

            let neg_psi = g_psi.clone().negate(cs.ns(|| "neg inputs"))?;
            let neg_psi_prep = P::prepare_g1(cs.ns(|| "inputs prep"), &neg_psi)?;

            let neg_alpha = pvk.alpha_g1.clone().negate(cs.ns(|| "neg alpha"))?;
            let neg_alpha_prep = P::prepare_g1(cs.ns(|| "neg alpah prep"), &neg_alpha)?;

            P::miller_loop(
                cs.ns(|| "Miller loop"),
                &[
                    a_prep,         // done
                    neg_psi_prep,   // done
                    neg_c_prep,     // done
                    neg_alpha_prep, // done
                ],
                &[
                    b_prep,                  // done
                    pvk.gamma_g2_pc.clone(), // done
                    pvk.delta_g2_pc.clone(), // done
                    pvk.beta_g2_pc.clone(),  // done
                ],
            )?
        };

        let test = P::final_exponentiation(cs.ns(|| "Final Exp 1"), &test1_exp).unwrap();

        let one = P::GTGadget::one(cs.ns(|| "GT One"))?;
        test.enforce_equal(cs.ns(|| "Test 1"), &one)?;

        Ok(())
    }
}

impl<PairingE, ConstraintE, P> AllocGadget<VerifyingKey<PairingE>, ConstraintE>
    for VerifyingKeyGadget<PairingE, ConstraintE, P>
where
    PairingE: PairingEngine,
    ConstraintE: PairingEngine,
    P: PairingGadget<PairingE, ConstraintE>,
{
    #[inline]
    fn alloc<FN, T, CS: ConstraintSystem<ConstraintE>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<VerifyingKey<PairingE>>,
    {
        value_gen().and_then(|vk| {
            let VerifyingKey {
                alpha_g1,
                beta_g1,
                beta_g2,
                gamma_g2,
                delta_g1,
                delta_g2,
                ic,
            } = vk.borrow().clone();

            let alpha_g1 =
                P::G1Gadget::alloc(cs.ns(|| "alpha_g1"), || Ok(alpha_g1.into_projective()))?;
            let beta_g1 =
                P::G1Gadget::alloc(cs.ns(|| "beta_g1"), || Ok(beta_g1.into_projective()))?;
            let beta_g2 =
                P::G2Gadget::alloc(cs.ns(|| "beta_g2"), || Ok(beta_g2.into_projective()))?;
            let gamma_g2 =
                P::G2Gadget::alloc(cs.ns(|| "gamma_g2"), || Ok(gamma_g2.into_projective()))?;
            let delta_g1 =
                P::G1Gadget::alloc(cs.ns(|| "delta_g1"), || Ok(delta_g1.into_projective()))?;
            let delta_g2 =
                P::G2Gadget::alloc(cs.ns(|| "delta_g2"), || Ok(delta_g2.into_projective()))?;

            let ic = ic
                .into_iter()
                .enumerate()
                .map(|(i, query_i)| {
                    P::G1Gadget::alloc(cs.ns(|| format!("query_{}", i)), || {
                        Ok(query_i.into_projective())
                    })
                })
                .collect::<Vec<_>>()
                .into_iter()
                .collect::<Result<_, _>>()?;

            Ok(Self {
                alpha_g1,
                beta_g1,
                beta_g2,
                gamma_g2,
                delta_g1,
                delta_g2,
                ic,
            })
        })
    }
    #[inline]
    fn alloc_input<FN, T, CS: ConstraintSystem<ConstraintE>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<VerifyingKey<PairingE>>,
    {
        value_gen().and_then(|vk| {
            let VerifyingKey {
                alpha_g1,
                beta_g1,
                beta_g2,
                gamma_g2,
                delta_g1,
                delta_g2,
                ic,
            } = vk.borrow().clone();

            let alpha_g1 =
                P::G1Gadget::alloc_input(cs.ns(|| "alpha_g1"), || Ok(alpha_g1.into_projective()))?;
            let beta_g1 =
                P::G1Gadget::alloc_input(cs.ns(|| "beta_g1"), || Ok(beta_g1.into_projective()))?;
            let beta_g2 =
                P::G2Gadget::alloc_input(cs.ns(|| "beta_g2"), || Ok(beta_g2.into_projective()))?;
            let gamma_g2 =
                P::G2Gadget::alloc_input(cs.ns(|| "gamma_g2"), || Ok(gamma_g2.into_projective()))?;
            let delta_g1 =
                P::G1Gadget::alloc_input(cs.ns(|| "delta_g1"), || Ok(delta_g1.into_projective()))?;
            let delta_g2 =
                P::G2Gadget::alloc_input(cs.ns(|| "delta_g2"), || Ok(delta_g2.into_projective()))?;

            let ic = ic
                .into_iter()
                .enumerate()
                .map(|(i, query_i)| {
                    P::G1Gadget::alloc_input(cs.ns(|| format!("query_{}", i)), || {
                        Ok(query_i.into_projective())
                    })
                })
                .collect::<Vec<_>>()
                .into_iter()
                .collect::<Result<_, _>>()?;

            Ok(Self {
                alpha_g1,
                beta_g1,
                beta_g2,
                gamma_g2,
                delta_g1,
                delta_g2,
                ic,
            })
        })
    }
}

impl<PairingE, ConstraintE, P> AllocGadget<Proof<PairingE>, ConstraintE>
    for ProofGadget<PairingE, ConstraintE, P>
where
    PairingE: PairingEngine,
    ConstraintE: PairingEngine,
    P: PairingGadget<PairingE, ConstraintE>,
{
    #[inline]
    fn alloc<FN, T, CS: ConstraintSystem<ConstraintE>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<Proof<PairingE>>,
    {
        value_gen().and_then(|proof| {
            let Proof { a, b, c } = proof.borrow().clone();
            let a = P::G1Gadget::alloc_checked(cs.ns(|| "a"), || Ok(a.into_projective()))?;
            let b = P::G2Gadget::alloc_checked(cs.ns(|| "b"), || Ok(b.into_projective()))?;
            let c = P::G1Gadget::alloc_checked(cs.ns(|| "c"), || Ok(c.into_projective()))?;
            Ok(Self { a, b, c })
        })
    }

    #[inline]
    fn alloc_input<FN, T, CS: ConstraintSystem<ConstraintE>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<Proof<PairingE>>,
    {
        value_gen().and_then(|proof| {
            let Proof { a, b, c } = proof.borrow().clone();
            // We don't need to check here because the prime order check can be performed
            // in plain.
            let a = P::G1Gadget::alloc_input(cs.ns(|| "a"), || Ok(a.into_projective()))?;
            let b = P::G2Gadget::alloc_input(cs.ns(|| "b"), || Ok(b.into_projective()))?;
            let c = P::G1Gadget::alloc_input(cs.ns(|| "c"), || Ok(c.into_projective()))?;
            Ok(Self { a, b, c })
        })
    }
}

impl<PairingE, ConstraintE, P> ToBytesGadget<ConstraintE>
    for VerifyingKeyGadget<PairingE, ConstraintE, P>
where
    PairingE: PairingEngine,
    ConstraintE: PairingEngine,
    P: PairingGadget<PairingE, ConstraintE>,
{
    #[inline]
    fn to_bytes<CS: ConstraintSystem<ConstraintE>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&self.alpha_g1.to_bytes(&mut cs.ns(|| "h_g2 to bytes"))?);
        bytes.extend_from_slice(
            &self
                .beta_g1
                .to_bytes(&mut cs.ns(|| "g_alpha_g1 to bytes"))?,
        );
        bytes.extend_from_slice(&self.beta_g2.to_bytes(&mut cs.ns(|| "h_beta_g2 to bytes"))?);
        bytes.extend_from_slice(
            &self
                .gamma_g2
                .to_bytes(&mut cs.ns(|| "g_gamma_g2 to bytes"))?,
        );
        bytes.extend_from_slice(&self.delta_g1.to_bytes(&mut cs.ns(|| "deltag1 to bytes"))?);
        bytes.extend_from_slice(&self.delta_g2.to_bytes(&mut cs.ns(|| "deltag2 to bytes"))?);
        for (i, q) in self.ic.iter().enumerate() {
            let mut cs = cs.ns(|| format!("Iteration {}", i));
            bytes.extend_from_slice(&q.to_bytes(&mut cs.ns(|| "q"))?);
        }
        Ok(bytes)
    }

    fn to_bytes_strict<CS: ConstraintSystem<ConstraintE>>(
        &self,
        cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        self.to_bytes(cs)
    }
}

#[cfg(test)]
mod test {
    use snark::{groth16::*, Circuit, ConstraintSystem, SynthesisError};

    use super::*;
    use algebra::{
        curves::{bls12_377::Bls12_377, sw6::SW6},
        fields::bls12_377::Fr,
        BitIterator, PrimeField,
    };
    use rand::{thread_rng, Rng};
    use snark_gadgets::{
        boolean::Boolean, pairing::bls12_377::PairingGadget as Bls12_377PairingGadget,
        test_constraint_system::TestConstraintSystem, utils::AllocGadget,
    };

    type TestProofSystem = Groth16<Bls12_377, Bench<Bls12_377>, Fr>;
    type TestVerifierGadget = Groth16VerifierGadget<Bls12_377, SW6, Bls12_377PairingGadget>;
    type TestProofGadget = ProofGadget<Bls12_377, SW6, Bls12_377PairingGadget>;
    type TestVkGadget = VerifyingKeyGadget<Bls12_377, SW6, Bls12_377PairingGadget>;

    struct Bench<E: PairingEngine> {
        inputs:          Vec<Option<E::Fr>>,
        num_constraints: usize,
    }

    impl<E: PairingEngine> Circuit<E> for Bench<E> {
        fn synthesize<CS: ConstraintSystem<E>>(self, cs: &mut CS) -> Result<(), SynthesisError> {
            assert!(self.inputs.len() >= 2);
            assert!(self.num_constraints >= self.inputs.len());

            let mut variables: Vec<_> = Vec::with_capacity(self.inputs.len());
            for (i, input) in self.inputs.into_iter().enumerate() {
                let input_var = cs.alloc_input(
                    || format!("Input {}", i),
                    || input.ok_or(SynthesisError::AssignmentMissing),
                )?;
                variables.push((input, input_var));
            }

            for i in 0..self.num_constraints {
                let new_entry = {
                    let (input_1_val, input_1_var) = variables[i];
                    let (input_2_val, input_2_var) = variables[i + 1];
                    let result_val = input_1_val
                        .and_then(|input_1| input_2_val.map(|input_2| input_1 * &input_2));
                    let result_var = cs.alloc(
                        || format!("Result {}", i),
                        || result_val.ok_or(SynthesisError::AssignmentMissing),
                    )?;
                    cs.enforce(
                        || format!("Enforce constraint {}", i),
                        |lc| lc + input_1_var,
                        |lc| lc + input_2_var,
                        |lc| lc + result_var,
                    );
                    (result_val, result_var)
                };
                variables.push(new_entry);
            }
            Ok(())
        }
    }

    #[test]
    fn groth16_verifier_test() {
        let num_inputs = 4;
        let num_constraints = num_inputs;
        let rng = &mut thread_rng();
        let mut inputs: Vec<Option<Fr>> = Vec::with_capacity(num_inputs);
        for _ in 0..num_inputs {
            inputs.push(Some(rng.gen()));
        }
        let params = {
            let c = Bench::<Bls12_377> {
                inputs: vec![None; num_inputs],
                num_constraints,
            };

            generate_random_parameters(c, rng).unwrap()
        };

        {
            let proof = {
                // Create an instance of our circuit (with the
                // witness)
                let c = Bench {
                    inputs: inputs.clone(),
                    num_constraints,
                };
                // Create a groth16 proof with our parameters.
                create_random_proof(c, &params, rng).unwrap()
            };

            // assert!(!verify_proof(&pvk, &proof, &[a]).unwrap());
            let mut cs = TestConstraintSystem::<SW6>::new();

            let inputs: Vec<_> = inputs.into_iter().map(|input| input.unwrap()).collect();
            let mut input_gadgets = Vec::new();

            {
                let mut cs = cs.ns(|| "Allocate Input");
                for (i, input) in inputs.into_iter().enumerate() {
                    let mut input_bits = BitIterator::new(input.into_repr()).collect::<Vec<_>>();
                    // Input must be in little-endian, but BitIterator outputs in big-endian.
                    input_bits.reverse();

                    println!("Input lenghts {:?}", input_bits.len());
                    let input_bits =
                        Vec::<Boolean>::alloc_input(cs.ns(|| format!("Input {}", i)), || {
                            Ok(input_bits)
                        })
                        .unwrap();
                    input_gadgets.push(input_bits);
                }
            }
            println!("Len of input gadgets {:?}", input_gadgets.len());
            println!("Constraints before Vk {:?}", cs.num_constraints());

            let vk_gadget = TestVkGadget::alloc_input(cs.ns(|| "Vk"), || Ok(&params.vk)).unwrap();
            println!("Constraints after Vk Gadget {:?}", cs.num_constraints());

            //            checks if vk elements on curve  ... why is that necessary ..
            // public input?
            let proof_gadget =
                TestProofGadget::alloc(cs.ns(|| "Proof"), || Ok(proof.clone())).unwrap();
            //           checks if on curve and if in correct subgroup
            <TestVerifierGadget as NIZKVerifierGadget<TestProofSystem, SW6>>::check_verify(
                cs.ns(|| "Verify"),
                &vk_gadget,
                input_gadgets.iter(),
                &proof_gadget,
            )
            .unwrap();

            println!(
                "Constraints after Verfier Gadget {:?}, for {:?} inputs ",
                cs.num_constraints(),
                num_inputs
            );

            if !cs.is_satisfied() {
                println!("=========================================================");
                println!("Unsatisfied constraints:");
                println!("{:?}", cs.which_is_unsatisfied().unwrap());
                println!("====================================================c=====");
            }

            // cs.print_named_objects();
            assert!(cs.is_satisfied());
        }
    }
}
