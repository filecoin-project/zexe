use algebra::{
    biginteger::BigInteger,
    fields::{Field, FpParameters, PrimeField},
    groups::Group,
    BitIterator8,
};
use blake2s_simd::Params;
use rand::Rng;
use std::{
    fmt::{Debug, Formatter, Result as FmtResult},
    marker::PhantomData,
};

use super::FixedLengthCRH;
use crate::Error;

/// First 64 bytes of the BLAKE2s input during group hash.
/// This is chosen to be some random string that we couldn't have anticipated
/// when we designed the algorithm, for rigidity purposes.
/// We deliberately use an ASCII hex string of 32 bytes here.
pub const GH_FIRST_BLOCK: &'static [u8; 64] =
    b"096b36a5804bfacef1691e173c366a47ff5ba84a44f26ddd7e8d9f79d5b42df0";

/// BLAKE2s Personalization for Pedersen hash generators.
pub const PEDERSEN_HASH_GENERATORS_PERSONALIZATION: &'static [u8; 8] = b"Zcash_PH";

pub trait PedersenWindow: Clone {
    const WINDOW_SIZE: usize;
    const NUM_WINDOWS: usize;
}

#[derive(Clone, Default)]
pub struct PedersenParameters<G: Group> {
    pub generators: Vec<G>,
    pub exp_table:  Vec<Vec<Vec<G>>>,
}

impl<G: Group> PedersenParameters<G> {
    pub fn chunks_per_generator(&self) -> usize {
        63
    }
}

pub struct PedersenCRH<G: Group, W: PedersenWindow> {
    group:  PhantomData<G>,
    window: PhantomData<W>,
}

fn find_group_hash<G: Group>(m: &[u8], personalization: &[u8; 8]) -> G {
    let mut tag = m.to_vec();
    let i = tag.len();
    tag.push(0u8);

    loop {
        let gh = group_hash(&tag, personalization);
        println!("group hash: {:?}", &gh);
        // We don't want to overflow and start reusing generators
        assert!(tag[i] != u8::max_value());
        tag[i] += 1;

        if let Some(gh) = gh {
            break gh;
        }
    }
}

fn group_hash<G: Group>(tag: &[u8], personalization: &[u8]) -> Option<G> {
    assert_eq!(personalization.len(), 8);

    // TODO: do we need to adjust for the different modulus bits in here? This check
    // fails for the different curves.
    // Check to see that scalar field is 255 bits
    // assert!(<G::ScalarField as PrimeField>::Params::MODULUS_BITS == 255);

    let mut p = Params::new();
    p.hash_length(32);
    p.personal(personalization);
    p.key(&[]);
    p.salt(&[]);
    let mut h = p.to_state();
    h.update(GH_FIRST_BLOCK);
    h.update(tag);
    let h = h.finalize();
    assert!(h.as_ref().len() == 32);

    println!("hash: {:?}", h.as_ref());

    // G::read reads x and y
    // edwards::Point::read in sapling-crypto reads y and constructs x from it
    // This is why the following read fails.
    let res = G::read(h.as_ref());
    println!("read res: {:?}", &res);

    match res {
        Ok(p) => {
            let p = p.mul_by_cofactor();

            if !p.is_zero() {
                Some(p)
            } else {
                None
            }
        },
        Err(_) => None,
    }
}

impl<G: Group, W: PedersenWindow> PedersenCRH<G, W> {
    // Create the bases for the Pedersen hashes
    pub fn create_generators() -> Vec<G> {
        let mut generators: Vec<G> = vec![];

        for m in 0..5 {
            use byteorder::{LittleEndian, WriteBytesExt};

            let mut segment_number = [0u8; 4];
            (&mut segment_number[0..4])
                .write_u32::<LittleEndian>(m)
                .unwrap();

            generators.push(find_group_hash(
                &segment_number,
                PEDERSEN_HASH_GENERATORS_PERSONALIZATION,
            ));
        }

        // Check for duplicates, far worse than spec inconsistencies!
        for (i, p1) in generators.iter().enumerate() {
            if p1.is_zero() {
                panic!("Neutral element!");
            }

            for p2 in generators.iter().skip(i + 1) {
                if p1 == p2 {
                    panic!("Duplicate generator!");
                }
            }
        }

        generators
    }

    // Create the exp table for the Pedersen hash generators
    pub fn create_exp_table(generators: &[G]) -> Vec<Vec<Vec<G>>> {
        let mut exp = vec![];
        let window = W::WINDOW_SIZE;

        for g in generators {
            let mut g = g.clone();
            let mut tables = vec![];

            let mut num_bits = 0;
            while num_bits <= <G::ScalarField as PrimeField>::Params::MODULUS_BITS {
                let mut table = Vec::with_capacity(1 << window);

                let mut base = G::zero();

                for _ in 0..(1 << window) {
                    table.push(base.clone());
                    base += &g;
                }

                tables.push(table);
                num_bits += window as u32;

                for _ in 0..window {
                    g.double_in_place();
                }
            }

            exp.push(tables);
        }

        exp
    }
}

impl<G: Group, W: PedersenWindow> FixedLengthCRH for PedersenCRH<G, W> {
    const INPUT_SIZE_BITS: usize = W::WINDOW_SIZE * W::NUM_WINDOWS;
    type Output = G;
    type Parameters = PedersenParameters<G>;

    fn setup<R: Rng>(_rng: &mut R) -> Result<Self::Parameters, Error> {
        let time = timer_start!(|| format!(
            "PedersenCRH::Setup: {} {}-bit windows; {{0,1}}^{{{}}} -> G",
            W::NUM_WINDOWS,
            W::WINDOW_SIZE,
            W::NUM_WINDOWS * W::WINDOW_SIZE
        ));
        let generators = Self::create_generators();
        let exp_table = Self::create_exp_table(&generators);

        timer_end!(time);
        Ok(Self::Parameters {
            generators,
            exp_table,
        })
    }

    fn evaluate(params: &Self::Parameters, input: &[u8]) -> Result<Self::Output, Error> {
        let eval_time = timer_start!(|| "PedersenCRH::Eval");

        let mut bits = BitIterator8::new(input).into_iter();

        let mut result = G::zero();
        let mut generators = params.exp_table.iter();

        loop {
            let mut acc = G::ScalarField::zero();
            let mut cur = G::ScalarField::one();
            let mut chunks_remaining = params.chunks_per_generator();
            let mut encountered_bits = false;

            // Grab three bits from the input
            while let Some(a) = bits.next() {
                encountered_bits = true;

                let b = bits.next().unwrap_or(false);
                let c = bits.next().unwrap_or(false);

                // Start computing this portion of the scalar
                let mut tmp = cur;
                if a {
                    tmp += &cur;
                }
                cur.double_in_place(); // 2^1 * cur
                if b {
                    tmp += &cur;
                }

                // conditionally negate
                if c {
                    tmp = -tmp;
                }

                acc += &tmp;

                chunks_remaining -= 1;

                if chunks_remaining == 0 {
                    break;
                } else {
                    cur.double_in_place(); // 2^2 * cur
                    cur.double_in_place(); // 2^3 * cur
                    cur.double_in_place(); // 2^4 * cur
                }
            }

            if !encountered_bits {
                break;
            }

            let mut table: &[Vec<G>] = &generators.next().expect("we don't have enough generators");
            let window = W::WINDOW_SIZE;
            let window_mask = (1 << window) - 1;

            let mut acc = acc.into_repr();

            let mut tmp = G::zero();

            while !acc.is_zero() {
                let i = (acc.as_ref()[0] & window_mask) as usize;

                tmp += &table[0][i];

                acc.divn(window as u32);
                table = &table[1..];
            }

            result += &tmp;
        }
        timer_end!(eval_time);

        Ok(result)
    }
}

pub fn bytes_to_bits(bytes: &[u8]) -> Vec<bool> {
    let mut bits = Vec::with_capacity(bytes.len() * 8);
    for byte in bytes {
        for i in 0..8 {
            let bit = (*byte >> i) & 1;
            bits.push(bit == 1)
        }
    }
    bits
}

impl<G: Group> Debug for PedersenParameters<G> {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        write!(f, "Pedersen Hash Parameters {{\n")?;
        for (i, g) in self.generators.iter().enumerate() {
            write!(f, "\t  Generator {}: {:?}\n", i, g)?;
        }
        // TODO: exp_table
        write!(f, "}}\n")
    }
}
