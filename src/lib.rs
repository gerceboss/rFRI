use lambdaworks_crypto::fiat_shamir::transcript::Transcript;
use lambdaworks_crypto::merkle_tree::backends::types::Keccak256Tree;
use lambdaworks_crypto::merkle_tree::merkle::MerkleTree;
use lambdaworks_crypto::merkle_tree::proof::Proof;
use lambdaworks_math::{
    // fft::polynomial::FFTPoly,
    fft::polynomial::FFTPoly,
    field::{
        element::FieldElement,
        traits::{IsFFTField, IsField, IsPrimeField},
    },
    polynomial::Polynomial,
    traits::{ByteConversion, Deserializable, Serializable},
};

// Merkle Trees configuration

// Security of both hashes should match
pub type FriMerkleTreeBackend<F> = Keccak256Tree<F>; // should be tree
pub type FriMerkleTree<F> = MerkleTree<FriMerkleTreeBackend<F>>;

// for 256 bits it is 32, for 512 it should be 64
pub const COMMITMENT_SIZE: usize = 32;
pub type Commitment = [u8; COMMITMENT_SIZE];

/// Uses randomness from the transcript to create a FieldElement
/// One bit less than the max used by the FieldElement is used as randomness. For StarkFields, this would be 251 bits randomness.
/// Randomness is interpreted as limbs in BigEndian, and each Limb is ordered in BigEndian
pub fn transcript_to_field<F: IsPrimeField, T: Transcript>(transcript: &mut T) -> FieldElement<F>
where
    FieldElement<F>: lambdaworks_math::traits::ByteConversion,
{
    let mut randomness = transcript.challenge();
    randomness_to_field(&mut randomness)
}

/// Transforms some random bytes to a field
/// Slicing the randomness to one bit less than what the max number of the field is to ensure each random element has the same probability of appearing
fn randomness_to_field<F: IsPrimeField>(randomness: &mut [u8; 32]) -> FieldElement<F>
where
    FieldElement<F>: ByteConversion,
{
    let random_bits_required = F::field_bit_size() - 1;
    let random_bits_created = randomness.len() * 8;
    let mut bits_to_clear = random_bits_created - random_bits_required;

    let mut i = 0;
    while bits_to_clear >= 8 {
        randomness[i] = 0;
        bits_to_clear -= 8;
        i += 1;
    }

    let pre_mask: u8 = 1u8.checked_shl(8 - bits_to_clear as u32).unwrap_or(0);
    let mask: u8 = pre_mask.wrapping_sub(1);
    randomness[i] &= mask;

    FieldElement::from_bytes_be(randomness).unwrap()
}

pub fn transcript_to_usize<T: Transcript>(transcript: &mut T) -> usize {
    const CANT_BYTES_USIZE: usize = (usize::BITS / 8) as usize;
    let value = transcript.challenge()[..CANT_BYTES_USIZE]
        .try_into()
        .unwrap();
    usize::from_be_bytes(value)
}
// folding divides the polynomail in even and odd coefficientsa and then combines them using a constant beta to reduce the final degree of the polynomial by 2.
// Hence this helps achieve reduce the polynomail to a constant value in log2(N) setps
pub fn fold_polynomial<F>(
    poly: &Polynomial<FieldElement<F>>,
    beta: &FieldElement<F>,
) -> Polynomial<FieldElement<F>>
where
    F: IsField,
{
    let coeffs = poly.coefficients();
    let even_coeffs: Vec<FieldElement<F>> = coeffs.iter().step_by(2).cloned().collect();
    let odd_coeffs_multiply_beta: Vec<FieldElement<F>> = coeffs
        .iter()
        .skip(1)
        .step_by(2)
        .map(|x| (x.clone() * beta))
        .collect();

    let (even_poly, odd_poly) = Polynomial::pad_with_zero_coefficients(
        &Polynomial::new(&even_coeffs),
        &Polynomial::new(&odd_coeffs_multiply_beta),
    );
    even_poly + odd_poly
}

// create layers
#[derive(Clone)]
pub struct FRILayer<F>
where
    F: IsField,
    FieldElement<F>: ByteConversion,
{
    pub evals: Vec<FieldElement<F>>,
    pub merkle_tree: FriMerkleTree<F>, //for committing our evals
    pub coset_offset: FieldElement<F>, //Generally added for evaluation in coset FFT
    pub domain_size: usize,            // to get the size to start with
}

// consider adding how the merkle tree is actually created
impl<F> FRILayer<F>
where
    F: IsField + IsFFTField,
    FieldElement<F>: ByteConversion,
{
    pub fn new(
        poly: &Polynomial<FieldElement<F>>,
        coset_offset: &FieldElement<F>,
        domain_size: usize,
    ) -> Self {
        let evals = poly
            .evaluate_offset_fft(1, Some(domain_size), coset_offset)
            .unwrap();

        let merkle_tree = FriMerkleTree::build(&evals);

        Self {
            evals,
            merkle_tree,
            coset_offset: coset_offset.clone(),
            domain_size,
        }
    }
}

// Commitment phase

// The commit phase will give us a vector of layers and the final value of the FRI protocol (when we get to a degree zero polynomial)
pub fn commit_phase<F: IsField + IsFFTField, T: Transcript>(
    number_of_layers: usize,
    p_0: Polynomial<FieldElement<F>>,
    transcript: &mut T,
    coset_offset: &FieldElement<F>,
    domain_size: usize,
) -> (Vec<FRILayer<F>>, FieldElement<F>)
where
    FieldElement<F>: ByteConversion,
{
    let mut domain_size = domain_size;
    let mut fri_layer_list = Vec::with_capacity(number_of_layers);

    let mut current_layer = FRILayer::new(&p_0, coset_offset, domain_size);
    let mut coset_offset = coset_offset.clone();
    let mut current_poly = p_0;

    // Send commitment: [p₀]
    transcript.append(&current_layer.merkle_tree.root);

    for _ in 1..number_of_layers {
        // Receive challenge (beta,k-1) (beta k-1 th)
        let beta = transcript_to_field(transcript);
        coset_offset = coset_offset.square();
        domain_size /= 2;

        // Compute polynomial for the new layer and domain
        current_poly = fold_polynomial(&current_poly, &beta);
        // NOTE:  folding can be done any number of times in powers of 2.
        // polynomial can be foled by 4 , by 8 , by 16 ...
        current_layer = FRILayer::new(&current_poly, &coset_offset, domain_size);
        fri_layer_list.push(current_layer.clone());

        // Send commitment: [pₖ]
        transcript.append(&current_layer.merkle_tree.root);
    }

    // Receive challenge: beta n-1
    let beta = transcript_to_field(transcript);

    let final_poly = fold_polynomial(&current_poly, &beta);

    let last_value = final_poly
        .coefficients()
        .get(0)
        .unwrap_or(&FieldElement::zero())
        .clone();

    // Send value: pₙ
    transcript.append(&last_value.to_bytes_be());

    (fri_layer_list, last_value)
}

// decommitment phase
// implement serialisation and deserialisation
pub struct FriDecommitment<F: IsPrimeField> {
    pub layers_auth_paths_sym: Vec<Proof<Commitment>>,
    pub layers_evaluations_sym: Vec<FieldElement<F>>,
    pub layers_auth_paths: Vec<Proof<Commitment>>,
    pub layers_evaluations: Vec<FieldElement<F>>,
}

// query phase and verification
pub fn fri_query_phase<F, T>(
    queries: i32, // should come from number of constraints in an AIR
    domain_size: usize,
    fri_layers: &Vec<FRILayer<F>>,
    transcript: &mut T,
) -> (Vec<FriDecommitment<F>>, Vec<usize>)
where
    F: IsFFTField,
    T: Transcript,
    FieldElement<F>: ByteConversion,
{
    if !fri_layers.is_empty() {
        let number_of_queries = queries;
        let iotas = (0..number_of_queries)
            .map(|_| transcript_to_usize(transcript) % domain_size)
            .collect::<Vec<usize>>();
        let query_list = iotas
            .iter()
            .map(|iota_s| {
                // Receive challenge  (iota_s)
                let mut layers_auth_paths_sym = vec![];
                let mut layers_evaluations_sym = vec![];
                let mut layers_evaluations = vec![];
                let mut layers_auth_paths = vec![]; //paths for authentication of inclusion

                for layer in fri_layers {
                    // symmetric element
                    let index = iota_s % layer.domain_size;
                    let index_sym = (iota_s + layer.domain_size / 2) % layer.domain_size;
                    let evaluation_sym = layer.evals[index_sym].clone();
                    let auth_path_sym = layer.merkle_tree.get_proof_by_pos(index_sym).unwrap();
                    let evaluation = layer.evals[index].clone();
                    let auth_path = layer.merkle_tree.get_proof_by_pos(index).unwrap();
                    layers_auth_paths_sym.push(auth_path_sym);
                    layers_evaluations_sym.push(evaluation_sym);
                    layers_evaluations.push(evaluation);
                    layers_auth_paths.push(auth_path);
                }

                FriDecommitment {
                    layers_auth_paths_sym,
                    layers_evaluations_sym,
                    layers_evaluations,
                    layers_auth_paths,
                }
            })
            .collect();

        (query_list, iotas)
    } else {
        (vec![], vec![])
    }
}

#[cfg(test)]
mod tests {
    use super::fold_polynomial;
    use lambdaworks_math::field::element::FieldElement;
    use lambdaworks_math::field::fields::u64_prime_field::U64PrimeField;
    use lambdaworks_math::polynomial::Polynomial;
    const MODULUS: u64 = 293;
    type F = FieldElement<U64PrimeField<MODULUS>>;

    #[test]
    fn test_fold() {
        let p0 = Polynomial::new(&[
            F::new(3),
            F::new(1),
            F::new(2),
            F::new(7),
            F::new(3),
            F::new(5),
        ]);
        let beta = F::new(4);
        let p1 = fold_polynomial(&p0, &beta);
        assert_eq!(p1, Polynomial::new(&[F::new(7), F::new(30), F::new(23),]));

        let gamma = F::new(3);
        let p2 = fold_polynomial(&p1, &gamma);
        assert_eq!(p2, Polynomial::new(&[F::new(97), F::new(23),]));

        let delta = F::new(2);
        let p3 = fold_polynomial(&p2, &delta);
        assert_eq!(p3, Polynomial::new(&[F::new(143)]));
        assert_eq!(p3.degree(), 0);
    }
}
