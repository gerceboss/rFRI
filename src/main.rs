use lambdaworks_crypto::fiat_shamir::transcript::Transcript;
use lambdaworks_crypto::merkle_tree::backends::types::Keccak256Tree;
use lambdaworks_crypto::merkle_tree::merkle::MerkleTree;
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

// query phase for verification

fn main() {
    println!("Hello, world!");
}
