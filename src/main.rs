use lambdaworks_crypto::merkle_tree::backends::types::Keccak256Tree;
use lambdaworks_crypto::merkle_tree::merkle::MerkleTree;
use lambdaworks_math::{
    // fft::polynomial::FFTPoly,
    fft::polynomial::FFTPoly,
    field::{
        element::FieldElement,
        traits::{IsFFTField, IsField},
    },
    polynomial::Polynomial,
    traits::ByteConversion,
};

// Merkle Trees configuration

// Security of both hashes should match
pub type FriMerkleTreeBackend<F> = Keccak256Tree<F>; // should be tree
pub type FriMerkleTree<F> = MerkleTree<FriMerkleTreeBackend<F>>;

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
// pub commit_phase(number_of_layers:usize,
//     p_0: Polynomial<FieldElement<F>>,
//     transcript: &mut T,
//     coset_offset: &FieldElement<F>,
//     domain_size: usize,)->(Vec<FRILayer<F>>,FieldElement<F>){

// }

// decommitment phase

// query phase for verification

fn main() {
    println!("Hello, world!");
}
