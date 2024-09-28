use lambdaworks_math::field::{element::FieldElement, traits::IsField};
use lambdaworks_math::polynomial::Polynomial;



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

// commitment phase

// decommitment phase

// query phase for verification

fn main() {
    println!("Hello, world!");
}