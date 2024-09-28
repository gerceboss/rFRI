# rFRI

FRI implementation in Rust for learning purposes. Verifier related function are not yet implemented but for verification he needs to perform both the inclusion proofs (to see all values belong to Merkle trees) and that each layer is obtained from the previous one until we reach degree zero.

Implementation includes:

- Folding of Polynomial
- Commitment phase
- Decommitment
- Query phase
- Test for folding of Polynomial
