use crate::algorithms::inexact_synth::grid_problem;


use crate::structs::rings::Float;
use num_traits::One;
type Comp = num_complex::Complex<Float>;

#[test]
pub fn inexact_synth_testing() 
{
    grid_problem(Comp::one(), 0.99,3);
    // TO be removed
    panic!();
}
