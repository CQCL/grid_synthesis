use crate::structs::rings::zroot2::Zroot2;
use crate::structs::rings::local_ring::Local;
use num_complex::Complex;
use crate::structs::sunimat::UniMat;
use crate::structs::rings::special_values::{mu_8, onebyroot2comp, sqrt2, sqrtminus1};

type Comp = Complex::<Local::<Zroot2>>;
type Mat = UniMat::<Comp>;


use num_traits::Zero;

#[test]
// Found on 12-12-2022
// A bit thorny
pub fn bens_matrix_test()
{
    let hadamard = Mat{
        u: sqrtminus1() / sqrt2(),
        t: sqrtminus1() / sqrt2(),
    };

    let t_gate = Mat{
        u: mu_8(),
        t: Comp::zero(),
    };

    let ht =  hadamard * t_gate;
    let hththt = hadamard * t_gate * hadamard * t_gate * hadamard * t_gate;

    println!("{}", hadamard);
    println!("{}", t_gate);



}
