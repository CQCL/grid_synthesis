// Ganesh Puja: 2022-09-13
//
// Cautions: Comments could look pretentious


pub mod structs;
pub mod tests;
pub mod algorithms;


// Num traits
use num_traits::Num;
use num_traits::Zero;
use num_traits::One;
use num_traits::NumCast;
use num_traits::FromPrimitive;

// type Angle = f64;
// type Error = f64;

//the compiler suggested this (and wouldn't compile otherwise)
// use crate::structs::rings::dyad::Dyad; 
// use crate::structs::rings::domega::DOmega; 
use crate::structs::sunimat::UniMat; 
use crate::structs::rings::complex::Complex; 
use crate::structs::rings::local_ring::Local; 
use crate::structs::rings::zroot2::Zroot2; 
use crate::structs::rings::quaternion::Quaternion;
// use crate::structs::rings::Fixable;
// use crate::structs::rings::Localizable;
use crate::structs::rings::Int;
// use crate::structs::rings::Conj;
// use crate::structs::rings::cyclotomic::Cyclotomic; 
use crate::structs::rings::pow;

use crate::algorithms::exact_synth::exact_synth_given_norm_1;

// use crate::tests::basic_identities;
// use crate::tests::basic_identities_with_conj;
// use crate::tests::basic_identities_with_unimat_over;
// use crate::tests::testing_complex_riags_vs_quaternions_over;
// use crate::tests::should_break_arithmetic_26_10_2022;
// use crate::tests::broke_arithmetic_until_26_10_2022;
// use crate::tests::doesnt_break_matrices_27_10_2022;
// use crate::tests::apply_gate_string_to_states_and_check_output;

// Better looking code
type Loc = Local<Zroot2>;
type Comp = Complex<Loc>;
type Mat = UniMat<Comp>;
// type Quat = Quaternion<Loc>;


fn main() {

    // Print text to the console
    println!("------------------------------------------");
    println!("-------------CODE IS RUNNING--------------");
    println!("------------------------------------------");

    let omega = Comp::mu_8();
    let onebyroot2 = Comp::onebyroot2();
    let root2 = Comp::root2();
    let one = Comp::one();
    let zero = Comp::zero();
    
    let u1 = ( one+omega )*onebyroot2*onebyroot2;
    let t1 = ( one-omega )*onebyroot2*onebyroot2; 

    let mut g = Mat{
        u : u1,
        t : t1
    };


    // g=g*g*g*g*g*g; //*g*g*g*g*g;

    // let (gate_sequence,output)  = exact_synth_given_norm_1(g);

    // println!("{}", gate_sequence);

}
