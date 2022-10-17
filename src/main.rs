// Ganesh Puja: 2022-09-13
//
// Cautions: Comments could look pretentious


pub mod structs;
pub mod tests;
pub mod algorithms;

// type Angle = f64;
// type Error = f64;

//the compiler suggested this (and wouldn't compile otherwise)
// use crate::structs::rings::dyad::Dyad; 
// use crate::structs::rings::domega::DOmega; 
// use crate::structs::rings::unimat::UniMat; 
use crate::structs::rings::complex::Complex; 
use crate::structs::rings::local_ring::Local; 
use crate::structs::rings::zroot2::Zroot2; 
use crate::structs::rings::quaternion::Quaternion;
// use crate::structs::rings::Fixable;
// use crate::structs::rings::Localizable;
// use crate::structs::rings::Int;
// use crate::structs::rings::cyclotomic::Cyclotomic; 


use crate::algorithms::exact_synth::exact_synth_given_norm_1;

// use crate::tests::basic_identities;
// use crate::tests::basic_identities_with_conj;
// use crate::tests::basic_identities_with_unimat_over;
// use crate::tests::testing_complex_rings_vs_quaternions_over;

// Better looking code
type Loc = Local<Zroot2>;
type Quat = Quaternion<Loc>;
type Comp = Complex<Loc>;

fn main() {

    // Print text to the console
    println!("------------------------------------------");
    println!("-------------CODE IS RUNNING--------------");
    println!("------------------------------------------");

    let h = Quat::h_gate();
    let t = Quat::t_gate();
    let omega = Comp::mu_8();
    let one = Comp::from(1);

    let g = Comp::quat_conj_transpose_second(one+omega,one-omega);
    println!("{}", g);
    println!("{}", g.w().sqnorm());
    

    let input = h*t*t*h*t*t*t*t*h*t*t*h*t*t*t*t*t*t*h*t*t*t*t*t*h*t*t*t*h*t*t*t*h*t*t*t*t*t*h*t*t*t*h;

    exact_synth_given_norm_1(input);

}


