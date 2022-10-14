// Ganesh Puja: 2022-09-13
//
// Cautions: Comments could look pretentious


pub mod structs;

pub mod tests;

// type Angle = f64;
// type Error = f64;

//the compiler suggested this (and wouldn't compile otherwise)
// use crate::structs::rings::dyad::Dyad; 
// use crate::structs::rings::domega::DOmega; 
// use crate::structs::rings::unimat::UniMat; 
// use crate::structs::rings::complex::Complex; 
use crate::structs::rings::local_ring::Local; 
use crate::structs::rings::zroot2::Zroot2; 
use crate::structs::rings::quaternion::Quaternion;
use crate::structs::rings::Fixable;
// use crate::structs::rings::Int;
// use crate::structs::rings::cyclotomic::Cyclotomic; 
  
  
// Obviously false code
// fn grid_synth(theta: Angle,epsilon: Error) -> f64 {
//     return 0.0*theta*epsilon
// }
  
// use crate::tests::basic_identities;
// use crate::tests::basic_identities_with_conj;
// use crate::tests::basic_identities_with_unimat_over;
use crate::tests::testing_complex_rings_vs_quaternions_over;

fn main() {

    // Print text to the console
    println!("------------------------------------------");
    println!("-------------CODE IS RUNNING--------------");
    println!("------------------------------------------");

    // Gotta make some tests
    //
    // basic_identities::<Complex>();
    // basic_identities::<Quaternion<Local<Zroot2>>>();
    // basic_identities_with_conj::<Quaternion<Local<Zroot2>>>();
    //
    // It seems to be going good with integers.
    // Let's implement gates now
    testing_complex_rings_vs_quaternions_over::<Zroot2>();
    
}


