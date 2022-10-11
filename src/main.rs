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
use crate::structs::rings::complex::Complex; 
use crate::structs::rings::local_ring::Local; 
use crate::structs::rings::zroot2::Zroot2; 
use crate::structs::rings::quaternion::Quaternion;
// use crate::structs::rings::Int;
// use crate::structs::rings::cyclotomic::Cyclotomic; 
  
  
// Obviously false code
// fn grid_synth(theta: Angle,epsilon: Error) -> f64 {
//     return 0.0*theta*epsilon
// }
  
use crate::tests::basic_identities;
use crate::tests::basic_identities_with_conj;
// use crate::tests::basic_identities_with_unimat_over;

fn main() {

    // Print text to the console
    println!("------------------------------------------");
    println!("-------------CODE IS RUNNING--------------");
    println!("------------------------------------------");

    // basic_identities::<Complex>();
    // basic_identities::<Zroot2>();
    // basic_identities::<Local<Zroot2>>();
    // basic_identities_with_conj::<Complex>();
    // basic_identities_with_conj::<Zroot2>();
    // basic_identities_with_unimat_over::<Complex>();
    // basic_identities_with_unimat_over::<Local<Cyclotomic>>();
    basic_identities::<Quaternion<Local<Zroot2>>>();
    basic_identities_with_conj::<Quaternion<Local<Zroot2>>>();

}


