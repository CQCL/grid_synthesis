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
  
  
// Obviously false code
// fn grid_synth(theta: Angle,epsilon: Error) -> f64 {
//     return 0.0*theta*epsilon
// }
  



use crate::tests::basic_identities;

fn main() {

    // Print text to the console
    println!("------------------------------------------");
    println!("-------------CODE IS RUNNING--------------");
    println!("------------------------------------------");


    basic_identities();

}


