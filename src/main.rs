// Ganesh Puja: 2022-09-13
//
// Cautions: Comments could look pretentious


pub mod structs;

type Angle = f64;
type Error = f64;

//the compiler suggested this (and wouldn't compile otherwise)
// use crate::structs::rings::Dyad; 
use crate::structs::rings::unimat::UniMat; 

// Obviously false code
fn grid_synth(theta: Angle,epsilon: Error) -> f64 {
    return 0.0*theta*epsilon
}


fn main() {

    println!("{:?}",d.conj());
    // Print text to the console
    println!("------------------------------------------");
    println!("-------------CODE IS RUNNING--------------");
    println!("------------------------------------------");


    let mut u=UniMat{
        u: 1,
        t: 0
    };


}


