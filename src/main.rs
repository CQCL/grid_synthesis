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
// use crate::structs::rings::Int;
// use crate::structs::rings::cyclotomic::Cyclotomic; 
  
  
// Obviously false code
// fn grid_synth(theta: Angle,epsilon: Error) -> f64 {
//     return 0.0*theta*epsilon
// }
  
// use crate::tests::basic_identities;
// use crate::tests::basic_identities_with_conj;
// use crate::tests::basic_identities_with_unimat_over;
// use crate::tests::testing_that_localizable_rings_work;

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

    
    let h = Quaternion::<Local::<Zroot2>>
    {
        0:  Local::<Zroot2>::from(0),
        1:  Local::<Zroot2>
            {
                num: Zroot2::from(1),
                log_den:1
            },
        2:  Local::<Zroot2>::from(0),
        3:  Local::<Zroot2>
            {
                num: Zroot2::from(1),
                log_den:1
            },
    };
    println!("{}")


    let h = Quaternion::<Local::<Zroot2>>
    {
        0:  Local::<Zroot2>::from(0),
        1:  Local::<Zroot2>
            {
                num: Zroot2::from(1),
                log_den:1
            },
        2:  Local::<Zroot2>::from(0),
        3:  Local::<Zroot2>
            {
                num: Zroot2::from(1),
                log_den:1
            },
    };
    
    
}


