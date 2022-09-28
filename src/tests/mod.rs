

// Importing some ring elements
use crate::structs::rings::unimat::UniMat; 
use crate::structs::rings::complex::Complex; 
use crate::structs::rings::dyad::Dyad; 
use crate::structs::rings::domega::DOmega; 

use crate::structs::rings::Constructs;

pub fn basic_identities()-> () {
    
    println!("Running some tests. Results are below!");
    println!("--------------------------------------");

    println!("Test 1:");
    let u: UniMat<Complex> = UniMat::one();
    assert_eq!(u,u*u,"Failed to check that {} =\n {}\n*\n{}",u,u,u);
    
    println!("Test 2:");
    let u: UniMat<Complex> = UniMat::one();
    assert_eq!(u,u.inv(),"Failed to check that {} =\n {}",u,u);

    println!("Test 3:");
    let u: DOmega = Constructs::<DOmega>::one();
    // println!("{}",u);
    // println!("{}",u+u);
    assert_eq!(u,u+u,"Failed to check that {} =\n {}\n+\n{}",u,u,u);
    
    println!("test 4:");
    let u: UniMat<DOmega> = UniMat::one();
    println!("{}",u*u);
    assert_eq!(u,u*u," \n--------------\n Failed to check that \n{} \n= \n{} X \n{}",u,u,u);

    println!("test 5:");
    let u: UniMat<Complex> = UniMat::one();
    assert_eq!(u,u*u," \n--------------\n Failed to check that \n{} \n= \n{} X \n{}",u,u,u);
    
    println!("test 6:");
    let u: DOmega = Constructs::<DOmega>::one();
    // println!("{}",u*u);
    // println!("{}",u);
    assert_eq!(u,u*u," \n--------------\n Failed to check that \n{} \n= \n{} X \n{}",u,u,u);
    
    println!("test 7:");
    let u: Dyad = Constructs::<Dyad>::one();
    // println!("{}",u*u);
    // println!("{}",u);
    assert_eq!(u,u*u," \n--------------\n Failed to check that \n{} \n= \n{} X \n{}",u,u,u);
    
    println!("test 8:");
    let u: Dyad = Constructs::<Dyad>::one();
    // println!("{}",u+u);
    // println!("{}",u);
    assert_eq!(u,u+u," \n--------------\n Failed to check that \n{} \n= \n{} + \n{}",u,u,u);
}
