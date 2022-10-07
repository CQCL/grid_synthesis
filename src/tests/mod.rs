

// Importing some ring elements
use crate::structs::rings::unimat::UniMat; 
use crate::structs::rings::complex::Complex; 
use crate::structs::rings::dyad::Dyad; 
use crate::structs::rings::domega::DOmega; 

// use crate::structs::rings::Constructs;

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
    let u: DOmega = DOmega::from(1);
    // println!("{}",u);
    // println!("{}",u+u);
    assert_eq!(u,u*u," \n--------------\n Failed to check that \n{} \n= \n{} + \n{}",u,u,u);
    
    println!("test 4: TODO");
    // let u: UniMat<DOmega> = UniMat::one();
    // println!("{}",u*u);
    // assert_eq!(u,u*u," \n--------------\n Failed to check that \n{} \n= \n{} X \n{}",u,u,u);

    println!("test 5:");
    let u: UniMat<Complex> = UniMat::one();
    assert_eq!(u,u*u," \n--------------\n Failed to check that \n{} \n= \n{} X \n{}",u,u,u);
    
    println!("test 6:");
    let u: DOmega = DOmega::from(1);
    // println!("{}",u*u);
    // println!("{}",u);
    assert_eq!(u,u*u," \n--------------\n Failed to check that \n{} \n= \n{} X \n{}",u,u,u);
    
    println!("test 7:");
    let u: Dyad = Dyad::from(1);
    // println!("{}",u*u);
    // println!("{}",u);
    assert_eq!(u,u*u," \n--------------\n Failed to check that \n{} \n= \n{} X \n{}",u,u,u);
    
    println!("test 8:");
    let u: Dyad = Dyad::from(1);
    let v: Dyad = Dyad::from(2);
    // println!("{}",u+u);
    // println!("{}",u);
    assert_eq!(v,u+u," \n--------------\n Failed to check that \n{} \n= \n{} + \n{}",v,u,u);


    println!("Test 9:");
    let u: DOmega = DOmega::from(1);
    let z: DOmega = DOmega::from(0);
    // println!("{}",z);
    // println!("{}",u-u);
    assert_eq!(z,u-u,"Failed to check that {} =\n {}\n-\n{}",z,u,u);

    println!("Test 10:");
    let u: DOmega = DOmega::from(2);
    let v: DOmega = DOmega::from(3);
    let w: DOmega = DOmega::from(5);
    assert_eq!(w,v+u,"Failed to check that {} =\n {}\n+\n{}",w,v,u);

}
