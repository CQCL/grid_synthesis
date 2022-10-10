
// Importing some ring elements
// use crate::structs::rings::unimat::UniMat; 
// use crate::structs::rings::complex::Complex; 
use std::ops::Neg; 
use std::ops::Add; 
use std::ops::Sub; 
use std::ops::Mul; 
use std::cmp::PartialEq; 
use std::fmt::Display;
use std::fmt::Debug;

use crate::structs::rings::Conj;
use crate::structs::rings::Int;

// use crate::structs::rings::Constructs;
pub fn basic_identities<T>() -> () 
    where T: Copy+Debug+Display,
          T: Add<Output=T>+Mul<Output=T>+Sub<Output=T>+Neg,
          T: PartialEq+From<Int>
{

    println!("Testing basic identities");
    println!("Type is {}", std::any::type_name::<T>());
    println!("--------------------------------------");

    println!("Test 1: 1*1 == 1");
    let u=T::from(1);
    assert_eq!(u,u*u,"Failed to check that {} =\n {}\n*\n{}",u,u,u);


    println!("Test 2: 0+0 == 0");
    let u= T::from(0);
    // println!("{}",u);
    // println!("{}",u+u);
    assert_eq!(u,u+u," \n--------------\n Failed to check that \n{} \n= \n{} + \n{}",u,u,u);

    println!("Test 3: 2+3 == 5");
    let u=T::from(2);
    let v=T::from(3);
    let w=T::from(5);
    assert_eq!(w,v+u,"Failed to check that {} =\n {}\n+\n{}",w,v,u);

    println!("test 4: 1-1 == 0");
    let u= T::from(1);
    let z= T::from(0);
    // // println!("{}",z);
    // // println!("{}",u-u);
    assert_eq!(z,u-u,"Failed to check that {} =\n {}\n-\n{}",z,u,u);

}


// use crate::structs::rings::Constructs;
pub fn basic_identities_with_conj<T>() -> () 
    where T: Copy+Debug+Display,
          T: Add<Output=T>+Mul<Output=T>+Sub<Output=T>+Neg+Conj<T>,
          T: PartialEq+From<Int>
{
    
    println!("Testing conjugation on {}", std::any::type_name::<T>());
    println!("--------------------------------------");
    println!("Test 1: u.conj.conj = u");
    println!("        When u = 2");
    let u= T::from(2);
    assert_eq!(u,u.conj().conj(),"Failed to check that {} =\n {}",u,u);

}
