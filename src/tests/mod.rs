

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

// use crate::structs::rings::Constructs;
pub trait CheckRingIdentities<T>
where T: Copy+Debug+Display,
      T: Add<Output=T>+Conj<T>+Mul<Output=T>+Sub<Output=T>+Neg,
      T: PartialEq+From<i32>
{
    fn basic_identities()-> () {

        println!("Running some tests. Results are below!");
        println!("--------------------------------------");

        println!("Test 1:");
        let u=T::from(1);
        assert_eq!(u,u*u,"Failed to check that {} =\n {}\n*\n{}",u,u,u);

        println!("Test 2:");
        let u= T::from(2);
        assert_eq!(u,u.conj().conj(),"Failed to check that {} =\n {}",u,u);

        println!("Test 3:");
        let u= T::from(0);
        // println!("{}",u);
        // println!("{}",u+u);
        assert_eq!(u,u+u," \n--------------\n Failed to check that \n{} \n= \n{} + \n{}",u,u,u);

        println!("Test 4:");
        let u=T::from(2);
        let v=T::from(3);
        let w=T::from(5);
        assert_eq!(w,v+u,"Failed to check that {} =\n {}\n+\n{}",w,v,u);

        // println!("test 5:");
        let u= T::from(1);
        let z= T::from(0);
        // // println!("{}",z);
        // // println!("{}",u-u);
        assert_eq!(z,u-u,"Failed to check that {} =\n {}\n-\n{}",z,u,u);

        // println!("test 6:");
        // let u: DOmega = DOmega::from(1);
        // // println!("{}",u*u);
        // // println!("{}",u);
        // assert_eq!(u,u*u," \n--------------\n Failed to check that \n{} \n= \n{} X \n{}",u,u,u);

        // println!("test 7:");
        // let u: Dyad = Dyad::from(1);
        // // println!("{}",u*u);
        // // println!("{}",u);
        // assert_eq!(u,u*u," \n--------------\n Failed to check that \n{} \n= \n{} X \n{}",u,u,u);

        // println!("test 8:");
        // let u: Dyad = Dyad::from(1);
        // let v: Dyad = Dyad::from(2);
        // // println!("{}",u+u);
        // // println!("{}",u);
        // assert_eq!(v,u+u," \n--------------\n Failed to check that \n{} \n= \n{} + \n{}",v,u,u);


        // println!("Test 9:");
        // let u: DOmega = DOmega::from(1);
        // let z: DOmega = DOmega::from(0);
        // // println!("{}",z);
        // // println!("{}",u-u);
        // assert_eq!(z,u-u,"Failed to check that {} =\n {}\n-\n{}",z,u,u);


    }
}
