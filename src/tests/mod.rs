// Importing some ring elements
use crate::structs::unimat::UniMat; 
use crate::structs::rings::Localizable; 
// use crate::structs::rings::Fixable; 
use std::ops::Neg; 
use std::ops::Add; 
use std::ops::Sub; 
use std::ops::Mul; 
use std::cmp::PartialEq; 
use std::fmt::Display;
use std::fmt::Debug;

use crate::structs::rings::Conj;
use crate::structs::rings::local_ring::Local;
use crate::structs::rings::quaternion::Quaternion;
use crate::structs::rings::Int;

// use crate::structs::rings::Constructs;
pub fn basic_identities<T>() -> () 
    where T: Copy+Debug+Display,
          T: Add<Output=T>+Mul<Output=T>+Sub<Output=T>+Neg<Output=T>,
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

    println!("Test 5: 2-3 == -1");
    let u=T::from(2);
    let v=T::from(3);
    let w=T::from(-1);
    assert_eq!(w,u-v,"Failed to check that {} =\n {}\n-\n{}",w,u,v);
    
    println!("Test 6: -2*2 == -4");
    let u=T::from(2);
    let v=T::from(-2);
    let w=T::from(-4);
    assert_eq!(w,u*v,"Failed to check that {} =\n{}*\n{}",w,v,u);
    
    println!("Test 6: -3*-3*-3 == -27");
    let u=T::from(-3);
    let w=T::from(-27);
    assert_eq!(w,u*u*u,"Failed to check that {} ={}*\n{}*\n{}",w,u,u,u);
}


// use crate::structs::rings::Constructs;
pub fn basic_identities_with_conj<T>() -> () 
    where T: Copy+Debug+Display,
          T: Add<Output=T>+Mul<Output=T>+Sub<Output=T>+Neg<Output=T>+Conj<T>,
          T: PartialEq+From<Int>
{
    
    println!("Testing conjugation on {}", std::any::type_name::<T>());
    println!("--------------------------------------");
    println!("Test 1: u.conj.conj = u");
    println!("        When u = 2");
    let u= T::from(2);
    assert_eq!(u,u.conj().conj(),"Failed to check that {} =\n {}",u,u);

}



pub fn basic_identities_with_unimat_over<T>() -> () 
    where T: Copy+Debug+Display,
          T: Add<Output=T>+Mul<Output=T>+Sub<Output=T>+Neg<Output=T>+Conj<T>,
          T: PartialEq+From<Int>
{
    
    println!("Testing UniMat over {}", std::any::type_name::<T>());
    println!("--------------------------------------");
    println!("Test 1: Id.Id = Id");
    let u = UniMat::<T>::one();
    assert_eq!(u,u*u,"Failed to check that {} =\n {}\n*\n{}",u,u,u);

}


pub fn testing_localizable_rings<T>() -> ()
where T: Copy+Debug+Display,
      T: Add<Output=T>+Mul<Output=T>+Sub<Output=T>+Neg<Output=T>+Conj<T>,
      T: PartialEq+From<Int>+Localizable
{

    println!("--------------------");
    let mut b = T::from(4);
    println!("{}",b);
    b= b.perform_n_multiplications(1);
    println!("{}",b);
    b= b.perform_n_multiplications(1);
    println!("{}",b);
    b= b.perform_n_multiplications(1);
    println!("{}",b);
    b= b.perform_n_multiplications(1);
    println!("{}",b);
    let i :Int;
    (b , i) = b.reduce_by_dividing();
    println!("{}, divided {} times",b,i);

}

pub fn testing_complex_rings_vs_quaternions_over<T>() 
where T: Copy+Debug+Display,
      T: Add<Output=T>+Mul<Output=T>+Sub<Output=T>+Neg<Output=T>+Conj<T>,
      T: PartialEq+From<Int>+Localizable
{
    let h = Quaternion::<Local::<T>>
    {
        0:  Local::<T>::from(0),
        1:  Local::<T>
            {
                num: T::from(1),
                log_den:1
            },
        2:  Local::<T>::from(0),
        3:  Local::<T>
            {
                num: T::from(1),
                log_den:1
            },
    };

    // println!("{}",h);
    // h=h*h*h;

    // println!("{}",h);
    // println!("{}", h.rsqnorm());

    // Warning: THis is the T gate only when 
    // the type T used here is Zroot2
    let t = Quaternion::<Local::<T>>
    {
        0:  Local::<T>
            {
                num: T::from(1)+(T::from(1)).perform_n_multiplications(1),
                log_den:1
            },
        1:  Local::<T>
            {
                num: T::from(1),
                log_den:1
            },
        2:  Local::<T>::from(0),
        3:  Local::<T>::from(0)
    };


    let t2 = t*t;
    let m= h*t*h*t*t*t*h*t*h*t*t*t*h*t*t*t*t*t;

    println!("Test 1: Checking norm is equal to sum of sqnorms");
    assert_eq!(h.rsqnorm(),h.z().sqnorm()+h.w().sqnorm());
    assert_eq!(t.rsqnorm(),t.z().sqnorm()+t.w().sqnorm());
    assert_eq!(t2.rsqnorm(),t2.z().sqnorm()+t2.w().sqnorm());
    assert_eq!(m.rsqnorm(),m.z().sqnorm()+m.w().sqnorm());
}
