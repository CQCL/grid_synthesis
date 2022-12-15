#![allow(warnings)]

pub type Int = i64;
pub type LogDepInt = i32;
pub type Float = f64;

// Supertrait
// They become various conjugation operators on our rings
pub trait Conj{
    fn conj(&self) -> Self;
}

// Supertrait
// They become various conjugation operators on our rings
pub trait LocalizableNorm{
    fn norm(&self) -> Int;
}


// Supertrait
// Allows you to define a possible prime ideal for localization
// See Zroot2 for an implementation
// Also see Int for an implementation
pub trait Localizable{
    // one should be able to check if the divisibility by ideal exists
    fn is_divisible(self) -> bool;
    // if it is divisible, we perform the division
    // We return the number of times we divided
    fn reduce_by_dividing(self) -> (Self,LogDepInt)
    where Self: Sized;

    // Multiply the ideal generator by the number of times given by n
    // WARNING: Bit Overflow may occur here if we perform too much multiplication
    fn perform_n_multiplications(self,_:LogDepInt) -> Self;
    

}
    // The following two functions should have their own trait.
    // TODO
    
    // Should return norm
    // fn norm(self) -> Int;
    
    // Should check if norm is a unit 
    // fn is_unit(self) -> bool;





pub mod zroot2;
pub mod zomega;
pub mod int_localization;
pub mod local_ring;
pub mod special_values;
pub mod quaternion;

use std::ops::Mul;

pub fn pow<T>(t: T, p: Int) -> T
where T: Mul<Output=T>+Copy
{
    if p<=0
    {
        panic!("Power on a non-positive integer");
    }

    let mut out = t;
    // TODO: Smarter implementation using 
    // binary expansion of p is possible here
    for i in 1..p
    {
        out=out*t;
    }
    return out;
}


use std::ops::Neg;
use num_traits::Num;

// Implementing Conjuation for complex numbers
impl<T> Conj for num_complex::Complex<T>
where T: Neg<Output=T>,
      T: Copy,
      T: Num
{
    fn conj(&self) -> Self
    {
        return Self
        {
            re: self.re,
            im: -self.im,
        };
    }
}
