// This rust module contains all our rings etc.
// Some of these O(1) operations about multiplication and addition can be optimized a little bit
// The DOmega ring can be stored directly without any underlying Dyad implentation, I feel
pub type Int = i64;
pub type Float = f64;

// Supertrait
// They become various conjugation operators on our rings
pub trait Conj{
    fn conj(self) -> Self;
}

// Supertrait (deprecated now, use from::(0) etc instead)
// They return zero and one elements of our rings
// pub trait Constructs<T>{
//     fn zero() -> Self;
//     fn one() -> Self;
// }

// Supertrait
// Allows you to define a possible prime ideal for localization
// See Zroot2 for an implementation
pub trait Localizable{
    // one should be able to check if the divisibility by ideal exists
    fn is_divisible(self) -> bool;
    // if it is divisible, we perform the division
    // We return the number of times we divided
    fn reduce_by_dividing(self) -> (Self,Int)
    where Self: Sized;

    // Multiply the ideal generator by the number of times given by n
    // WARNING: Bit Overflow may occur here if we perform too much multiplication
    fn perform_n_multiplications(self,_:Int) -> Self;
    

    // The following two functions should have their own trait.
    // TODO
    
    // Should return norm
    fn norm(self) -> Int;
    
    // Should check if norm is a unit 
    fn is_unit(self) -> bool;
}

// Supertrait
// Allows addition and multiplication in a local ring
// See local_ring.rs for an implementation
pub trait Fixable{
    
    // Reduce a number to it's lowest form, so that we can extract its p-adic valuation
    fn fix(self) -> Self;
    
    // Return the log_den value to be used by something else
    fn logden(self) -> Int;
    
    // Return sqrt2 for zroot2
    // Or in case of any other localizable ring
    // return something other than sqrt2
    fn local_gen() -> Self;
}

pub mod zroot2;
pub mod local_ring;
pub mod complex;
// pub mod quaternion;

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
    for i in 1..(p-1)
    {
        out=out*t;
    }
    return out;
}



