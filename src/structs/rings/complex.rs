// Ring of complex numbers
// This implementation is agnostic of what the complex numbers are defined over
// Basically for any ring R, it makes numbers of the form a+bI, where a and b are in R


// type of floats
// use crate::structs::rings::Float;
use crate::structs::rings::Int;


use crate::structs::rings::Conj; //Conjugation trait

// We bring them in so that we can overload the operators
// Rust must learn how to do arithmetics in our rings
use std::ops::Neg; 
use std::ops::Add; 
use std::ops::Sub; 
use std::ops::Mul; 
use std::cmp::PartialEq; 


// For display
use std::fmt::Result;
use std::fmt::Display;
use std::fmt::Formatter;

#[derive(Copy, Clone, Debug)]
pub struct Complex<T>(pub T,pub T);
// the real part is variable.0
// the imaginary part is variable.1



// Conjugate Complex elements
impl<T> Conj<T> for Complex<T>
where T: Neg<Output=T>
{
    fn conj(self) -> Self {
        Complex(self.0,-self.1)
    }
}


// Negatation on Complex
impl<T> Neg for Complex<T>
where T: Neg<Output=T>
{
    type Output = Self;
    fn neg(self) -> Self {
        Complex(-self.0,-self.1)
    }
}



// Teaching rust how to add Complex elements
impl<T> Add for Complex<T> 
where T: Add<Output=T>
{
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Complex(self.0+other.0,self.1+other.1)
    }
}



// Teaching rust how to subtract Complex elements
impl<T> Sub for Complex<T> 
where T: Sub<Output=T> 
{
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self
        {
            0: self.0-other.0,
            1: self.1-other.1
        }
    }
}


// Teaching rust how to multiply Complex elements
impl<T> Mul for Complex<T> 
where T:Mul<Output=T> + Add<Output=T> + Sub<Output=T>+Copy
{

    type Output = Self;

    fn mul(self, other: Self) -> Self 
    {
        // TODO: Implement the faster multiplication
        // Given here:
        // https://www.embedded.com/digital-signal-processing-tricks-fast-multiplication-of-complex-numbers/
        Self
        {
            0: other.0*self.0 - other.1*self.1,
            1: other.0*self.1 + other.1*self.0 
        }
    }
}

// Nicely display Complex Matrices
impl<T> Display for Complex<T>
where T: Display
{
    fn fmt(&self, f: &mut Formatter) -> Result{
        write!(f,"{}+{}i",self.0, self.1)
    }
}

// Get zero and one as Complex numbers
// NOW DEPRECATED
// impl<T> Constructs<T> for Complex
// {
//     fn zero() -> Self {
//         return Complex(0.0,0.0);
//     }

//     fn one() -> Self {
//         return Complex(1.0,0.0);
//     }
// }


// Teaching rust how to compare these ring elements
impl<T> PartialEq for Complex<T>
where T: PartialEq
{
    fn eq(&self, other: &Self) -> bool {
        return self.0==other.0 && self.1==other.1;
    }
}


impl<T> From<Int> for Complex<T>
where T: From<Int>
{
    fn from(int: Int) -> Self {
        Complex(T::from(int),T::from(0))
    }
}

impl<T> Complex<T>
where T: Mul<Output=T> + Add<Output=T>,
      T: Copy
{
    // reduced square norm of our quaternion algebra
    // Equal to whatever is the output here
    pub fn sqnorm(self) -> T
    {
        return self.0*self.0+self.1*self.1;
    }
}

