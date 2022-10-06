// Ring of complex numbers
// Over floats


// type of floats
type Float = f64;


use crate::structs::rings::Conj; //Conjugation trait
// use crate::structs::rings::Constructs; //Conjugation trait

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
pub struct Complex(pub Float,pub Float);
// the real part is variable.0
// the imaginary part is variable.1



// Conjugate Complex elements
impl<T> Conj<T> for Complex {
    fn conj(self) -> Complex {
        Complex(self.0,-self.1)
    }
}


// Negatation on Complex
impl Neg for Complex{
    type Output = Complex;
    fn neg(self) -> Complex {
        Complex(-self.0,-self.1)
    }
}



// Teaching rust how to add Complex elements
impl Add for Complex {
    type Output = Complex;

    fn add(self, other: Complex) -> Complex {
        Complex(self.0+other.0,self.1+other.1)
    }
}



// Teaching rust how to subtract Complex elements
impl Sub for Complex {
    type Output = Complex;

    fn sub(self, other: Complex) -> Complex {
        self+(-other) //subtraction is just adding the additive inverse
    }
}


// Teaching rust how to multiply Complex elements
impl Mul for Complex {
    type Output = Complex;

    fn mul(self, other: Complex) -> Complex {
        Complex(
            other.0*self.0 - other.1*self.1,
            other.0*self.1 + other.1*self.0 
            )
    }
}

// Nicely display Complex Matrices
impl Display for Complex
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
impl PartialEq for Complex
{
    fn eq(&self, other: &Self) -> bool {
        return self.0==other.0 && self.1==other.1;
    }
}


impl From<i32> for Complex {
    fn from(int: i32) -> Self {
        Complex(int.try_into().unwrap(),0.0)
    }
}

impl From<Float> for Complex {
    fn from(int: Float) -> Self {
        Complex(int,0.0)
    }
}
