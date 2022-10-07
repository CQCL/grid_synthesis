// Quadratic number field with root 2
// This struct assumes that you are going to localize it 
// at sqrt(2)
pub struct Zroot2(pub i64,pub i64); //a+b\sqrt(2)

                
// use crate::structs::rings::Conj; //Conjugation trait
// use crate::structs::rings::Constructs; //Conjugation trait
use crate::structs::rings::Localizable; //Conjugation trait


// We bring them in so that we can overload the operators
// Rust must learn how to do arithmetics in our rings
// use std::ops::Neg; 
use std::ops::Add; 
// use std::ops::Sub; 
// use std::ops::Mul; 
// use std::cmp::PartialEq; 

// For display
use std::fmt::Result;
use std::fmt::Display;
use std::fmt::Formatter;
                
// Rust must know how to diplay elements of this ring
impl Display for Zroot2{
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        if self.1<0 {
            write!(f, "{}{}\\sqrt{{2}}", self.0, self.1)
        }
        else
        {
            write!(f, "{}+{}\\sqrt{{2}}", self.0, self.1)
        }
    }
}

// Teaching rust how to add Zroot2 elements
impl Add for Zroot2 {
    type Output = Zroot2;

    fn add(self, other: Zroot2) -> Zroot2 {
        Zroot2(self.0 + other.0, self.1 + other.1)
    }
}

// Conjugate Zroot2 elements
// The only non-trivial element in the Galois group
impl Zroot2 {
    pub fn conj(self) -> Zroot2 {
        Zroot2(self.0,-self.1)
    }
}

// Allows localization
// Hence, we can have numbers like a+bsqrt(2)/sqrt(2)^k
// See code in local_ring.rs
// And the comments in mod.rs
impl Localizable for Zroot2
{
    fn is_divisible(self) -> bool{
        if self.0%2==0 { true }
        else { false }
    }
    fn reduce_by_dividing(mut self) -> u32
    {
        let trail0 = self.0.trailing_zeros();
        let trail1 = self.1.trailing_zeros();

        if trail0<=trail1
        {
            // a*2^t0 + b*2^(t1+1/2) =  2^t0( a +b*2^(t1-t0+1/2))
            // if t0=t1 then
            // a*2^t0 + b*2^(t1+1/2) =  2^(t0)( a +b*sqrt(2))
            self.0 << (trail0);
            self.1 << (trail0);
            return 2*trail0;
        }
        else 
        {
            // a*2^t0 + b*2^(t1+1/2) =  2^(t1+1/2)( a*2^(t0-t1-1+1/2) +b)
            self.0 >> (trail1+1);
            self.1 >> (trail1);
            (self.0,self.1) = (self.1,self.0);
            return 2*trail1+1;
        }
    }


    fn perform_n_multiplications(self, n: u32) -> ()
    {
        // (  a+bsqrt(2) )*2^(n/2) = (a*2^(n/2) +b*( k/2 + 1/2) ) 
        let ntemp=n >> 1;
        self.0 << ntemp;
        self.1 << ntemp;
        
        if n%2==1
        {
            (self.0,self.1) = (self.1 << 1,self.0);
        }
    }
}

// Get zero and one as ring elements
// impl<T> Constructs<T> for Zroot2
// {
//     fn zero() -> Self {
//         return Self(0,0);

//     }

//     fn one() -> Self {
//         return Self(1,0);
//     }
// }
