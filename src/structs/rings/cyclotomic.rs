// Trying to implement a localizable cyclotomic ring
// UNFINISHED


use crate::structs::rings::Conj; //Conjugation trait
use crate::structs::rings::Localizable; //Localizing trait


// We bring them in so that we can overload the operators
// Rust must learn how to do arithmetics in our rings
use std::ops::Neg; 
use std::ops::Add; 
use std::ops::Sub; 
use std::ops::Mul; 

// For display
use std::fmt::Result;
use std::fmt::Display;
use std::fmt::Formatter;

// Globally defined integers
use crate::structs::rings::Int;

// Ring of numbers of the form (a+b\omega+c\omega^2+d\omega^3)
// where \omega is the eighth root of unity
// There will be a localization implementation below
#[derive(Debug, Copy, Clone)]
pub struct Cyclotomic(Int, Int, Int, Int);
// variable.0, variable.1 and so on are intgers


// Conjugate Cyclotomic elements
impl<T> Conj<T> for Cyclotomic {
    fn conj(self) -> Cyclotomic{
        Cyclotomic{ 
            0: self.0,
            1: -self.3,
            2: -self.2,
            3: -self.1
        }
    }
}

// Nicely display Cyclotomic elements
impl Display for Cyclotomic
{
    fn fmt(&self, f: &mut Formatter) -> Result{
        write!(f,"({}+{}w+{}w^2+{}w^3)",self.0, self.1, self.2,self.3)
    }
}

// Negatation on Cyclotomic
impl Neg for Cyclotomic{
    type Output = Cyclotomic;

    fn neg(self) -> Cyclotomic {
        Cyclotomic{ 
            0: -self.0,
            1: -self.1,
            2: -self.2,
            3: -self.3,
        }
    }
}


// Teaching rust how to add Cyclotomic elements
impl Add for Cyclotomic {
    type Output = Cyclotomic;

    fn add(self, other: Cyclotomic) -> Cyclotomic {
        // WRONG CODE!!!! 
        // NEED TO WORK ON THIS
        Cyclotomic{ 
        0: other.0+self.0,
        1: other.1+self.1,
        2: other.2+self.2,
        3: other.3+self.3,
        }
    }
}

// Teaching rust how to subtract Cyclotomic elements
impl Sub for Cyclotomic {
    type Output = Cyclotomic;

    fn sub(self, other: Cyclotomic) -> Cyclotomic {
        self+(-other) //subtraction is just adding the additive inverse
    }
}


// Teaching rust how to multiply Cyclotomic elements
impl Mul for Cyclotomic {
    type Output = Cyclotomic;

    fn mul(self, other: Cyclotomic) -> Cyclotomic {
        Cyclotomic(
            other.0*self.0 - other.1*self.3 - other.2*self.2 - other.3*self.1,
            other.0*self.1 + other.1*self.0 - other.2*self.3 - other.3*self.2,
            other.0*self.2 + other.1*self.1 + other.2*self.0 - other.3*self.3,
            other.0*self.3 + other.1*self.2 + other.2*self.1 + other.3*self.0
            )
    }
}



// Teaching rust how to compare these ring elements
impl PartialEq for Cyclotomic
{
    fn eq(&self, other: &Self) -> bool {
        // println!("{},{},{},{}",self.0,self.1,self.2,self.3);
        // println!("{},{},{},{}",other.0,other.1,other.2,other.3);
        return self.0==other.0 && self.1==other.1 && self.2==other.2 && self.3==other.3;
    }
}



// To construct Cyclotomic directly from integers
impl From<Int> for Cyclotomic {
    fn from(int: Int) -> Self {
        Cyclotomic(int,0,0,0)
    }
}



// Allows localization
// Hence, we can have numbers like a+bsqrt(2)/sqrt(2)^k
// See code in local_ring.rs
// And the comments in mod.rs
impl Localizable for Cyclotomic
{
    fn is_divisible(self) -> bool{
        if (self.0-self.1+self.2-self.3)%2==0 { true }
        else { false }
    }
    fn reduce_by_dividing(mut self) -> Int
    {
        // TODO
        panic!("Not yet implemented");
        // return 0;
    }


    fn perform_n_multiplications(mut self, n: Int) -> ()
    {
        // TODO
        panic!("Not yet implemented");
    }
}

