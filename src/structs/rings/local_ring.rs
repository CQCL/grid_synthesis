// The main purpose of this is to define 
// localization with respect to a prime ideal in a number field
//
//
// For non-mathematicians, what I want to really use it for is
// To define numbers of the form (a+sqrt(2)b)/sqrt(2)^k for integers a,b,k
// However, I will do a very general version which involves no 
// description of sqrt(2)


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

use crate::structs::rings::Constructs; // Construction trait

// Ring of numbers of the form a/2^n for integers a and n
// These are dyadic integers
// Each integer looks something like
#[derive(Debug, Copy, Clone)]
pub struct Local<T> 
{
    pub num: T,
    // The integer a above
    pub log_den: u32,
    
    // Perhaps I could save an bool
    // That remembers weather or not a+bsqrt(2)/sqrt(2)^k i
    // is in the simplest form
    //
    // pub is_reduced: bool
}

// Internal function
// Used to make the numerator independent of denominator
impl Local<T> 
where T: Localizable
{
    fn fix(mut self) -> T {

        if self.num != T::zero(){ 

            while self.is_divisible
            {
                self.perform_one_division();
                self.log_den= self.log_den-1;
            }


        }
        else
        {
            self.log_den=0;
        }
    }
}

impl Neg for Local<T>{
    // type Output = T;
    fn neg(self) -> Self {
        Self{
            num: -self.num, 
            log_den: self.log_den,
        }
    }

}


// Teaching rust how to add Local<T> elements
impl Add for Local<T> 
where 
{
    type Output = Self;

    fn add(self, other: Self) -> Self {
        if self.log_den>other.log_den 
        {
            let temp = Dyad{
                num: self.num + other.num*(2_i64.pow(self.log_den-other.log_den)),
                log_den: self.log_den
            };

            return temp.fix()
        }
        else
        {
            let temp = Dyad{
                num: other.num + self.num*( 2_i64.pow(other.log_den-self.log_den) ),
                log_den: other.log_den
            };

            // println!("Panicking with");
            // println!("{}-{}",self.log_den,other.log_den);
            return temp.fix()
        }
    }
}


// Teaching rust how to multiply Dyad elements
impl Mul for Dyad {
    type Output = Dyad;

    fn mul(self, other: Dyad) -> Dyad {
        let temp = Dyad{
            num: self.num*other.num,
            log_den: self.log_den+other.log_den
        };
        return temp.fix()
    }
}


// Teaching rust how to subtract Dyad elements
impl Sub for Dyad {
    type Output = Dyad;

    fn sub(self, other: Dyad) -> Dyad {
        self+(-other) //subtraction is just adding the additive inverse
    }

}



// Nicely display Local<T> elements
impl Display for Dyad
{
    fn fmt(&self, f: &mut Formatter) -> Result{
        write!(f,"{}*2^(-{})",self.num, self.log_den)
    }
}


// Teaching rust how to compare these ring elements
impl PartialEq for Dyad {
    fn eq(&self, other: &Self) -> bool {
        if self.num%2==0 && self.num!=0
        {
            self.fix();
        }
        if other.num%2==0 && other.num!=0
        {
            other.fix();
        }
        return self.num==other.num && self.log_den==other.log_den;
    }
}


// Get zero and one as Complex numbers
impl<T> Constructs<T> for Dyad
{
    fn zero() -> Self {
        return Dyad{
            num: 0,
            log_den: 0
        };
    }
    
    fn one() -> Self {
        return Dyad{
            num: 1,
            log_den: 0
        };
    }
}

// To construct Local<T> directly from integers
impl From<u32> for Dyad {
    fn from(int: u32) -> Self {
        Dyad{
            num: int as i64,
            log_den: 0
        }
    }
}
