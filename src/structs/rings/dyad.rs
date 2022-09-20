// We bring them in so that we can overload the operators
// Rust must learn how to do arithmetics in our rings
use std::ops::Neg; 
use std::ops::Add; 
use std::ops::Mul; 

// use std::fmt; // To teach rust how to display our ring elements


// Ring of numbers of the form a/2^n for integers a and n
// These are dyadic integers
// Each integer looks something like
pub struct Dyad {
    pub num: i64, // The integer a above
    pub log_den: i64, // The integer n above 
                      // log_den could be made i8 or something, since 2^(2^8) >2^64
}

// Internal function
// Used to make the numerator odd
impl Dyad {
    fn fix(self) -> Dyad {
        while self.num%2==0
        {
            self.num=self.num/2;
            log_den=log_den+1;
        }
    }
}

impl Neg for Dyad{
    type Output = Dyad;

    fn neg(self) -> Dyad {
        Dyad{
            num: -self.num, log_den: self.log_den
        }
    }

}


// Teaching rust how to add Dyad elements
impl Add for Dyad {
    type Output = Dyad;

    fn add(self, other: Dyad) -> Dyad {
        if self.log_den>other.log_den 
        {
            temp = Dyad{
                num: self.num + other.num*2^(self.log_den-other.log_den),
                log_den = self.log_den
            };

            return temp.fix()
        }
        else
        {
            temp = Dyad{
                num: other.num + self.num*2^(other.log_den-self.log_den),
                log_den = other.log_den
            };

            return temp.fix()
        }
    }
}


