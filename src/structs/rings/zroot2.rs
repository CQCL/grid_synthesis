use crate::structs::rings::Localizable; //Localizing trait
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


// Setting up my boundaries
use crate::structs::rings::Int;



// Quadratic number field with root 2
// This struct assumes that you are going to localize it 
// at sqrt(2)
#[derive(Copy,Debug,Clone)]
pub struct Zroot2(pub Int,pub Int); //a+b\sqrt(2)

// Rust must know how to diplay elements of this ring
// Rust could learn some latex
impl Display for Zroot2{
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        if self.0!=0
        {
            if self.1<0 
            {
                write!(f, "{}{}\\sqrt{{2}}", self.0, self.1)
            }
            else if self.1>0
            {
                if self.1==1
                {
                    write!(f, "{}+\\sqrt{{2}}", self.0)
                }
                else
                {
                    write!(f, "{}+{}\\sqrt{{2}}", self.0, self.1)
                }
            }
            else 
            {
                write!(f, "{}",self.0)
            }
        }
        else
        {
            if self.1!=0
            {
                write!(f, "{}\\sqrt{{2}}", self.1)
            }
            else
            {
                write!(f, "{}", self.1)
            }
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
    fn reduce_by_dividing(mut self) -> ( Self, Int )
    {
        // println!("Hi, I'm reducing self");

        // Remember, trailing zeroes can be size(Int)
        // when Int==0
        
        let trail0 = self.0.trailing_zeros();
        let trail1 = self.1.trailing_zeros();

        // println!("Here are trailing zeroes: {} {}", trail0,trail1);
        // println!("      of our variables  : {} {}", self.0,self.1);

        if trail0<=trail1
        {
            // a*2^t0 + b*2^(t1+1/2) =  2^t0( a +b*2^(t1-t0+1/2))
            // if t0=t1 then
            // a*2^t0 + b*2^(t1+1/2) =  2^(t0)( a +b*sqrt(2))
            self.0 = self.0  >> (trail0);
            self.1 = self.1  >> (trail0);
            // println!("We modified them to be  : {} {}", self.0,self.1);
            // println!("Returning self as {}",self);
            return (self, (2*trail0).try_into().unwrap());
        }
        else 
        {
            // a*2^t0 + b*2^(t1+1/2) =  2^(t1+1/2)( a*2^(t0-t1-1+1/2) +b)
            self.0 = self.0 >> (trail1+1);
            self.1 = self.1 >> (trail1);
            (self.0,self.1) = (self.1,self.0);
            // println!("We modified them to be  : {} {}", self.0,self.1);
            // println!("Returning self as {}",self);
            return (self, ( 2*trail1+1 ).try_into().unwrap());
        }
    }


    fn perform_n_multiplications(mut self, n: Int) -> Self 
    {
        // (  a+bsqrt(2) )*2^(n/2) = (a*2^(n/2) +b*( k/2 + 1/2) ) 
        // println!("Will multiply {} times",n);
        let ntemp=n >> 1;
        self.0 = self.0 << ntemp;
        self.1 = self.1 << ntemp;

        if n%2==1
        {
            (self.0,self.1) = (self.1 << 1,self.0);
        }

        return self;
    }
}

impl From<Int> for Zroot2 {
    fn from(int: Int) -> Self {
        Zroot2(int.try_into().unwrap(),0)
    }
}



impl PartialEq for Zroot2
{
    fn eq(&self, other: &Self) -> bool {
        return self.0==other.0 && self.1==other.1;
    }
}



// Negatation on Zroot2
impl Neg for Zroot2{
    type Output = Zroot2;
    fn neg(self) -> Zroot2 {
        Zroot2(-self.0,-self.1)
    }
}

// Teaching rust how to subtract Zroot2 elements
impl Sub for Zroot2 {
    type Output = Zroot2;

    fn sub(self, other: Zroot2) -> Zroot2 {
        self+(-other) //subtraction is just adding the additive inverse
    }
}



// Teaching rust how to multiply Zroot2 elements
impl Mul for Zroot2 {
    type Output = Zroot2;

    fn mul(self, other: Zroot2) -> Zroot2 {
        Zroot2(
            other.0*self.0 + 2*other.1*self.1,
            other.0*self.1 + other.1*self.0 
            )
    }
}

// Conjugate Zroot2 elements
impl<T> Conj<T> for Zroot2 {
    fn conj(self) -> Zroot2 {
        Zroot2(self.0,-self.1)
    }
}
