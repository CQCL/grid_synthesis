use crate::structs::rings::Localizable; //Localizing trait
use crate::structs::rings::Conj; //Conjugation trait


// We bring them in so that we can overload the operators
// Rust must learn how to do arithmetics in our rings
use std::ops::Neg; 
use std::ops::Add; 
use std::ops::Sub; 
use std::ops::Mul; 
use std::ops::Div; 
use std::ops::Rem; 
use std::cmp::PartialEq; 

// For display
use std::fmt::Result;
use std::fmt::Display;
use std::fmt::Formatter;


// Setting up my boundaries
use crate::structs::rings::Int;
use crate::structs::rings::LogDepInt;

// Num traits
use num_traits::Num;
use num_traits::Zero;
use num_traits::One;
use num_traits::NumOps;
use num_traits::NumCast;
use num_traits::ToPrimitive;


// Quadratic number field with root 2
// This struct assumes that you are going to localize it 
// at sqrt(2)
#[derive(Copy,Debug,Clone,PartialEq,PartialOrd)]
pub struct Zroot2(pub Int,pub Int); //a+b\sqrt(2)

// Rust must know how to diplay elements of this ring
// Rust could learn some latex
impl Display for Zroot2{
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        if self.0!=0
        {
            if self.1<0 
            {
                write!(f, "{}{}√2", self.0, self.1)
            }
            else if self.1>0
            {
                if self.1==1
                {
                    write!(f, "{}+√2", self.0)
                }
                else
                {
                    write!(f, "{}+{}√2", self.0, self.1)
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
                write!(f, "{}√2", self.1)
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
    fn reduce_by_dividing(mut self) -> ( Self, LogDepInt )
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


    fn perform_n_multiplications(mut self, n: LogDepInt) -> Self 
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
impl Zroot2
{
    pub fn norm(self) -> Int
    {
        return self.0*self.0-2*self.1*self.1;
    }

    // Should check if norm is a unit 
    fn is_unit(self) -> bool
    {
        if self.norm()==1 || self.norm()==-1
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}

impl From<Int> for Zroot2 {
    fn from(int: Int) -> Self {
        Zroot2(int.try_into().unwrap(),0)
    }
}




impl Zero for Zroot2
{
    fn zero() -> Self
    {
        return Zroot2(0,0);
    }

    fn is_zero(&self) -> bool
    {
        return self.0==0 && self.1==0;
    }

}

impl One for Zroot2
{
    fn one() -> Self
    {
        return Zroot2(1,0);
    }


}


// Negatation on Zroot2
impl Neg for Zroot2
{
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
        // Naive multiplication, for now
        Zroot2(
            other.0*self.0 + 2*other.1*self.1,
            other.0*self.1 + other.1*self.0 
            )
    }
}

// Conjugate Zroot2 elements
impl Conj for Zroot2 {
    fn conj(&self) -> Self 
    {
        Zroot2(self.0,-self.1)
    }
}


impl Rem for Zroot2 
{
    type Output = Self;
    fn rem(self, other: Self) -> Self
    {
        panic!("TODO");
    }

}

impl Div for Zroot2 
{
    type Output = Self;
    fn div(self, other: Self) -> Self
    {
        // Zroot2 is norm-Euclidean
        // Therefore a division algorithm is very much possible for Zroot2
        // Will look in to this: TODO
        let mut nor = other.norm();
        if nor==1 || nor==-1
        {
            let norself = Self{
                0: nor,
                1: 0
            };
            return self*other.conj()*norself;
        }
        else if nor!=0
        {
            // WARNING: This is bad mathematics is bad mathematics 
            //          because when self/other is not exactly in Zroot2
            //          the output of this will be self = q*other + r
            //          where absolute value of r.norm() is less than
            //          absolute value of other.norm().
            //
            //          Basically, this does Euclidean division
            //          and finds q given above
            let numerator = self*other.conj();
            let denominator = nor;


            return Zroot2( nearest_integer(numerator.0,nor),  nearest_integer(numerator.1, nor) );

            todo!();

        }
        else 
        {
            println!("Wanted to divide {} with {}",self,other);
            panic!("Division impossible");
        }

    }
}


// Return closest integer to top/bottom
// Will send halfs to the ceiling
// Probably shouldn't matter
// TODO: Optimize this function
pub fn nearest_integer(top :Int, bottom: Int) -> Int
{
    if bottom == 0
    {
        panic!("What do you think?");
    }
    else if top ==0
    {
        return 0;
    }
    else if ( top > 0 && bottom > 0 )     
    { 
        let twice = (top << 1)/bottom;
        if twice%2==0
        {
            return twice>>1;
        }
        else 
        {
            return (twice+1)>>1  ;
        }

    }
    else if ( top < 0 && bottom < 0 ) 
    {
        return nearest_integer(-top, -bottom);
    }
    else if ( top > 0 && bottom < 0 ) 
    {
        return -nearest_integer(top, -bottom);
    }
    else 
    {
        return -nearest_integer(-top, bottom);
    }
}

impl Num for Zroot2
{
    // Some rust requirements to make a Complex compatible local
    type FromStrRadixErr = std::num::ParseIntError  ;
    fn from_str_radix(_: &str, _: u32) -> core::result::Result<Self, Self::FromStrRadixErr>
    {
        // Meant to parse a string and take the input as what is in the string
        // Kinda complicated to do for this local ring
        // Will skip it
        panic!("Unimplimented or impossible");
        //
    }
}

impl ToPrimitive for Zroot2
{
    // Rust needs this. 
    // The NumCast library does not work otherwise
    // Honestly, I cannot convert Zroot2 to i64 or u64
    fn to_i64(&self) -> Option<i64> 
    { 
        return None;
    }

    fn to_u64(&self) -> Option<u64> 
    { 
        return None;
    }

}

impl NumCast for Zroot2
{
    fn from<T>(given: T) -> Option<Self>
        where T: ToPrimitive
        {
            return Some(Self
                        {
                            0: given.to_i64().unwrap(),
                            1: 0
                        });
        }
}
