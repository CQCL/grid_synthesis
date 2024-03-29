// Consider the eighth cyclotomic field
// and consider the ring of integers in this field
// THis is what we have here
//
//


use crate::structs::rings::Conj; //Conjugation trait
use crate::structs::rings::Localizable;
use crate::structs::rings::Int; //Conjugation trait
use crate::structs::rings::zroot2::Zroot2;
use crate::structs::rings::local_ring::Local;

use crate::structs::rings::special_values::onebyroot2loc;

use crate::algorithms::near_int::nearest_integer;


use num_traits::NumCast;
use num_traits::ToPrimitive;
use num_traits::Zero;
use num_traits::One;

// We bring them in so that we can overload the operators
// Rust must learn how to do arithmetics in our rings
use std::ops::Neg; 
use std::ops::Add; 
use std::ops::Sub; 
use std::ops::Mul; 
use std::ops::Div; 
use std::ops::Rem; 

// For display
use std::fmt::Result;
use std::fmt::Display;
use std::fmt::Formatter;


// Ring of numbers of the form a+b\omega+c\omega^2+d\omega^3
// where \omega is the eighth root of unity
#[derive(Debug,Copy, Clone)]
pub struct Zomega(pub Int, pub Int, pub Int, pub Int); 



// Conjugate Zomega elements
impl Conj for Zomega {
    fn conj(&self) -> Zomega {
        return Zomega(self.0,-self.3, -self.2, -self.1);
    }
}

// Conjugate sqrt2 in Zomega elements 
impl Zomega 
{
    fn conj_rt2(&self) -> Zomega 
    {
        return Zomega(self.0,-self.1, self.2, -self.3);
    }
}

// Nicely display Zomega Matrices
impl Display for Zomega
{
    fn fmt(&self, f: &mut Formatter) -> Result{
        write!(f,"{}+{}w+{}w^2+{}w^3",self.0, self.1, self.2,self.3)
    }
}

// Negatation on Zomega
impl Neg for Zomega{
    type Output = Zomega;

    fn neg(self) -> Zomega {
        Zomega(-self.0,-self.1,-self.2,-self.3)
    }
}


// Teaching rust how to add Zomega elements
impl Add for Zomega {
    type Output = Zomega;

    fn add(self, other: Zomega) -> Zomega {
        Zomega(self.0+other.0,self.1+other.1,self.2+other.2,self.3+other.3)
    }
}

// Teaching rust how to subtract Zomega elements
impl Sub for Zomega {
    type Output = Zomega;

    fn sub(self, other: Zomega) -> Zomega {
        self+(-other) //subtraction is just adding the additive inverse
    }
}


// Teaching rust how to multiply Zomega elements
impl Mul for Zomega {
    type Output = Zomega;

    fn mul(self, other: Zomega) -> Zomega {
        Zomega(
            other.0*self.0 - other.1*self.3 - other.2*self.2 - other.3*self.1,
            other.0*self.1 + other.1*self.0 - other.2*self.3 - other.3*self.2,
            other.0*self.2 + other.1*self.1 + other.2*self.0 - other.3*self.3,
            other.0*self.3 + other.1*self.2 + other.2*self.1 + other.3*self.0
            )
    }
}

// impl Localizable for Zomega{

//         todo!();
// }

// Get zero and one as Zomega numbers
// impl<T> Constructs<T> for Zomega
// {
//     fn one() -> Self {
//         let a0=Constructs::<Dyad>::zero();
//         let a1=Constructs::<Dyad>::one();
//         return Zomega(a1,a0,a0,a0);
//     }
    
//     fn zero() -> Self {
//         let a0=Constructs::<Dyad>::zero();
//         return Zomega(a0,a0,a0,a0);
//     }
// }


// Teaching rust how to compare these ring elements
impl PartialEq for Zomega
{
    fn eq(&self, other: &Self) -> bool {
        // println!("{},{},{},{}",self.0,self.1,self.2,self.3);
        // println!("{},{},{},{}",other.0,other.1,other.2,other.3);
        return self.0==other.0 && self.1==other.1 && self.2==other.2 && self.3==other.3;
    }
}



impl ToPrimitive for Zomega
{
    // Rust needs this. 
    // The NumCast library does not work otherwise
    // Honestly, I cannot convert Zomega to i64 or u64
    fn to_i64(&self) -> Option<i64> 
    { 
        return None;
    }

    fn to_u64(&self) -> Option<u64> 
    { 
        return None;
    }

}

impl NumCast for Zomega
{
    fn from<T>(given: T) -> Option<Self>
        where T: ToPrimitive
        {
            return Some(Self
                        {
                            0: given.to_i128().unwrap(),
                            1: 0,
                            2: 0,
                            3: 0
                        });
        }
}


impl Zero for Zomega
{
    fn zero() -> Self
    {
        return Zomega(0,0,0,0);
    }

    fn is_zero(&self) -> bool
    {
        return self.1==0 && self.0 == 0 && self.2==0 && self.3==0;
    }
}

impl One for Zomega
{
    fn one() -> Self
    {
        return Zomega(1,0,0,0);
    }
}

impl Zomega
{
    pub fn norm(&self) -> Int
    {
        let term1 = (self.0*self.0 + self.1*self.1 + self.2*self.2 + self.3*self.3);
        let term2 = (self.0*self.1 + self.1*self.2 - self.0*self.3 + self.2*self.3);
        return term1*term1 - 2*term2*term2;
    }

}


// This code here has been the culprit for a lot of integer overflows
// This is basically because calculating norm is like calculating the 
// fourth powers of some integers; 
impl Div for Zomega 
{

    type Output = Self;
    fn div(self, other: Self) -> Self
    {

        // Zomega is norm-Euclidean
        let nor = other.norm();
        if nor==1 
        {
            return self*other.conj()*other.conj_rt2()*other.conj_rt2().conj();
        }
        else if nor == -1
        {
            return -self*other.conj()*other.conj_rt2()*other.conj_rt2().conj();
        }
        else if nor!=0
        {
            // WARNING: This is bad mathematics 
            //          because when self/other is not exactly in Zroot2
            //          the output of this will be self = q*other + r
            //          where absolute value of r.norm() is less than
            //          absolute value of other.norm().
            //
            //          Basically, this does Euclidean division
            //          and finds q given above
            let numerator =  self*other.conj()*other.conj_rt2()*other.conj_rt2().conj();
            let denominator = nor;


            return Zomega( nearest_integer(numerator.0,nor),  nearest_integer(numerator.1, nor) , nearest_integer( numerator.2, nor) , nearest_integer( numerator.3, nor) );

        }
        else 
        {
            println!("Wanted to divide {} with {}",self,other);
            panic!("Division impossible");
        }

    }
}



impl Rem for Zomega 
{
    type Output = Self;
    fn rem(self, other: Self) -> Self
    {
        let q = self/other;
        return self-q*other;
    }

}


type Loc = Local::<Zroot2>;
impl Zomega
{
    // Spit out real and imaginary parts as Local::Zroot2 elements
    // a + b*omega + c*omega^2 + d*omega^3 
    // =  ( b-d + a*sqrt2 )/ sqrt2 + i (b+d + c * sqrt2) / sqrt2
    pub fn real_part(&self) -> Loc
    {
        let zrt = Zroot2( self.1 - self.3 , self.0 );
        let term = Loc::from_base(zrt);
        return term*onebyroot2loc();
    }
    pub fn imag_part(&self) -> Loc
    {
        let zrt = Zroot2( self.1 + self.3 , self.2 );
        let term = Loc::from_base(zrt);
        return term*onebyroot2loc();
    }

    pub fn from_zroot2(input: Zroot2) -> Zomega 
    {
        return Zomega(input.0, input.1, 0 ,-input.1);
    }

}
