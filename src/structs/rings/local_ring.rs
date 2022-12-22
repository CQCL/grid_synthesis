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
use std::ops::Div; 
use std::ops::Rem; 
use std::cmp::PartialEq; 

// For display
use std::fmt::Result;
use std::fmt::Display;
use std::fmt::Formatter;

// Num traits
use num_traits::Num;
use num_traits::Zero;
use num_traits::One;
use num_traits::ToPrimitive;
use num_traits::NumCast;


// // Localization trait. 
// See mod.rs for some info
// Or Zroot2
use crate::structs::rings::Localizable; 
use crate::structs::rings::Int; 
use crate::structs::rings::Conj; //Conjugation trait
use crate::structs::rings::LocalizableNorm; //Norm trait

//Integer type is set globally
// use crate::structs::rings::Int; 
use crate::structs::rings::LogDepInt; 

// Ring of numbers of the form a/2^n for integers a and n
// These are dyadic integers
// Each integer looks something like
#[derive(Debug, Copy, Clone, PartialEq, PartialOrd, Hash, Eq)]
pub struct Local<T> 
{
    pub num: T,
    // The integer a above
    pub log_den: LogDepInt,

    // Perhaps I could save an bool
    // That remembers weather or not a+bsqrt(2)/sqrt(2)^k i
    // is in the simplest form
    //
    // pub is_reduced: bool
}

// Internal function
// Used to make the numerator independent of denominator
impl<T> Local<T> 
where T: Localizable+PartialEq+Copy,
      T: Num

{
    fn fix(&mut self) -> Self 
    {
        // println!("Fixing height {}",self.log_den);
        if self.num != T::zero()
        { 
            // println!("Number is non-zero");
            if self.num.is_divisible()
            {
                // println!("Number is divisible");
                let pow: LogDepInt;
                (self.num, pow) = self.num.reduce_by_dividing();
                // println!("Self obtained after it performed {} divisions",pow);
                // println!("log_den is {}",self.log_den);
                self.log_den = self.log_den - pow;
                //
            }
            // println!("log_den is {}",self.log_den);

            // Return ownership
            return *self;
        }
        else
        {
            // println!("This should be minus infinity");
            self.log_den=0;

            // Return ownership
            return *self;
        }
    }

    fn logden(self) -> LogDepInt
    {
        return self.log_den;
    }

    fn local_gen() -> Self
    {
        return Self
        {
            num: T::one(),
            log_den: -1
        }
    }
}


impl<T> Neg for Local<T>
where T:Neg<Output=T>
{
    type Output = Self;
    fn neg(self) -> Self {
        Self{
            num: -self.num, 
            log_den: self.log_den,
        }
    }
}


// Teaching rust how to add Local<T> elements
impl<T> Add for Local<T> 
where T: Localizable+PartialEq+Copy,
      T: Num
{
    type Output = Self;
    fn add(self, other: Self) -> Self {
        if other.num==T::zero()
        {
            return self;
        }
        if self.num==T::zero()
        {
            return other;
        }
        // println!("Want to add height {} and height {}",self.log_den,other.log_den);
        if self.log_den>other.log_den 
        {
            let mut other_temp = other.num;
            // println!("Will invoke perform_n_multiplications with n={}",self.log_den-other.log_den);
            other_temp = other_temp.perform_n_multiplications(self.log_den-other.log_den);

            let temp = 
                Self{
                    num: self.num + other_temp,
                    log_den: self.log_den
                };

            return temp;
        }
        else
        {
            let mut self_temp = self.num;
            self_temp = self_temp.perform_n_multiplications(other.log_den-self.log_den);
            let mut temp = Self
            {
                num: other.num + self_temp,
                log_den: other.log_den
            };
            // Possible speedup by returning without fixing??
            // p-adic valuation of unequal
            if other.log_den == self.log_den
            {
                // println!("Fix routine");
                temp= temp.fix();
                // println!("Fixed");
            }
            return temp;
        }

    }
}


// Teaching rust how to multiply Local elements
impl<T> Mul for Local<T> 
where T: Mul<Output=T>+PartialEq+Copy,
      T: Localizable,
      T: Num
{
    type Output = Self;

    fn mul(self, other: Self) -> Self 
    {

        if self.num == T::zero() || other.num == T::zero()
        {
            return Self::zero();
        }
        else
        {
            let mut temp = Self
            {
                num: self.num*other.num,
                log_den: self.log_den+other.log_den
            };
            // return temp;
            temp.fix();
            return temp;
        }

    }
}



// Teaching rust how to subtract Dyad elements
impl<T> Sub for Local<T> 
where T: Num,
      T: Localizable,
      T: Copy
{

    type Output = Self;
    fn sub(self, other: Self) -> Self 
    {
        let other_neg = Self
        {
            num: T::zero() - other.num,
            log_den: other.log_den,
        };
        return self+other_neg;
    }

}


// Nicely display Local<T> elements
impl<T> Display for Local<T>
where T: Display+PartialEq,
      T: Zero+One
{
    fn fmt(&self, f: &mut Formatter) -> Result{
        if self.log_den>0  
        {   
            if self.num== T::one()
            {
                write!(f,"√2^(-{})", self.log_den)
            }
            else if self.num == T::zero()
            {
                write!(f,"0")
            }
            else 
            {
                write!(f,"({})*√2^(-{})",self.num, self.log_den)
            }
        }        
        else if self.log_den<0  
        { 
            if self.num== T::zero()
            {
                write!(f,"√2^({})",-self.log_den)
            }
            else if self.num == T::zero()
            {
                write!(f,"0")
            }
            else 
            {
                write!(f,"({})*√2^({})",self.num, -self.log_den)
            }
        }
        else 
        {   
            write!(f,"({})",self.num) 
        }
    }
}





impl<T> ToPrimitive for Local<T>
{
    // Rust needs this. 
    // The NumCast library does not work otherwise
    // Honestly, I cannot convert Local<T> to i64 or u64
    fn to_i64(&self) -> Option<i64> 
    { 
        return None;
    }

    fn to_u64(&self) -> Option<u64> 
    { 
        return None;
    }

}

impl<T> NumCast for Local<T>
where T: NumCast,
      T: Localizable,
      T: PartialEq,
      T: Copy,
      T: Num
{
    fn from<E>(given: E) -> Option<Self>
        where E: ToPrimitive
        {
            let mut temp = Self
            {
                num: T::from(given).unwrap(),
                log_den:0
            };
            temp = temp.fix();
            return Some(temp);
        }
}



// To construct Local<T> directly from integers
// impl<T> From<Int> for Local<T> 
// where T: Num,
//       T: Localizable+PartialEq,
//       T: Copy
// {
//     fn from(input : Int) -> Self 
//     {
//         panic!("Need to deprecate");
//         // let out = Self{
//         //     num: T::from(input.try_into().unwrap()),
//         //     log_den: 0
//         // };

//         // if input!=1 && input!=0
//         // {
//         //     out.fix();
//         // }

//         // return out;
//     }
// }


// To construct Local<T> from T
impl<T> Local<T> 
where T: Localizable,
      T: PartialEq+Copy,
      T: Num
{
    pub fn from_base(input : T) -> Self 
    {

        let mut out = Self{
            num: input,
            log_den: 0
        };

        if input!=T::one() && input!=T::zero()
        {
            out.fix();
        }

        return out;
    }
}


// Conjugate Complex elements
impl<T> Conj for Local<T> 
where T: Conj+Copy
{
    fn conj(&self) -> Self {
        let temp_num = self.num;
        temp_num.conj();
        Self{
            num: temp_num,
            log_den: self.log_den
        }
    }
}


// Teaching rust how to multiply Local elements
impl<T> Div for Local<T> 
where T: Num,
      T: PartialEq+Copy,
      T: Localizable,
{
    type Output = Self;

    fn div(self, other: Self) -> Self 
    {
        if self.num==T::zero()
        {
            return Self::zero();
        }
        else if other.num.is_zero()
        {
            panic!("Division by zero")
        }
        // println!("Dividing height {} with height {}",self.log_den,other.log_den);
        else 
        {
            let temp = Self
            {
                num: self.num/other.num,
                log_den: self.log_den-other.log_den
            };

            //fixing might not be needed
            // temp.fix(); 
            return temp;
        }
    }
}

impl<T> One for Local<T>
where T: Copy+PartialEq,
      T: Num,
      T: Localizable,
{
    fn one() -> Self
    {
        return Self
        {
            num: T::one(),
            log_den: 0,
        }
    }


}

impl<T> Zero for Local<T>
where T: Copy+PartialEq,
      T: Num,
      T: Localizable,
{
    fn zero() -> Self
    {
        return Self
        {
            num: T::zero(),
            log_den: 0,
        }
    }

    fn is_zero(&self) -> bool
    {
        return self.num == T::zero();
    }

}

impl<T> Rem for Local<T>
{

    type Output = Self;
    fn rem(self, _ : Self) -> Self
    {
        panic!("Unimplimented or impossible");
        // return self;
    }
}

// impl<T> FromStrRadixErr for Local<T>
// {

// }

impl<T> Num for Local<T>
where T: Num,
      T: Copy,
      T: Localizable+Conj,

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



impl<T> Local<T>
where T: LocalizableNorm
{
    pub fn norm(&self) -> Local::<Int>
    {
        return Local::<Int>
        {
            num: self.num.norm(),
            log_den: 2*self.log_den
        }
    }

}


