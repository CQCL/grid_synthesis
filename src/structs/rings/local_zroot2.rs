// The main purpose of this is to define 
// localization with respect to a prime ideal in a number field
//
//
// For non-mathematicians, what I want to really use it for is
// To define numbers of the form (a+sqrt(2)b)/sqrt(2)^k for integers a,b,k

// This is now being written directly without using dependent types

// We bring them in so that we can overload the operators
// Rust must learn how to do arithmetics in our rings
use std::ops::Neg; 
use std::ops::Add; 
use std::ops::Sub; 
use std::ops::Mul; 
use std::ops::Div; 
use std::cmp::PartialEq; 

// For display
use std::fmt::Result;
use std::fmt::Display;
use std::fmt::Formatter;

// // Localization trait. 
// See mod.rs for some info
// Or Zroot2
use crate::structs::rings::Localizable; 
use crate::structs::rings::Fixable; 
use crate::structs::rings::Conj; //Conjugation trait

//Integer type is set globally
use crate::structs::rings::Int; 

// Ring of numbers of the form a/2^n for integers a and n
// These are dyadic integers
// Each integer looks something like
#[derive(Debug, Copy, Clone)]
pub struct Local<T> 
{
    pub num: T,
    // The integer a above
    pub log_den: Int,

    // Perhaps I could save an bool
    // That remembers weather or not a+bsqrt(2)/sqrt(2)^k i
    // is in the simplest form
    //
    // pub is_reduced: bool
}

// Internal function
// Used to make the numerator independent of denominator
impl<T> Fixable for Local<T> 
where T: Localizable+PartialEq+From<Int>+Copy
{
    fn fix(mut self) -> Self 
    {
        // println!("Fixing height {}",self.log_den);
        if self.num != T::from(0)
        { 
            // println!("Number is non-zero");
            if self.num.is_divisible()
            {
                // println!("Number is divisible");
                let pow: Int;
                (self.num, pow) = self.num.reduce_by_dividing();
                // println!("Self obtained after it performed {} divisions",pow);
                // println!("log_den is {}",self.log_den);
                self.log_den = self.log_den - pow;
                // println!("log_den is {}",self.log_den);
                //
            }
            // Return ownership
            return self;
        }
        else
        {
            self.log_den=0;

            // Return ownership
            return self;
        }
    }

    fn logden(self) -> Int
    {
        return self.log_den;
    }

    fn local_gen() -> Self
    {
        return Self
        {
            num: T::from(1),
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
where T:Localizable+Add<Output=T>+PartialEq+From<Int>+Copy
{
    type Output = Self;
    fn add(self, other: Self) -> Self {
        if other.num==T::from(0)
        {
            return self;
        }
        if self.num==T::from(0)
        {
            return other;
        }
        // println!("Want to add height {} and height {}",self.log_den,other.log_den);
        if self.log_den>other.log_den 
        {
            let mut other_temp = self.num;
            // println!("Will invoke perform_n_multiplications with n={}",self.log_den-other.log_den);
            other_temp = other_temp.perform_n_multiplications(self.log_den-other.log_den);
            let temp = 
                Self{
                    num: self.num + other_temp,
                    log_den: self.log_den
                };

            // return temp; // Possible speedup?
            // temp.fix();
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
where T:Mul<Output=T>+Localizable+PartialEq+From<Int>+Copy
{
    type Output = Self;

    fn mul(self, other: Self) -> Self 
    {
        let temp = Self
        {
            num: self.num*other.num,
            log_den: self.log_den+other.log_den
        };
        // return temp;
        temp.fix();
        return temp;
    }
}



// Teaching rust how to subtract Dyad elements
impl<T> Sub for Local<T> 
where T:Neg<Output=T>+Localizable+PartialEq+Add<Output=T>+From<Int>+Copy
{

    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self+(-other) //subtraction is just adding the additive inverse
    }

}


// Nicely display Local<T> elements
impl<T> Display for Local<T>
where T: Display+From<Int>+PartialEq
{
    fn fmt(&self, f: &mut Formatter) -> Result{
        if self.log_den>0  
        {   
            if self.num== T::from(1)
            {
                write!(f,"\\sqrt(2)^(-{})", self.log_den)
            }
            else if self.num == T::from(0)
            {
                write!(f,"0")
            }
            else 
            {
                write!(f,"({})*\\sqrt(2)^(-{})",self.num, self.log_den)
            }
        }        
        else if self.log_den<0  
        { 
            if self.num== T::from(1)
            {
                write!(f,"\\sqrt(2)^({})",self.log_den)
            }
            else if self.num == T::from(0)
            {
                write!(f,"0")
            }
            else 
            {
                write!(f,"({})*\\sqrt(2)^({})",self.num, self.log_den)
            }
        }
        else 
        {   
            write!(f,"({})",self.num) 
        }
    }
}


// Teaching rust how to compare these ring elements
impl<T> PartialEq for Local<T> 
where T: Localizable+PartialEq+From<Int>+Copy
{
    fn eq(&self, other: &Self) -> bool {
        if other.num.is_divisible() && other.num != T::from(0)
        {
            other.fix();
        }
        if self.num.is_divisible() && self.num != T::from(0)
        {
            self.fix();
        }
        return self.num==other.num && self.log_den==other.log_den;
    }
}



// To construct Local<T> directly from integers
impl<T> From<Int> for Local<T> 
where T: From<Int>+Localizable+PartialEq+Copy
{
    fn from(input : Int) -> Self {

        let out = Self{
            num: T::from(input.try_into().unwrap()),
            log_den: 0
        };

        if input!=1 && input!=0
        {
            out.fix();
        }

        return out;
    }
}


// To construct Local<T> from T
impl<T> From<T> for Local<T> 
where T: Localizable+PartialEq+Copy+From<Int>
{
    fn from(input : T) -> Self {

        let out = Self{
            num: input,
            log_den: 0
        };

        if input!=T::from(1) && input!=T::from(0)
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
    fn conj(self) -> Self {
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
where T: Mul<Output=T>+Conj+Neg<Output=T>,
      T: Localizable+PartialEq+From<Int>+Copy
{
    type Output = Self;

    fn div(self, other: Self) -> Self 
    {
        let norm = other.num.norm();
        if norm!=1 && norm!=-1
        {
            panic!("Division by a non-unit")
        }
        else if norm==1
        {
            let temp = Self
            {
                num: self.num*other.num.conj(),
                log_den: self.log_den-other.log_den
            };
            // return temp;
            temp.fix();
            return temp;
        }
        else
        {
            let temp = Self
            {
                num: -self.num*other.num.conj(),
                log_den: self.log_den-other.log_den
            };
            // return temp;
            temp.fix();
            return temp;
        }
    }
}

