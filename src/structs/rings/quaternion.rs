// Take a ring
// And add i, j, k etc to it
// This will create a Quaternion algebra over anything

use crate::structs::rings::Conj; //Conjugation trait


use crate::structs::rings::Int; //Integer type standard

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


// Quaternions.
// They are of the form a+ib+cj+dk
// where a,b,c,d are some ring elements
// /         \
// | u  -t^* |
// | t   u^* |
// \         /
// where u and t are in giventype
// It is required that giventype has implementations of Add,Sub,Mult,Neg,conj
#[derive(Debug,Copy,Clone)]
pub struct Quaternion<T>(T,T,T,T);


// Nicely display Unitary Matrices
impl<T> Display for Quaternion<T>
where T:Display
{
    fn fmt(&self, f: &mut Formatter) -> Result{
        // write!(f,"/       \\");
        // write!(f,"| {} {} |", self.u,-self.t.conj());
        write!(f,"{}+{}*I+{}*J+{}*K", self.0, self.1 , self.2,self.3)
    }
}

// Conjugate-transpose Quaternion<T> elements
// Same as taking an inverse
impl<T> Neg for Quaternion<T> 
where T: Neg<Output=T>
{
    type Output = Self;
    fn neg(self) -> Self {
        Self{
            0: -self.0,
            1: -self.1,
            2: -self.2,
            3: -self.3
        }
    }
}

// Conjugate Complex elements
impl<T,E> Conj<E> for Quaternion<T> 
where T: Neg<Output=T>
{
    fn conj(self) -> Self {
        Self{
            0: self.0,
            1: -self.1,
            2: -self.2,
            3: -self.3,
        }
    }
}




// Teaching rust how to multiply Quaternion<T> elements
// Also see this for why it looks so weird:
// https://stackoverflow.com/questions/39169795/error-when-using-operators-with-a-generic-type
impl<T> Mul for Quaternion<T> 
where T: Copy+Mul<Output=T>+Neg+Add<Output=T>+Sub<Output=T>
{
    type Output = Quaternion<T>;
    fn mul(self, other: Quaternion<T>) -> Quaternion<T>
    {
        // NAIVE ALGORITHM:
        //
        // One can reduce the number of multiplications:
        // https://math.stackexchange.com/questions/1103399/alternative-quaternion-multiplication-method
        Self{
            0: self.0*other.0-self.1*other.1-self.2*other.0-self.3*other.2,
            1: self.0*other.1+self.1*other.0+self.2*other.3-self.3*other.1,
            2: self.0*other.2-self.1*other.3+self.2*other.2+self.3*other.0,
            3: self.0*other.3+self.1*other.2-self.2*other.1+self.3*other.3,
        }
    }
}



// Get zero and one as Unitary matrices
impl<T> From<Int> for Quaternion<T> 
where T: From<Int> 
{
    fn from(int: Int) -> Self 
    {
        Self
        {
            0: T::from(int),
            1: T::from(0),
            2: T::from(0),
            3: T::from(0)
        }
    }
}


// Teaching rust how to compare these ring elements
impl<T> PartialEq for Quaternion<T> 
where T: PartialEq
{
    fn eq(&self, other: &Self) -> bool {
        return self.0==other.0 && self.1==other.1 && self.2==other.2 && self.3==other.3;
    }
}

// Teaching rust how to add Quaternion elements
impl<T> Add for Quaternion<T>
where T: Add<Output=T>
{
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Quaternion(self.0+other.0,self.1+other.1,self.2+other.2,self.3+other.3)
    }
}


// Teaching rust how to subtract Quaternion elements
impl<T> Sub for Quaternion<T> 
where T: Sub<Output=T>
{
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self{
            0: self.0-other.0,
            1: self.1-other.1,
            2: self.2-other.2,
            3: self.3-other.3,
        }
    }
}



//Will work only if T has a division implementation
impl<T> Quaternion<T>
where T: Mul<Output=T> + Add<Output=T>,
      T: Copy+PartialEq
{
    // reduced square norm of our quaternion algebra
    // Equal to whatever is the output here
    pub fn rsqnorm(self) -> T
    {
        return self.0*self.0+self.1*self.1+self.2*self.2+self.3*self.3
    }
}



//Will work only if T has a division implementation
impl<T> Quaternion<T>
where T: Mul<Output=T> + Add<Output=T> + Neg<Output=T> + Div<Output=T>,
      T: Copy+PartialEq+From<Int>
{

    pub fn inv(self) -> Self
    {
        let norm = self.rsqnorm();
        if norm==T::from(0)
        {
            panic!("Division by zero in Quaternions!!!!");
        }
        else
        {
            Self{
               0:  self.0/norm,
               1: -self.1/norm,
               2: -self.2/norm,
               3: -self.3/norm,
            }
        }


    }

}

//Will work only if T has a division implementation
impl<T> Div for Quaternion<T>
where T: Mul<Output=T> + Add<Output=T> + Neg<Output=T> + Div<Output=T>+Sub<Output=T>,
      T: Copy+PartialEq+From<Int>
{

    type Output = Self;
    fn div(self, other: Self ) -> Self
    {
        return self.mul(other.inv())
    }
}
