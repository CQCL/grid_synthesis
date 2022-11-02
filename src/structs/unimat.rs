// Unitary matrices


use crate::structs::rings::Conj; //Conjugation trait
// use crate::structs::rings::Constructs; //Construction trait

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


// Unitary matrices. They are the form 
// /                \
// | u  -t^*(omega) |
// | t   u^*(omega) |
// \                /
// where u and t are in giventype
// and omega is a unit norm object
#[derive(Debug,Copy,Clone)]
pub struct UniMat<T>
{
    pub u: T,
    pub t: T,
}


// Nicely display Unitary Matrices
impl<T> Display for UniMat<T>
where T:Neg<Output=T>+Conj<T>+Display+Copy
{
    fn fmt(&self, f: &mut Formatter) -> Result{
        // write!(f,"/       \\");
        // write!(f,"| {} {} |", self.u,-self.t.conj());
        write!(f,"|\t{}\t{}\t|\n|\t{}\t{}\t|\n",self.u, -self.t.conj() , self.t,self.u.conj())
    }
}

// Conjugate-transpose UniMat<T> elements
// Same as taking an inverse
impl<T> UniMat<T> 
where T: Neg<Output=T>+Conj<T>
{
    pub fn inv(self) -> UniMat<T> {
        Self{
            u : self.u.conj(),
            t : -self.t
        }
    }
}


// Teaching rust how to multiply UniMat<T> elements
// Also see this for why it looks so weird:
// https://stackoverflow.com/questions/39169795/error-when-using-operators-with-a-generic-type
impl<T> Mul for UniMat<T> 
where T: Copy+Mul<Output=T>+Conj<T>+Neg+Add<Output=T>+Sub<Output=T>
{
    type Output = UniMat<T>;
    fn mul(self, other: UniMat<T>) 
        -> UniMat<T>
    {
        Self{
            u: self.u*other.u-self.t.conj()*other.t,
            t: self.t*other.u+self.u.conj()*other.t
        }
    }
}



// Get zero and one as Unitary matrices
impl<T> UniMat<T>
where T: From<i32>
{
    // WARNING: zero is possible to construct, but avoid using it
    // It is not a unitary matrix
    pub fn zero() -> Self {
        return Self{ u: T::from(0), t: T::from(0)}

    }

    pub fn one() -> Self {
        return Self{ u: T::from(1), t: T::from(0)};
    }

    pub fn h_gate() -> Self 
    {
        return Self{ u: T::from(0), t: T::from(1)};
    }

    pub fn t_gate() -> ()
    {
        // TODO
        // These gates may not have determinant 1
    }
}


// Teaching rust how to compare these ring elements
impl<T> PartialEq for UniMat<T> 
where T: PartialEq
{
    fn eq(&self, other: &Self) -> bool {
        return self.u==other.u && self.t==other.t;
    }
}
