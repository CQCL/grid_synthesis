// Unitary matrices
// That's what quantum gates are
// It is a Lie group over the reals of dimension 3
// It is equivalently, also the space of unit quaternions
// We are interested in unitary matrices over the ring D[\omega]


use crate::structs::rings::Conj; //Conjugation trait

// We bring them in so that we can overload the operators
// Rust must learn how to do arithmetics in our rings
use std::ops::Neg; 
use std::ops::Add; 
use std::ops::Sub; 
use std::ops::Mul; 


// use crate::structs::rings::domega::DOmega; //Dyads form DOmega ring elements. This might change.
// use crate::structs::rings::complex::Complex; //Dyads form DOmega ring elements. This might change.


// Unitary matrices. They are of the form
// /         \
// | u  -t^* |
// | t   u^* |
// \         /
// where u and t are in giventype
// It is required that giventype has implementations of Add,Sub,Mult,Neg,conj
pub struct UniMat<T>{
    pub u: T,
    pub t: T,
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

