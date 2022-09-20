// Unitary matrices
// That's what quantum gates are
// It is a Lie group over the reals of dimension 3
// It is equivalently, also the space of unit quaternions
// We are interested in unitary matrices over the ring D[\omega]


// We bring them in so that we can overload the operators
// Rust must learn how to do arithmetics in our rings
use std::ops::Mul; 


use crate::structs::rings::domega::DOmega; //Dyads form DOmega ring elements. This might change.
                                           //
                                           //

// Unitary matrices. They are of the form
// /         \
// | u  -t^* |
// | t   u^* |
// \         /
// where u and t are in D_omega
pub struct UniMat{
    pub u: DOmega,
    pub t: DOmega,
}


// Conjugate-transpose UniMat elements
// Same as taking an inverse
impl UniMat {
    pub fn inv(self) -> UniMat {
        UniMat{
        u : self.u.conj(),
        t : -self.t
        }
    }
}


// Teaching rust how to multiply UniMat elements
impl Mul for UniMat {
    type Output = UniMat;

    fn mul(self, other: UniMat) -> UniMat {
        UniMat{
            u: self.u*other.u - self.t.conj()*other.t,
            t: self.t*other.u+self.u.conj()*other.t
        }
    }
}



