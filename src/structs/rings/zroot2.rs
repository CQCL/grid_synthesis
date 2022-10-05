// Quadratic number field with root 2
pub struct Zroot2(pub i64,pub i64); //a+b\sqrt(2)

// use std::ops::Neg; 
use std::ops::Add; // We bring them in so that we can overload the operators
                   // Rust must learn how to do arithmetics in our rings

use std::fmt; // To teach rust how to display our ring elements
                
                
// Rust must know how to diplay elements of this ring
impl fmt::Debug for Zroot2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.1<0 {
            write!(f, "{}{}\\sqrt{{2}}", self.0, self.1)
        }
        else
        {
            write!(f, "{}+{}\\sqrt{{2}}", self.0, self.1)
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
impl Localizable for Zroot2
{
    fn is_divisible(self) -> bool{
        if self.0%2==0:
            return true
    }
    fn perform_one_division(self) -> bool{
        // divide by sqrt2 if it is indeed divisible by sqrt2

    }

}


// Get zero and one as ring elements
impl Constructs for Zroot2
where T: Constructs<T>
{
    
    // WARNING: zero is possible to construct, but avoid using it
    // It is not a unitary matrix
    fn zero() -> Self {
        return Self{ u: T::zero(), t:T::zero()}

    }

    fn one() -> Self {
        return Self{ u: T::one(), t:T::zero()}
    }
}
