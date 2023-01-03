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
use num_traits::Pow;
use num_traits::Zero;
use num_traits::One;

// For display
use std::fmt::Result;
use std::fmt::Display;
use std::fmt::Formatter;


// Base matrix
use crate::structs::sunimat::SUniMat;
use crate::structs::rings::local_ring::Local;
use crate::structs::rings::zroot2::Zroot2;
use num_complex::Complex;
use crate::structs::rings::special_values::mu_8;
use crate::structs::rings::special_values::onebyroot2comp;
type KMMring = Complex<Local<Zroot2>>;

// Unitary matrices. They are the form 
// /                \
// | u  -t^* * (omega) |
// | t   u^* * (omega) |
// \                /
// where u and t are in giventype
// and omega is such that abs(omega) = 1
//
// NOTE that the when u,t, omega are in the KMM ring
// that is, when the unitary is from Clifford+T gate set
// then we know that omega is an 8th root of unity
// so we can instead save a single u8 value mod 8


type ExpInt = u8;

#[derive(Debug,Copy,Clone)]
pub struct ExactUniMat
{
    pub mat: SUniMat<KMMring>,       // This will save the u and t
    pub omega_exp: ExpInt            // This will actually just be a number mod 8
}



// Nicely display Unitary Matrices
impl Display for ExactUniMat
{
    fn fmt(&self, f: &mut Formatter) -> Result{
        // write!(f,"/       \\");
        // write!(f,"| {} {} |", self.u,-self.t.conj());
        write!(f,"|\t{}\t{}\t|\n|\t{}\t{}\t|\n",self.mat.u, -self.mat.t.conj()*mu_8().pow(self.omega_exp) , self.mat.t,self.mat.u.conj()*mu_8().pow(self.omega_exp))
    }
}


// Auxilary function to help in multiplication, inverse, etc.
impl SUniMat<KMMring>
{
    pub fn twist(self, no_of_turns : ExpInt) -> SUniMat<KMMring> 
    {
        Self{
            u : self.u,
            t : self.t*mu_8().pow(no_of_turns)
        }
    }
}


// Conjugate-transpose UniMat<T> elements
// Same as taking an inverse
impl ExactUniMat
{
    pub fn inv(self) -> ExactUniMat 
    {
        return ExactUniMat
        {
            mat: self.mat.inv().twist(self.omega_exp),
            omega_exp: 8-self.omega_exp,
        }
    }
}


// Teaching rust how to multiply UniMat<T> elements
// Also see this for why it looks so weird:
// https://stackoverflow.com/questions/39169795/error-when-using-operators-with-a-generic-type
impl Mul for ExactUniMat
{
    type Output = ExactUniMat;
    fn mul(self, other: ExactUniMat) -> ExactUniMat
    {
        Self{
            mat: self.mat*other.mat.twist(8-self.omega_exp),
            omega_exp: (self.omega_exp+other.omega_exp)%8,
        }
    }
}



// // Get zero and one as Unitary matrices
impl ExactUniMat
{
    // WARNING: zero is possible to construct, but avoid using it
    // It is not a unitary matrix
    pub fn zero() -> Self 
    {
        return Self{
            mat: SUniMat::<KMMring>
                { u: KMMring::zero(), t: KMMring::zero()},
            omega_exp : 0
        }

    }

    pub fn one() -> Self 
    {
        return Self
        {
            mat: SUniMat::<KMMring>
                { u: KMMring::one(), t: KMMring::zero()},
            omega_exp : 0
        }
    }

    pub fn h_gate() -> Self 
    {
        return Self
        {
            mat: SUniMat::<KMMring>
                { u: onebyroot2comp(), t: onebyroot2comp()},
            omega_exp : 4
        }

    }

    pub fn t_gate() -> Self
    {
        return Self
        {
            mat: SUniMat::<KMMring>
                { u: KMMring::one(), t: KMMring::zero() },
            omega_exp : 1
        };
    }

    pub fn from_string(gate_string : &String) -> Self
    {

        let mut output = ExactUniMat::one();



        for i in gate_string.chars()
        {
            if i == 'H'
            {
                output = ExactUniMat::h_gate()*output;
            }

            else if i == 'T'
            {
                output = ExactUniMat::t_gate()*output;
            }
            
            else 
            {
                panic!("Gates other than H or T in the sequence");
            }

        }

        return output;
    }
}


// Teaching rust how to compare these ring elements
impl PartialEq for ExactUniMat
{
    fn eq(&self, other: &Self) -> bool {
        return self.mat == other.mat && self.omega_exp%8 == other.omega_exp%8;
    }
}





