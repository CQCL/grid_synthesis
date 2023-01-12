use crate::structs::rings::local_ring::Local;
use crate::structs::rings::zroot2::Zroot2;
use crate::structs::rings::zomega::Zomega;
use num_complex::Complex;

use num_traits::One;
use num_traits::Zero;


// Returning the eigth root of unity
type Loc = Local<Zroot2>;
type Comp = Complex<Loc>;


pub fn onebyroot2comp() -> Comp
{
    return Comp
    {
        re: Loc
        {
            num: Zroot2::one(),
            log_den: 1
        },

        im: Loc
        {
            num: Zroot2::zero(),
            log_den: 0
        }

    };

}

pub fn sqrt2loc() -> Loc
{
    return Loc
        {
            num: Zroot2::one(),
            log_den: -1
        };

}

pub fn iota_zomega() -> Zomega
{
    return Zomega(0,0,1,0);
}

pub fn omega() -> Zomega
{
    return Zomega(0,1,0,0);
}

pub fn onebyroot2loc() -> Loc
{
    return Loc
        {
            num: Zroot2::one(),
            log_den: 1
        };

}


pub fn mu_8() -> Comp
{
    return Comp
    {
        re: Loc
        {
            num: Zroot2::one(),
            log_den: 1
        },

        im: Loc
        {
            num: Zroot2::one(),
            log_den: 1
        }
    };
}

pub fn sqrt2() -> Comp
{
    return Comp
    {
        re: Loc{
            num: Zroot2::one(),
            log_den: -1,
        },
        im: Loc::zero(),
    };
}

pub fn sqrtminus1() -> Comp
{
    return Comp{
        re: Loc::zero(),
        im: Loc::one(),
    }
}
