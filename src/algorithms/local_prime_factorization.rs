//  What we really want is to prime factorize a Loc

// This can be done by either rewriting prime factorization algorithm
// for this ring. This would make me learn a lot, but would quite increase
// the duration of my internship. 
// The advantage would be to reduce external dependencies and perhaps speedup?
//
//
//
//
use num_traits::One;
use num_traits::pow;
use num_traits::NumCast;


use prime_factorization::Factorization;


use crate::structs::rings::Int; 
use crate::structs::rings::LogDepInt; 
use crate::structs::rings::local_ring::Local; 
use crate::structs::rings::zroot2::Zroot2; 
use crate::structs::rings::special_values::sqrt2loc;

type Loc = Local::<Zroot2>;

type FactorInt = u64;
type FactorPowerInt = u32;



pub fn attempt_to_write_this_number_as_sum_of_two_squares_in_loc(our_num: Loc)  -> Option::<(Loc,Loc)>
{
    let factorvec = prime_factorization_of_loc(our_num);

    let mut left = Loc::one();
    let mut right= Loc::one();

    for (prime, locprime,power) in factorvec
    {
        if power%2==1
        {
            if prime%8==7
            {
                return None;
            }
        }
        else
        {
            let powerby2 = (power >> 1);
            left = left*pow(locprime,powerby2.try_into().unwrap());
            right = right*pow(locprime,powerby2.try_into().unwrap());
        }
    }
    return None;


}

// Will keep track of base primes, prime above in Zroot2 and the multiplicty it has
pub fn prime_factorization_of_loc( input: Local::<Zroot2> ) -> Vec::<( FactorInt, Loc, LogDepInt)>
{
        let num = input.num.norm() as FactorInt;
        let factorvec = Factorization::run(num).prime_factor_repr();
        let mut factorvecloc = vec!((2, sqrt2loc() , -input.log_den));
        for (prime,power) in factorvec
        {
            if prime%7==3 || prime%7==5
            {
                let primeloc = <Loc as NumCast >::from(prime).unwrap();
                let powerby2 = (power >> 1);
                factorvecloc.push( (prime,primeloc,powerby2.try_into().unwrap() ) );
            }
            // if prime is not 3,5 modulo 8, then it is 1,7 modulo 8
            else
            {
                // As suggested in the preprint
                // we have primeloc = gcd(prime,x^2-2)
            }
        }

        return factorvecloc;

}
