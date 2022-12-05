//  What we really want is to prime factorize a Loc

// This can be done by either rewriting prime factorization algorithm
// for this ring. This would make me learn a lot, but would quite increase
// the duration of my internship. 
// The advantage would be to reduce external dependencies and perhaps speedup?
//
//
//

use prime_factorization::Factorization;


use crate::structs::rings::Int; 
use crate::structs::rings::local_ring::Local; 
use crate::structs::rings::zroot2::Zroot2; 




pub fn attempt_to_write_this_number_as_sum_of_two_squares_in_loc(our_num: Loc)  -> Option::<(Loc,Loc)>
{
    prime_factorization_of_loc(our_num);
    return None;
}


pub fn prime_factorization_of_loc( input: Local::<Zroot2> ) -> ()
{
        let num = input.num.norm() as u128;
        let factorvec = Factorization::run(num).prime_factor_repr();
    
        for (prime,power) in factorvec
        {
            todo!();
        }

}
