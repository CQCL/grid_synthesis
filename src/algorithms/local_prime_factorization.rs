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




pub fn prime_factorization_of_int_local( input : Local::<Int> ) -> ()
{
    let number  = input.num as u32;
    let factorvec = Factorization::run(number);
}


pub fn prime_factorization_of_loc( input: Local::<Zroot2> ) -> ()
{
    todo!();
}
