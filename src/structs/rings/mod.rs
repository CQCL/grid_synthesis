// This rust module contains all our rings etc.
// Some of these O(1) operations about multiplication and addition can be optimized a little bit
// The DOmega ring can be stored directly without any underlying Dyad implentation, I feel



// Supertrait
// They become various conjugation operators on our rings
pub trait Conj<T>{
    fn conj(self) -> Self;
}


// Supertrait (deprecated now, use from::(0) etc instead)
// They return zero and one elements of our rings
// pub trait Constructs<T>{
//     fn zero() -> Self;
//     fn one() -> Self;
// }

// Supertrait
// Allows you to define a possible prime ideal for localization
// See Zroot2 for an implementation
pub trait Localizable{
    // one should be able to check if the divisibility by ideal exists
    fn is_divisible(self) -> bool;
    
    // if it is divisible, we perform the division
    // We return the number of times we divided
    fn reduce_by_dividing(self) -> u32;

    // Multiply the ideal generator by the number of times given by n
    // WARNING: Bit Overflow may occur here
    fn perform_n_multiplications(self,_:i64) -> ();
}



pub mod zroot2;
pub mod local_ring;
pub mod domega;
pub mod dyad;
pub mod unimat;
pub mod complex;
