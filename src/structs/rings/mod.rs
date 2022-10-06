// This rust module contains all our rings etc.
// Some of these O(1) operations about multiplication and addition can be optimized a little bit
// The DOmega ring can be stored directly without any underlying Dyad implentation, I feel



// Supertrait
// They become various conjugation operators on our rings
pub trait Conj<T>{
    fn conj(self) -> Self;
}


// Supertrait
// They return zero and one elements of our rings
pub trait Constructs<T>{
    fn zero() -> Self;
    fn one() -> Self;
}

// Supertrait
// Allows you to define a possible prime ideal for localization
pub trait Localizable{
    // one should be able to check if the divisibility by ideal exists
    fn is_divisible(self) -> bool;
    
    // if it is divisible, we perform the division
    fn perform_one_division(self) -> ();

    // Multiply the ideal generator once
    fn perform_n_multiplications(self,u32 n) -> ();
}
pub mod zroot2;
pub mod local_ring;
pub mod domega;
pub mod dyad;
pub mod unimat;
pub mod complex;
