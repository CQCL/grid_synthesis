// This rust module contains all our rings etc.
// Some of these O(1) operations about multiplication and addition can be optimized a little bit
// The DOmega ring can be stored directly without any underlying Dyad implentation, I feel
//


// They become various conjugation operators on our rings
pub trait Conj<T>{
    fn conj(self) -> Self;
}


// They return zero and one elements of our rings
pub trait Constructs{
    fn zero() -> Self;
    fn one() -> Self;
}

pub mod zroot2;
pub mod domega;
pub mod dyad;
pub mod unimat;
pub mod complex;
