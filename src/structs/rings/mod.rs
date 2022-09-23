// This rust module contains all our rings etc.
// Some of these O(1) operations about multiplication and addition can be optimized a little bit
// The DOmega ring can be stored directly without any underlying Dyad implentation, I feel
//

pub trait Conj<T>{
    fn conj(self) -> Self;
}


pub mod zroot2;
pub mod domega;
pub mod dyad;
pub mod unimat;
pub mod complex;
