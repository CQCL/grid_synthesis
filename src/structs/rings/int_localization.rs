// In this code, I will make the Int localizable
// Mostly just to see that it works
// We are localizing integers with respect to the prime 2
//
// This should create some dyadic integers



//Integer type is set globally
use crate::structs::rings::Int; 
use crate::structs::rings::Localizable; 


impl Localizable for Int
{
    fn is_divisible(self) -> bool
    {
        return self%2==0
    }

    fn reduce_by_dividing(self) -> (Self,Int)
    {
        let trailing = self.trailing_zeros();
        return (self >> trailing, trailing.try_into().unwrap());
    }
    
    fn perform_n_multiplications(self, n: Int) -> Self
    {
        return self << n;
    }
}
