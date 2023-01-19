use crate::structs::rings::Int;
use num_traits::Zero;




// Return closest integer to top/bottom
// Will send halfs to the ceiling
// Probably shouldn't matter
pub fn nearest_integer(top :Int, bottom: Int) -> Int
{
    if bottom == Int::zero()
    {
        panic!("What do you think?");
    }
    else if top ==Int::zero()
    {
        return Int::zero();
    }
    else if  top > Int::zero() && bottom > Int::zero() 
    { 
        let twice = (top << 1)/bottom;
        if twice%2==Int::zero()
        {
            return twice>>1;
        }
        else 
        {
            return (twice+1)>>1  ;
        }

    }
    else if  top < Int::zero() && bottom < Int::zero()
    {
        return nearest_integer(-top, -bottom);
    }
    else if  top > Int::zero() && bottom < Int::zero() 
    {
        return -nearest_integer(top, -bottom);
    }
    else 
    {
        return -nearest_integer(-top, bottom);
    }
}
