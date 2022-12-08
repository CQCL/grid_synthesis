use crate::structs::rings::Int;




// Return closest integer to top/bottom
// Will send halfs to the ceiling
// Probably shouldn't matter
// TODO: Optimize this function
pub fn nearest_integer(top :Int, bottom: Int) -> Int
{
    if bottom == 0
    {
        panic!("What do you think?");
    }
    else if top ==0
    {
        return 0;
    }
    else if ( top > 0 && bottom > 0 )     
    { 
        let twice = (top << 1)/bottom;
        if twice%2==0
        {
            return twice>>1;
        }
        else 
        {
            return (twice+1)>>1  ;
        }

    }
    else if ( top < 0 && bottom < 0 ) 
    {
        return nearest_integer(-top, -bottom);
    }
    else if ( top > 0 && bottom < 0 ) 
    {
        return -nearest_integer(top, -bottom);
    }
    else 
    {
        return -nearest_integer(-top, bottom);
    }
}
