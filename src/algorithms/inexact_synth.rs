// Will implement the Ross Selinger algorithm here
// This is the content of 1403.2975v3 Section 7.3 on page 14
//
//
//

use crate::structs::rect::Rect;
use crate::structs::rings::Float;
use crate::structs::rings::Int;

type Comp = num_complex::Complex<Float>;
// type Mat = nalgebra::Matrix4<Float>;


pub fn ellipse_parameters_for_region_A(direction: Comp, epsilonA: Float, maxlogdep: Int) -> ()
{
    // Let direction = e_1 + i e_2
    // Want to have x^2 + y^2 <1 and 1-(xe_1 + ye_2) < \epsilon
    //
    // The ellipse is then the region
    // c(x^2+y^2) + (1-c) (  1- (xe_1+ye_2) )^2 < c + (1-c) \epsilon^2
    //
    // This ellipse has smallest volume when c= \epsilon 
    // 
    //
    // A two dimensional ellipsoid:
    // 1. has a center
    // Ans: 
    // It is at (e_1+ie_2)*(1-epsilon)
    //
    //
    // 2. has a 2x2 symmetric positive definite matrix giving an inner product
    // Ans: 
    // It is   
    // [ x y ] [ c+e_1^2(1-c)    e_1e_2c      ] [ x ]  
    //         [  e_1 e_2 c    c+e_2^2(1-c)   ] [ y ] 
    // where c= \epsilon
    //
    //
    // 3. Has a radius in terms of this inner product
    // Ans:  Radius is 2c^2 - c^3 where c = \epsilon
    //
    
    todo!();
}


// This is an implementation of Proposition 5.22
pub fn grid_problem( direction: Comp, epsilonA: Float  ,  maxlogdep: Int )-> ()
{


    // Trying to find the lattice points in the intersection of a disc and a half plane
    // Zoom out for the picture
    // +----------------------------------------------------------------------------------------------------+
    // |                                                                                                    |
    // |                                                                                                    |
    // |                                                                                                    |
    // |                                                                                                    |
    // |                                                                                                    |
    // |                                   /----------\                                                     |
    // |                            -------            --------                                             |
    // |                          -/                           \-                                           |
    // |                        -/                               \-                                         |
    // |                      -/                                   \-                                       |
    // |                    -|                                       |-                                     |
    // |                     |                                       |                                      |
    // |                    /                                         \                                     |
    // |                    |                                         |                                     |
    // |                    |                                         |                                     |
    // |                   /                                           \                   ------           |
    // |                   |                                           |         ---------/                 |
    // |                   \                                           ---------/                           |
    // |                    |                               ----------/                                     |
    // |                    \                     ---------/          /                                     |
    // |                     |         ----------/                   |                                      |
    // |                    ----------/                              |-                                     |
    // |           ---------/ -\                                   /-                                       |
    // |     -----/             -\                               /-                                         |
    // |                          -\                           /-                                           |
    // |                            -------            --------                                             |
    // |                                   \----------/                                                     |
    // |                                                                                                    |
    // |                                                                                                    |
    // |                                                                                                    |
    // +----------------------------------------------------------------------------------------------------+
                                                                                                      



    todo!();
}

