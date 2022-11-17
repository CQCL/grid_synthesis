// Will implement the Ross Selinger algorithm here
// This is the content of 1403.2975v3 Section 7.3 on page 14
//
//
//

use crate::structs::rect::Rect;
use crate::structs::rings::Float;
use crate::structs::rings::Int;


// PLEASE IMPROVE FOR BETTER ACCURACY!
const sqrt2:Float = 1.4142135623730951;

type Comp = num_complex::Complex<Float>;
type Mat2 = nalgebra::Matrix2<Float>;
type Mat4 = nalgebra::Matrix4<Float>;

// See comments for output idea
pub fn ellipse_parameters_for_region_A(direction: Comp, epsilonA: Float ) -> (Comp, Mat2, Float)
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
    // 2. has a 2x2 symmetric positive definite matrix giving an inner product
    // Ans: 
    // It is   
    // [ x y ] [ c+e_1^2(1-c)     e_1e_2(1-c)    ] [ x ]  
    //         [  e_1e_2(1-c)     c+e_2^2(1-c)   ] [ y ] 
    // where c= \epsilon
    // This matrix has determinant = c
    // This matrix has trace  = 1+c
    // the square root of this matrix should be (check)
    //  [ c+e_1^2(1-c)+sqrt(c)              e_1e_2(1-c) ]  * (1+sqrt(c))^-1
    //  [  e_1e_2(1-c)            c+e_2^2(1-c)+sqrt(c)  ] 
    //
    // 3. Has a radius in terms of this inner product
    // Ans:  Radius is 2c^2 - c^3 where c = \epsilon


    let c = epsilonA;
    let sqrtc = c.sqrt();
    let e_1 = direction.re;
    let e_2 = direction.im;
    // todo!();
    
    let mut outputmatrix =  nalgebra::matrix![
     c+e_1*e_1*(1.0-c)+sqrtc   , e_1*e_2*(1.0-c)   ;
     e_1*e_2*(1.0-c)  , c+e_2*e_2*(1.0-c)+sqrtc  ;
     ];

     outputmatrix *= 1.0/(1.0+sqrtc);

     return (direction*(1.0-c), outputmatrix, c*c*(2.0-c));



}



// This is an implementation of Proposition 5.22
pub fn grid_problem( direction: Comp, epsilonA: Float,  maxlogdep: Int )-> ()
{

    // Trying to find the lattice points in the intersection of a disc and a half plane
    // Somtimes this region is called the "miniscule" region
    // Zoom out for the picture
    // +------------------------------------------------------------+
    // |                                                            |
    // |                                                            |
    // |                      /----------\                          |
    // |               -------            --------                  |
    // |             -/                           \-                |
    // |           -/                               \-              |
    // |         -/                                   \-            |
    // |       -|                                       |-          |
    // |        |                                       |           |
    // |       /                                         \          |
    // |       |                                         |          |
    // |       |                                         |          |
    // |      /                                           \         |
    // |      |                                           |         |
    // |      \                                           --------- |
    // |       |                               ----------/          |
    // |       \                     ---------/          /          |
    // |        |         ----------/                   |           |
    // |       ----------/                              |-          |
    // |  -----/ -\                                   /-            |
    // |           -\           THIS REGION         /-              |
    // |             -\                           /-                |
    // |               -------            --------                  |
    // |                      \----------/                          |
    // |                                                            |
    // |                                                            |
    // +------------------------------------------------------------+
                                                                                                      
    // Bound the region in an ellipse
    let (center,mat,radius) =   ellipse_parameters_for_region_A(direction, epsilonA);

    //create a 4 dimensional ellipsoid binding the region above combined with unit disc 
    
    let c = 1.0/(1.0+radius);
    
    let mat_big = 
             Mat4::new(mat[(0,0)]*c, mat[(0,1)]*c, 0.0  , 0.0,
                       mat[(1,0)]*c, mat[(1,1)]*c, 0.0  , 0.0,
                               0.0 ,         0.0 , 1.0-c, 0.0,
                               0.0 ,         0.0 , 0.0  , 1.0-c);
    let radius_big = 1.0;

    let mat_lattice = 
             Mat4::new( 1.0 ,   sqrt2  , 0.0  , 0.0,
                        1.0 ,  -sqrt2  , 0.0  , 0.0,
                        0.0 ,      0.0 , 1.0  ,  sqrt2 ,
                        0.0 ,      0.0 , 1.0  , -sqrt2 );

    let mut reducable = mat_big*mat_lattice;
    lll_rs::lll::biglll::lattice_reduce(&mut reducable);

    todo!("List lattice points in this region");
}

