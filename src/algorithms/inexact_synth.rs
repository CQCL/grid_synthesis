// Will implement the Ross Selinger algorithm here
// This is the content of 1403.2975v3 Section 7.3 on page 14
//
//
//

use crate::structs::rings::Float;
use crate::structs::rings::Int;


// PLEASE IMPROVE FOR BETTER ACCURACY!
const sqrt2:Float = 1.4142135623730951;

type Comp = num_complex::Complex<Float>;
type Mat2 = nalgebra::Matrix2<Float>;
type Mat4 = nalgebra::Matrix4<Float>;
type Vec4 = nalgebra::Matrix4x1<Float>;

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

pub fn call_lll_on_nalgebra_matrix( m : Mat4) -> Mat4
{
    use lll_rs::{
        // l2::{bigl2, l2f},
        lll::lllf,
        matrix::Matrix,
        vector::VectorF,
    };

    // use rug::{Integer,Assign};

    // Init the matrix with Integer
    let mut basis: Matrix<VectorF> = Matrix::init(4, 4);

    // Populate the matix
    basis[0] = VectorF::from_vector(vec![
                                    m[(0,0)],
                                    m[(1,0)],
                                    m[(2,0)],
                                    m[(3,0)],
    ]);

    basis[1] = VectorF::from_vector(vec![
                                    m[(0,1)],
                                    m[(1,1)],
                                    m[(2,1)],
                                    m[(3,1)],
    ]);
    basis[2] = VectorF::from_vector(vec![
                                    m[(0,2)],
                                    m[(1,2)],
                                    m[(2,2)],
                                    m[(3,2)],
    ]);
    basis[3] = VectorF::from_vector(vec![
                                    m[(0,3)],
                                    m[(1,3)],
                                    m[(2,3)],
                                    m[(3,3)],
    ]);

    // Perfom the LLL basis redution
    lllf::lattice_reduce(&mut basis);

    // OR
    // Perfom the LLL basis redution
    // Specify the delta and eta coefficient for the reduction
    // bigl2::lattice_reduce(&mut basis, 0.5005, 0.999);
    //
    return Mat4::new( 
        basis[0][0], basis[0][1], basis[0][2], basis[0][3],
        basis[1][0], basis[1][1], basis[1][2], basis[1][3],
        basis[2][0], basis[2][1], basis[2][2], basis[2][3],
        basis[3][0], basis[3][1], basis[3][2], basis[3][3],
        );

}

pub fn find_lattice_points_in_a_box(red_matrix: Mat4, center: Vec4, r: Float) -> ()
{
    // Given a box of radius r (half of the side length = radius) 
    // This will list down all the points in the box
    // assuming that the lattice is LLL-reduced
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

    let mat_lattice = 
        Mat4::new( 1.0 ,   sqrt2  , 0.0  , 0.0,
                   1.0 ,  -sqrt2  , 0.0  , 0.0,
                   0.0 ,      0.0 , 1.0  ,  sqrt2 ,
                   0.0 ,      0.0 , 1.0  , -sqrt2 );

    let mut reducable = mat_big*mat_lattice;

    let reduced = call_lll_on_nalgebra_matrix(reducable);
    
    
    let center4d = reduced.try_inverse().unwrap()*Vec4::new(center.re,center.im,0.0,0.0);
    let radius_big = 1.0;

    // List lattice points in a four dimensional ball 
    // Given the the basis is LLL-reduced
    // and given a center and radius
    find_lattice_points_in_a_box(reduced, center4d, radius_big);

    

    // N should be calculated based on LLL parameters
    //
    let N = 1;

    for i1 in -N..N 
    {   for i2 in -N..N
        {
            for i3 in -N..N
            {
                for i4 in -N..N
                {
                    let point = reduced*Vec4::new(i1 as Float,i2 as Float,i3 as Float,i4 as Float);
                    if point.norm()<=1.0
                    {
                        println!("{}", point);
                        println!("{},{},{},{}",i1,i2,i3,i4 );
                    }
                }

            }
        }

    }

    //////////////////// DEBUG ZONE  /////////////////////
    //
    // println!("{}", reduced);
    //
    ///////////////////END OF DEBUG ZONE /////////////////

}

