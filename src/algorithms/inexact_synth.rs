// Will implement the Ross Selinger algorithm here
// This is the content of 1403.2975v3 Section 7.3 on page 14
//
//
//

//

use crate::structs::rings::Float;
use crate::structs::rings::local_ring::Local;
use crate::structs::rings::zroot2::Zroot2;
use crate::structs::rings::Int;
use crate::structs::rings::LogDepInt;
use crate::structs::rings::Localizable;


use crate::algorithms::local_prime_factorization::prime_factorization_of_loc;
use crate::structs::sunimat::SUniMat;


use num_traits::Pow;
use num_traits::One;
use nalgebra::linalg::QR;

// PLEASE IMPROVE FOR BETTER ACCURACY!
const SQRT2:Float = 1.4142135623730951;

type Comp = num_complex::Complex<Float>;
type Loc = Local::<Zroot2>;
type CompLoc = num_complex::Complex<Loc>;
type Mat2 = nalgebra::Matrix2<Float>;
type Mat4 = nalgebra::Matrix4<Float>;
type Vec4 = nalgebra::Matrix4x1<Float>;
type Vec4Int = nalgebra::Matrix4x1<Int>;
type Sunimat = SUniMat<Float>;
type Sunimatloc = SUniMat<Loc>;

// See comments for output idea
pub fn ellipse_parameters_for_region_a(direction: Comp, epsilon_a: Float ) -> (Comp, Mat2, Float)
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


    let c = epsilon_a;
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


pub fn find_operator_norm(input : Mat4) -> Float
{
    return input.singular_values()[0];
}



pub fn extract_gate_coordinate_in_local_ring(integer_coords: Vec4Int) -> Option::<(Loc,Loc)>
{
    let zrt_left = Zroot2(integer_coords[0],integer_coords[1]);
    let zrt_right = Zroot2(integer_coords[2],integer_coords[3]);

    // We throw away the points if they were covered while checking smaller values of exactlogdep
    if zrt_left.is_divisible() && zrt_right.is_divisible()
    {
        return None;
    }
    else 
    {
        return Some( (Loc::from_base(zrt_left), Loc::from_base(zrt_left)) );
    }
}


// Integer point should be useful
// That is, if it roughly lies in the region close to the required place,,
// Here we make sure that it is truly within the region we want
pub fn test_this_integer_point( this_point : Vec4Int, coordinate_basis: Mat4 , epsilon_a: Float, direction: Comp, exactlogdep: LogDepInt  ) -> bool
{

    // take the input and mutiply coordinate change matrix
    // and scale the vector by sqrt2^exactlogdep
    let actual_point = coordinate_basis*Vec4::new(this_point[0] as Float,this_point[1] as Float,this_point[2] as Float,this_point[3] as Float) * SQRT2.pow(-exactlogdep);

    let complex_point = Comp
    { 
        re: actual_point[0]+SQRT2*actual_point[1],
        im: actual_point[2]+SQRT2*actual_point[3],
    };

    let complex_point_dot_conj = Comp
    { 
        re: actual_point[0]-SQRT2*actual_point[1],
        im: actual_point[2]-SQRT2*actual_point[3],
    };

    if complex_point_dot_conj.norm_sqr() < 1.0
    {
        if complex_point.re*direction.re+complex_point.im*direction.im > 1.0- epsilon_a
        {
            return true;
        }
        else 
        {
            return false;
        }
    }
    else 
    {
        return false;
    }
}


pub fn attempt_to_write_this_number_as_sum_of_two_squares_in_loc(our_num: Loc)  -> Option::<(Loc,Loc)>
{
    prime_factorization_of_loc(our_num);
    return None;
}


// collection.push(Vec4Int::new( i1, i2, i3, i4));
// Return a complete unitary, if it can be returned
// Else return none
pub fn attempt_to_figure_out_gate_given_that_this_point_is_feasible( this_point: Vec4Int, exactlogdep: LogDepInt) -> Option::<Sunimatloc>
{
    
    let get_coord =  extract_gate_coordinate_in_local_ring(this_point);

    if get_coord == None
    {
        return None;
    }
    else
    {
        let (left,right) = get_coord.unwrap();

        let left_scaled = Loc
        {
            num: left.num,
            log_den: left.log_den + exactlogdep,
        };
        
        let right_scaled = Loc
        {
            num: right.num,
            log_den: right.log_den + exactlogdep,
        };

        let our_num = Loc::one() - left_scaled*left_scaled - right_scaled*right_scaled;

        attempt_to_write_this_number_as_sum_of_two_squares_in_loc(our_num);

        return None;

    }
    

}


pub fn consider(radius: Float, direction_of_rotation: Comp, epsilon_a: Float, exactlogdep: LogDepInt, int_center: Vec4Int, coordinate_basis: Mat4, this_point: Vec4Int) -> Option::<Sunimatloc>
{

    if test_this_integer_point( int_center+this_point, coordinate_basis, epsilon_a, direction_of_rotation, exactlogdep)
    {


        return attempt_to_figure_out_gate_given_that_this_point_is_feasible(int_center+this_point, exactlogdep);
    }
    else
    {
        return None;
    }
}



pub fn test_integer_points_in_ball_around_integer_center_of_radius(radius: Float, direction_of_rotation: Comp, epsilon_a: Float, coordinate_basis: Mat4, exactlogdep: LogDepInt, int_center: Vec4Int) -> Option<Sunimatloc>
{
    let mut collection = Vec::<Vec4Int>::new();

    let n = (radius*radius).floor() as Int + 1;



    let mut i1 : Int;
    let mut i2 : Int;
    let mut i3 : Int;
    let mut i4 : Int;

    i1 = 0;
    while i1*i1 < n 
    {
        i2 = 0;
        while i2*i2 < n - i1*i1 
        {
            i3 = 0;
            while i3*i3 < n - i1*i1 - i2*i2 
            {

                i4 = 0;
                while i4*i4 < n - i1*i1 - i2*i2 - i3*i3 
                {
                    if i1==0 && i2==0 && i3==0 && i4==0 
                    {
                        collection.push(Vec4Int::new( i1, i2, i3, i4));
                        collection.push(Vec4Int::new( i1, i2, i3, i4));
                    }
                    else if i1!=0 && i2==0 && i3==0 && i4==0 
                    {
                        collection.push(Vec4Int::new( i1, i2, i3, i4));
                        collection.push(Vec4Int::new(-i1, i2, i3, i4));
                    }
                    else if i1==0 && i2!=0 && i3==0 && i4==0 
                    {
                        collection.push(Vec4Int::new( i1, i2, i3, i4));
                        collection.push(Vec4Int::new( i1,-i2, i3, i4));
                    }
                    else if i1==0 && i2==0 && i3!=0 && i4==0 
                    {
                        collection.push(Vec4Int::new( i1, i2, i3, i4));
                        collection.push(Vec4Int::new( i1, i2,-i3, i4));
                    }
                    else if i1==0 && i2==0 && i3==0 && i4!=0 
                    {
                        collection.push(Vec4Int::new( i1, i2, i3, i4));
                        collection.push(Vec4Int::new( i1, i2, i3,-i4));
                    }
                    else if i1!=0 && i2!=0 && i3==0 && i4==0 
                    {
                        collection.push(Vec4Int::new( i1, i2, i3, i4));
                        collection.push(Vec4Int::new(-i1, i2, i3, i4));
                        collection.push(Vec4Int::new( i1,-i2, i3, i4));
                        collection.push(Vec4Int::new(-i1,-i2, i3, i4));
                    }
                    else if i1!=0 && i2==0 && i3!=0 && i4==0 
                    {
                        collection.push(Vec4Int::new( i1, i2, i3, i4));
                        collection.push(Vec4Int::new(-i1, i2, i3, i4));
                        collection.push(Vec4Int::new( i1, i2,-i3, i4));
                        collection.push(Vec4Int::new(-i1, i2,-i3, i4));
                    }
                    else if i1!=0 && i2==0 && i3==0 && i4!=0 
                    {
                        collection.push(Vec4Int::new( i1, i2, i3, i4));
                        collection.push(Vec4Int::new(-i1, i2, i3, i4));
                        collection.push(Vec4Int::new( i1, i2, i3,-i4));
                        collection.push(Vec4Int::new(-i1, i2, i3,-i4));
                    }
                    else if i1==0 && i2!=0 && i3!=0 && i4==0 
                    {
                        collection.push(Vec4Int::new( i1, i2, i3, i4));
                        collection.push(Vec4Int::new( i1,-i2, i3, i4));
                        collection.push(Vec4Int::new( i1, i2,-i3, i4));
                        collection.push(Vec4Int::new( i1,-i2,-i3, i4));
                    }
                    else if i1==0 && i2!=0 && i3==0 && i4!=0 
                    {
                        collection.push(Vec4Int::new( i1, i2, i3, i4));
                        collection.push(Vec4Int::new( i1,-i2, i3, i4));
                        collection.push(Vec4Int::new( i1, i2, i3,-i4));
                        collection.push(Vec4Int::new( i1,-i2, i3,-i4));
                    }
                    else if i1==0 && i2==0 && i3!=0 && i4!=0 
                    {
                        collection.push(Vec4Int::new( i1, i2, i3, i4));
                        collection.push(Vec4Int::new( i1, i2,-i3, i4));
                        collection.push(Vec4Int::new( i1, i2, i3,-i4));
                        collection.push(Vec4Int::new( i1, i2,-i3,-i4));
                    }
                    else if i1!=0 && i2!=0 && i3!=0 && i4==0 
                    {
                        collection.push(Vec4Int::new( i1, i2, i3, i4));
                        collection.push(Vec4Int::new(-i1, i2, i3, i4));
                        collection.push(Vec4Int::new( i1,-i2, i3, i4));
                        collection.push(Vec4Int::new(-i1,-i2, i3, i4));
                        collection.push(Vec4Int::new( i1, i2,-i3, i4));
                        collection.push(Vec4Int::new(-i1, i2,-i3, i4));
                        collection.push(Vec4Int::new( i1,-i2,-i3, i4));
                        collection.push(Vec4Int::new(-i1,-i2,-i3, i4));
                    }
                    else if i1!=0 && i2!=0 && i3==0 && i4!=0 
                    {
                        collection.push(Vec4Int::new( i1, i2, i3, i4));
                        collection.push(Vec4Int::new(-i1, i2, i3, i4));
                        collection.push(Vec4Int::new( i1,-i2, i3, i4));
                        collection.push(Vec4Int::new(-i1,-i2, i3, i4));
                        collection.push(Vec4Int::new( i1, i2, i3,-i4));
                        collection.push(Vec4Int::new(-i1, i2, i3,-i4));
                        collection.push(Vec4Int::new( i1,-i2, i3,-i4));
                        collection.push(Vec4Int::new(-i1,-i2, i3,-i4));
                    }
                    else if i1!=0 && i2==0 && i3!=0 && i4!=0 
                    {
                        collection.push(Vec4Int::new( i1, i2, i3, i4));
                        collection.push(Vec4Int::new(-i1, i2, i3, i4));
                        collection.push(Vec4Int::new( i1, i2,-i3, i4));
                        collection.push(Vec4Int::new(-i1, i2,-i3, i4));
                        collection.push(Vec4Int::new( i1, i2, i3,-i4));
                        collection.push(Vec4Int::new(-i1, i2, i3,-i4));
                        collection.push(Vec4Int::new( i1, i2,-i3,-i4));
                        collection.push(Vec4Int::new(-i1, i2,-i3,-i4));
                    }
                    else if i1==0 && i2!=0 && i3!=0 && i4!=0 
                    {
                        collection.push(Vec4Int::new( i1, i2, i3, i4));
                        collection.push(Vec4Int::new( i1,-i2, i3, i4));
                        collection.push(Vec4Int::new( i1, i2,-i3, i4));
                        collection.push(Vec4Int::new( i1,-i2,-i3, i4));
                        collection.push(Vec4Int::new( i1, i2, i3,-i4));
                        collection.push(Vec4Int::new( i1,-i2, i3,-i4));
                        collection.push(Vec4Int::new( i1, i2,-i3,-i4));
                        collection.push(Vec4Int::new( i1,-i2,-i3,-i4));
                    }
                    else if i1!=0 && i2!=0 && i3!=0 && i4!=0 
                    {
                        collection.push(Vec4Int::new( i1, i2, i3, i4));
                        collection.push(Vec4Int::new(-i1, i2, i3, i4));
                        collection.push(Vec4Int::new( i1,-i2, i3, i4));
                        collection.push(Vec4Int::new(-i1,-i2, i3, i4));
                        collection.push(Vec4Int::new( i1, i2,-i3, i4));
                        collection.push(Vec4Int::new(-i1, i2,-i3, i4));
                        collection.push(Vec4Int::new( i1,-i2,-i3, i4));
                        collection.push(Vec4Int::new(-i1,-i2,-i3, i4));
                        collection.push(Vec4Int::new( i1, i2, i3,-i4));
                        collection.push(Vec4Int::new(-i1, i2, i3,-i4));
                        collection.push(Vec4Int::new( i1,-i2, i3,-i4));
                        collection.push(Vec4Int::new(-i1,-i2, i3,-i4));
                        collection.push(Vec4Int::new( i1, i2,-i3,-i4));
                        collection.push(Vec4Int::new(-i1, i2,-i3,-i4));
                        collection.push(Vec4Int::new( i1,-i2,-i3,-i4));
                        collection.push(Vec4Int::new(-i1,-i2,-i3,-i4));
                    }

                    // Check all the points in the collection and empty it
                    while let Some(last_point) = collection.pop() 
                    {

                        // println!(" --------- HELLO --------------");
                        let mut possible_answer =  consider(radius, direction_of_rotation, epsilon_a, exactlogdep, int_center, coordinate_basis ,last_point);

                        if possible_answer != None
                        {
                            return possible_answer;
                        }

                    }

                    ///////////// DEBUG ZONE /////////////////
                    //
                    
                    if collection.len()>0
                    {
                        panic!("Something is wrong");
                    }

                    ///////////// END OF DEBUG ZONE /////////////


                    // Iterate
                    i4 = i4 +1
                }

                // Iterate
                i3 = i3+1;

            }

            //Iterate
            i2 = i2+1;
        }

        //Iterate 
        i1= i1+1;
    }

    return None;

    // return collection;
}

pub fn pseudo_nearest_neighbour_integer_coordinates(lattice_matrix_orthognal_part_inverse : Mat4, input_vector: Vec4 )  -> Vec4Int
{
    let floating_vector = lattice_matrix_orthognal_part_inverse*input_vector;
    let output = Vec4Int::new(  
        floating_vector[0].round() as Int,
        floating_vector[1].round() as Int,
        floating_vector[2].round() as Int,
        floating_vector[3].round() as Int,
        );
    return output;
}



// This is an implementation of Proposition 5.22
pub fn grid_problem_given_depth( direction: Comp, epsilon_a: Float,  exactlogdep: LogDepInt )-> ()
{

    // Trying to find the lattice points in the intersection of a disc and a half plane
    // Somtimes this region is called the "miniscule" region
    // Zoom out for the oversimilified picture
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
    //
    // The actual problem is a four dimensional lattice intersecting with a region like this

    // Bound the region in an ellipse
    let (center,mat,radius) =   ellipse_parameters_for_region_a(direction, epsilon_a);

    let centervec = Vec4::new(center.re,center.im,0.0,0.0);

    //create a 4 dimensional ellipsoid binding the region above combined with unit disc 

    let c = 1.0/(1.0+radius);

    let mat_big = 
        Mat4::new(mat[(0,0)]*c, mat[(0,1)]*c, 0.0  , 0.0,
        mat[(1,0)]*c, mat[(1,1)]*c, 0.0  , 0.0,
        0.0 ,         0.0         , 1.0-c, 0.0,
        0.0 ,         0.0         , 0.0  , 1.0-c);

    let radius_big = 1.0;

    let mat_lattice = 
        Mat4::new( 1.0 ,   SQRT2  , 0.0  , 0.0,
                   0.0 ,      0.0 , 1.0  ,  SQRT2 ,
                   1.0 ,  -SQRT2  , 0.0  , 0.0,
                   0.0 ,      0.0 , 1.0  , -SQRT2 );

    let mut reducable = mat_big*mat_lattice;

    let reduced = call_lll_on_nalgebra_matrix(reducable);


    let center4d = reduced.try_inverse().unwrap()*Vec4::new(center.re,center.im,0.0,0.0);
    let radius_big = 1.0;



    let operatornorm = find_operator_norm(reduced);
    let reduced_orthogonal_part_inverse = reduced.qr().q().try_inverse().unwrap();



    // Instead of this, we could have a rust thread
    // stream out lattice points by communicating messages.
    // and in parallel another one testing that it works
    let integer_center = pseudo_nearest_neighbour_integer_coordinates(reduced_orthogonal_part_inverse, centervec*SQRT2.pow(exactlogdep));
    test_integer_points_in_ball_around_integer_center_of_radius(operatornorm*SQRT2.pow(exactlogdep),  direction,epsilon_a, reducable , exactlogdep, integer_center);


    //////////////////// DEBUG ZONE  /////////////////////
    //
    //
    ///////////////////END OF DEBUG ZONE /////////////////

}



pub fn grid_problem( direction: Comp, epsilon_a: Float)-> ()
{
    for i in 0..9
    {
        grid_problem_given_depth(direction, epsilon_a, i);
    }
}

