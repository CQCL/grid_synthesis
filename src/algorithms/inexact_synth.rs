// Will implement the Ross Selinger algorithm here
// This is the content of 1403.2975v3 Section 7.3 on page 14


use crate::structs::rings::Float;
use crate::structs::rings::Int;
use crate::structs::rings::local_ring::Local;
use crate::structs::rings::zroot2::Zroot2;
use crate::structs::rings::LogDepInt;
use crate::structs::rings::Localizable;
use crate::structs::unimat::ExactUniMat;


use crate::structs::rings::special_values::sqrt2loc;


use crate::algorithms::local_prime_factorization::prime_factorization_of_loc;
use crate::algorithms::local_prime_factorization::attempt_to_write_this_number_as_sum_of_two_squares_in_loc;
use crate::algorithms::lll::nearest_plane;
use crate::algorithms::lll::lll_reduce;
use crate::algorithms::lll::gram_schmidt_orthogonalization;
use crate::algorithms::exact_synth::exact_synth;

use num_traits::Pow;
use num_traits::pow;
use num_traits::One;
use nalgebra::linalg::QR;


// Random number generators
use rand::thread_rng;
use rand::Rng;

// This will be rounded down to 64 or 128 bits
// For bigger float, add more digits please
pub const SQRT2:Float = 1.414213562373095048801688724209698078569671875376948073176679737990732478462;

type Comp = num_complex::Complex<Float>;
type Loc = Local::<Zroot2>;
type CompLoc = num_complex::Complex<Loc>;
type Mat2 = nalgebra::Matrix2<Float>;
type Mat4 = nalgebra::Matrix4<Float>;
type Mat4Int = nalgebra::Matrix4<Int>;
type Vec4 = nalgebra::Matrix4x1<Float>;
type Vec4Int = nalgebra::Matrix4x1<Int>;
use crate::structs::sunimat::SUniMat;
type Sunimat = SUniMat<Float>;
type Sunimatloc = SUniMat<Loc>;

type GridParams = (Comp, Float);

const integral_to_complex_mat: Mat4  = 
Mat4::new( 1.0 ,   SQRT2  , 0.0  , 0.0,
    0.0 ,      0.0 , 1.0  ,  SQRT2 ,
    1.0 ,  -SQRT2  , 0.0  , 0.0,
    0.0 ,      0.0 , 1.0  , -SQRT2 );

const complex_to_integral_mat: Mat4 = 
Mat4::new( 0.5      ,      0.0 , 0.5      ,       0.0 ,
    0.5/SQRT2,      0.0 ,-0.5/SQRT2,       0.0 ,
    0.0      ,      0.5 , 0.0      ,       0.5 ,
    0.0      ,0.5/SQRT2 , 0.0      ,-0.5/SQRT2 );


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
    // so membership condition of the ellipse is 
    // [ x y ] [ c+e_1^2(1-c)     e_1e_2(1-c)    ] [ x ]     <   2*c^2 - c^3
    //         [  e_1e_2(1-c)     c+e_2^2(1-c)   ] [ y ] 


    let c = epsilon_a;
    let sqrtc = c.sqrt();
    let e_1 = direction.re;
    let e_2 = direction.im;
    // todo!();

    let mut outputmatrix =  nalgebra::matrix![
        c+e_1*e_1*(1.0-c)+sqrtc   , e_1*e_2*(1.0-c) ;
    e_1*e_2*(1.0-c)  , c+e_2*e_2*(1.0-c)+sqrtc  ;
    ];

    outputmatrix *= 1.0/(1.0+sqrtc);

    return (direction*(1.0-c), outputmatrix, c*c*(2.0-c));

}






pub fn find_operator_norm(input : Mat4) -> Float
{
    return input.singular_values()[0];
}


pub fn extract_gate_coordinate_in_local_ring(integer_coords: Vec4Int) -> Option::<(Loc,Loc)>
{
    // println!("attempting_to_figure_out with {}", integer_coords );

    let zrt_left = Zroot2(integer_coords[0],integer_coords[1]);
    let zrt_right = Zroot2(integer_coords[2],integer_coords[3]);

    // We throw away the points if they were covered while checking smaller values of exactlogdep
    if zrt_left.is_divisible() && zrt_right.is_divisible()
    {
        return None;
    }
    else 
    {
        return Some( (Loc::from_base(zrt_left), Loc::from_base(zrt_right)) );
    }
}



pub fn get_comp_point_from_integer_coord( this_point: Vec4Int, exactlogdep : LogDepInt) -> (Comp, Comp)
{

    // println!("{}", this_point);
    // println!("{}", coordinate_basis);

    // take the input and mutiply coordinate change matrix
    // and scale the vector by sqrt2^exactlogdep
    // This will give the actual point in 4d-space
    let actual_point = vec4int_to_vec4float(this_point)*SQRT2.pow(-exactlogdep);

    // println!("{}", actual_point);

    let complex_point = Comp
    { 
        re: actual_point[0] + actual_point[1]*SQRT2,
        im: actual_point[2] + actual_point[3]*SQRT2,
    };

    // Is this correct?
    // Depends on the function 
    let complex_point_dot_conj = Comp
    { 
        re: actual_point[0] - actual_point[1]*SQRT2,
        im: actual_point[2] - actual_point[3]*SQRT2,
    };

    return (complex_point, complex_point_dot_conj);

}


// this tests if this_point is in the correct vicinity
// That is, if it roughly lies in the region close to the required place,,
// Here we make sure that it is truly within the region we want
// and only then we take the trouble of all the prime factorization
pub fn test_this_complex_pair_of_points( complex_point: Comp , complex_point_dot_conj: Comp, (direction, epsilon_a): GridParams) -> bool
{

    if complex_point.norm_sqr() <= 1.0
    {
        if complex_point_dot_conj.norm_sqr() <= 1.0
        {
            // println!("Norm squared in control ");
            if complex_point.re*direction.re+complex_point.im*direction.im > 1.0- epsilon_a
            {
                // println!("Test passed");
                return true;
            }
            else 
            {
                // println!("Half plane is not crossed");
                return false;
            }
        }
        else 
        {
            // println!("Second complex has norm bigger than 1");
            return false;
        }
    }
    else
    {
        // println!("First complex has norm bigger than 1");
        return false;
    }


}



pub fn make_exact_gate_from(left: Loc ,right: Loc, left_scaled: Loc, right_scaled: Loc) -> ExactUniMat
{
    // println!("sum_of_all_squares = {}", left*left + right*right +left_scaled*left_scaled + right_scaled*right_scaled );
    return ExactUniMat
    {
        mat: SUniMat{
            u: CompLoc{re: left_scaled, im: right_scaled},
            t: CompLoc{re: left, im:  right}
        },
        omega_exp : 0
    };

}


pub fn is_doubly_positive( num : Loc) -> bool
{
    let x = num.num;

    let x0 = x.0 as Float;
    let x1 = x.1 as Float;

    let r1 = x0 + SQRT2 * x1;
    let r2 = x0 - SQRT2 * x1;

    if r1 >= 0.0 && r2 >= 0.0
    {
        return true;
    }
    else 
    {
        return false
    }


}

// collection.push(Vec4Int::new( i1, i2, i3, i4));
// Return a complete unitary, if it can be returned
// Else return none
pub fn attempt_to_figure_out_gate_given_that_this_point_is_feasible( this_point: Vec4Int, exactlogdep: LogDepInt) -> Option::<ExactUniMat>
{

    let get_coord =  extract_gate_coordinate_in_local_ring(this_point);
    if get_coord == None
    {
        return None;
    }
    else
    {
        let (left_raw,right_raw) = get_coord.unwrap();


        let left_scaled = left_raw * pow(Loc::one() / sqrt2loc(), exactlogdep.try_into().unwrap() );
        let right_scaled = right_raw * pow(Loc::one() / sqrt2loc(), exactlogdep.try_into().unwrap() );

        let our_num = Loc::one() - left_scaled*left_scaled - right_scaled*right_scaled;


        // preliminary test before sending off to computationally expensive prime numbers
        // Based on Lemma 6.1 of 1403.2975
        if (! is_doubly_positive( our_num )) 
        {
            return None;
        }

        let sum_of_square = attempt_to_write_this_number_as_sum_of_two_squares_in_loc(our_num);

        if sum_of_square != None
        {
            let (left,right) = sum_of_square.unwrap();

            // DEBUG ZONE
            if ( left*left+right*right != our_num )
            {
                println!("left = {}", left);
                println!("right = {}", right);
                println!("our_num = {}", our_num);
                panic!("Sum of squares are wrong");
            }
            assert!( left*left+right*right + left_scaled * left_scaled + right_scaled * right_scaled == Loc::one() );
            // END OF DEBUG ZONE

            return Some(make_exact_gate_from(left,right,left_scaled,right_scaled));
        }
        else
        {   
            // Nothing to check
            return None;
        }

    }


}


pub fn consider( this_point: Vec4Int, exactlogdep: LogDepInt , ( direction_of_rotation, epsilon_a): GridParams ) -> Option::<ExactUniMat>
{

    let (complex_point, complex_point_dot_conj) = get_comp_point_from_integer_coord(this_point, exactlogdep);



    if test_this_complex_pair_of_points( complex_point, complex_point_dot_conj,( direction_of_rotation, epsilon_a ) )
    {

        // println!("Working Complex point is : {}", complex_point);
        return attempt_to_figure_out_gate_given_that_this_point_is_feasible( this_point , exactlogdep);
    }
    else
    {
        return None;
    }
}



pub fn mat4int_to_mat4(input: Mat4Int) -> Mat4
{
    return Mat4::new(
        input[(0,0)] as Float, input[(0,1)] as Float, input[(0,2)] as Float, input[(0,3)] as Float,
        input[(1,0)] as Float, input[(1,1)] as Float, input[(1,2)] as Float, input[(1,3)] as Float,
        input[(2,0)] as Float, input[(2,1)] as Float, input[(2,2)] as Float, input[(2,3)] as Float,
        input[(3,0)] as Float, input[(3,1)] as Float, input[(3,2)] as Float, input[(3,3)] as Float,
    );

}


pub fn round_mat4_to_int(input: Mat4) -> Mat4Int
{
    return Mat4Int::new(
        input[(0,0)].round() as Int, input[(0,1)].round() as Int, input[(0,2)].round() as Int, input[(0,3)].round() as Int,
        input[(1,0)].round() as Int, input[(1,1)].round() as Int, input[(1,2)].round() as Int, input[(1,3)].round() as Int,
        input[(2,0)].round() as Int, input[(2,1)].round() as Int, input[(2,2)].round() as Int, input[(2,3)].round() as Int,
        input[(3,0)].round() as Int, input[(3,1)].round() as Int, input[(3,2)].round() as Int, input[(3,3)].round() as Int,
    );

}

pub fn vec4float_to_vec4int(input: Vec4) -> Vec4Int
{
    return Vec4Int::new(input[0].round() as Int, input[1].round() as Int, input[2].round() as Int, input[3].round() as Int );
}

pub fn vec4int_to_vec4float(input: Vec4Int) -> Vec4
{
    return Vec4::new(input[0] as Float, input[1] as Float, input[2] as Float, input[3] as Float );
}



// Takes in the grid problem parameters and returns ellipsoid paramters of a 4d ellipsoid function
// This is a bit like ellipse_parameters_for_region_a, but it is 4 dimensional instead of two
// dimensional
pub fn generate_coordinates_and_center(direction : Comp , epsilon_a: Float) -> (Vec4, Mat4, Float)
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


    // First we will bound the region in an ellipse
    let (center,mat,radius) =   ellipse_parameters_for_region_a(direction, epsilon_a);


    //create a 4 dimensional ellipsoid binding the region above combined with unit disc \
    // The way to do this is through convex combination of ellipsoids
    //  If   [ x_1  x_2 ] [  a_11  a_12 ]  [  x_1  ]      <   r_1
    //                    [  a_21  a_22 ]  [  x_2  ]
    //
    // and   [ y_1  y_2 ] [  b_11  b_12 ]  [  y_1  ]      <   r_2
    //                    [  b_21  b_22 ]  [  y_2  ]
    //
    //  then for any real 0 < c < 1 we have 
    //
    //    [  x_1  x_2  y_1  y_2 ] [  c*a_11   c*a_12                           ] [ x_1 ]
    //                            [  c*a_21   c*a_22                           ] [ x_2 ]   < c*r_1  + (1-c)*r_2
    //                            [                   (1-c)*b_11   (1-c)*b_12  ] [ y_1 ]
    //                            [                   (1-c)*b_21   (1-c)*b_22  ] [ y_2 ] 
    //
    // Hence the cartesian product of these two dimensional ellipsoid is contained in the 4
    // dimensional convex-combination-like ellipsoid
    //
    //
    // Here we will set c=0.5, since this will minimize the ellipsoid volume
    //
    // However, the 4x4 matrix that you see above will be = 
    //
    //         (mat_big*mat_lattice).transpose()*(mat_big*mat_lattice)
    //
    // So setting below c = 1/SQRT2 and 1-c = 1/SQRT2 is actually the correct choice, because c gets
    // squared when the above multiplication happens.
    //
    // This was a pesky bug that I took some time to catch.
    //

    // convex combination parameter
    // the c above is sqrt(c_actual)
    let c_actual: Float = 0.5;

    let c = c_actual.sqrt();
    let one_minus_c = (1.0-c_actual).sqrt();

    let reducable = 
        Mat4::new(mat[(0,0)]*c, mat[(0,1)]*c  ,  0.0   , 0.0    ,
        mat[(1,0)]*c, mat[(1,1)]*c  ,  0.0   , 0.0    ,
        0.0     ,         0.0   , one_minus_c  , 0.0    ,
        0.0     ,         0.0   , 0.0    , one_minus_c  );


    let centervec = Vec4::new(center.re,center.im,0.0,0.0);

    let radius_big = c_actual*radius + (1.0-c_actual)*1.0;


    return (centervec,reducable, radius_big);

}




// This is an implementation of Proposition 5.22
// as in it does the job of Proposition 5.22
// the way it does it is mostly Nihar's invention
pub fn grid_problem_given_depth( exactlogdep: LogDepInt, (direction,epsilon) :GridParams ) -> Option::<ExactUniMat>
{
    // Second attempt to write a function based on LLL

    // first obtain the 4d_ellipse_matrix
    let ( ellipse_complex_coord_center, comp_to_4d_matrix, ellipse_4d_radius_squared) = generate_coordinates_and_center( direction, epsilon);

    // inflate the ellipse_4d_matrix
    // This takes integers directly to 4d-space without any SQRT2 multiplications
    let int_to_4d_space = comp_to_4d_matrix * integral_to_complex_mat * SQRT2.pow(-exactlogdep);

    // lll-reduce the int_to_4d_space_matrix
    // This creates a new basis for our integer lattice, and also carries the lattice change of
    // coordinates needed
    // Multiplying the two matrices here should give the int_to_4d_space
    let ( new_int_to_4d_space, standard_int_to_new_int ) = lll_reduce( int_to_4d_space );
    let new_int_to_standard_int = mat4int_inverse( standard_int_to_new_int );

    // find a center for the doing lattice search in new_int_coorinates
    let center_in_4d_space =  comp_to_4d_matrix * ellipse_complex_coord_center ;
    let new_int_to_4d_space_star = gram_schmidt_orthogonalization(new_int_to_4d_space);
    let (new_center_in_4d_space ,new_center_in_new_int ) = nearest_plane( new_int_to_4d_space, new_int_to_4d_space_star, center_in_4d_space );

    // update the radius because the centered in now changed
    let error_in_approximation_in_4d_space  = (new_center_in_4d_space - center_in_4d_space ).norm();
    let new_radius = ellipse_4d_radius_squared + error_in_approximation_in_4d_space * error_in_approximation_in_4d_space + 2.0 *error_in_approximation_in_4d_space * ellipse_4d_radius_squared.sqrt();



    //    // ----------------------- DEBUG ZONE 
    // Time to test this much
    // println!("::: Error_in_approx {}", error_in_approximation_in_4d_space );
    // let region_a_point = crate::tests::inexact_synth_tests::produce_random_point_inside_miniscule( direction, epsilon);
    // let mut rng = thread_rng();
    // let random_length : Float  = rng.gen_range(0.0..1.0);
    // let (x,y) = crate::tests::inexact_synth_tests::random_points_on_2d_circle();
    // let region_b_point = Comp{re:x*random_length,im:y*random_length};
    // assert!( test_this_complex_pair_of_points(region_a_point, region_b_point, ( direction, epsilon) )  );

    // // this random test_point should be in the ellipse
    // let new_center_comp =  integral_to_complex_mat * vec4int_to_vec4float( new_int_to_standard_int * new_center_in_new_int)*SQRT2.pow(-exactlogdep);
    // let new_complex_point = Comp
    // { 
    //     re: new_center_comp[0],
    //     im: new_center_comp[1],
    // };
    // println!("new_center is at {}",new_complex_point );
    // println!("This should be close to {}",Comp{re: ellipse_complex_coord_center[0],im: ellipse_complex_coord_center[1]} );
    //    let testvec =  Vec4::new(region_a_point.re,region_a_point.im,region_b_point.re,region_b_point.im);
    //    let offset = testvec - new_center_comp;
    //    let offset_skewed = comp_to_4d_matrix*offset;
    //    let skewed_distance = ( offset_skewed.transpose() * offset_skewed )[(0,0)];
    //    if ( skewed_distance > new_radius )
    //    {
    //        panic!("FIX THIS BUG");
    //    };
    //    // ------------------- END OF DEBUG ZONE



    //    // Instead of this, we could have a rust thread
    //    // stream out lattice points by communicating messages.
    //    // and in parallel another one testing that it works
    let possible_output =  test_integer_points_in_ball_around_integer_center_of_radius(new_radius, new_center_in_new_int, new_int_to_standard_int, exactlogdep,  (direction,epsilon));
    return possible_output;

}




pub fn grid_problem( direction: Comp, epsilon_a: Float)-> ExactUniMat
{
    let problem_parameters = ( direction, epsilon_a*epsilon_a/2.0);

    let mut answer: Option::<ExactUniMat>;

    let maxdepth = 60;

    for i in 0..maxdepth
    {

        // println!("--------------- \n \n ");
        // println!("Working with depth {}",i );
        let answer = grid_problem_given_depth(i, problem_parameters);
        if answer != None
        {
            let gate = answer.unwrap();
            
            println!("Found a candidate at depth {}",i );
            println!("which is \n {}", gate);

            return(gate);
        }
        // println!("--------------- \n \n ");
    }

    panic!("Nothing found under depth {}. Consider increasing maxdepth", maxdepth);
}


pub fn grid_problem_given_theta_and_epsilon( theta: Float, epsilon : Float ) -> ExactUniMat
{
    let x = Comp::new(0.0, theta);
    return grid_problem( x.exp() , epsilon);

}




// To be replaced by some multithreading implimentation
pub fn test_integer_points_in_ball_around_integer_center_of_radius(radius: Float, int_center: Vec4Int,  lattice_automorphism: Mat4Int,  exactlogdep: LogDepInt,  (direction_of_rotation, epsilon_a) : GridParams )  -> Option<ExactUniMat>
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

                        let int_standard = lattice_automorphism * (int_center + last_point);

                        let mut possible_answer =  consider( int_standard, exactlogdep,  (direction_of_rotation, epsilon_a));

                        if possible_answer != None
                        {
                            return possible_answer;
                        }

                    }


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

// This piece of code was generated by sagemath using matrices over symbolic ring
// It is basically the formula of the adjugate of a matrix
// What I really want to take matrix inverse of 4x4 matrix without losing integer precision
//
// I should find a better way to do this eventually
pub fn mat4int_inverse(x : Mat4Int) -> Mat4Int
{
    // When life gives you lemons, make them lemonade
    let det = x[(0,3)]*x[(1,2)]*x[(2,1)]*x[(3,0)] - x[(0,2)]*x[(1,3)]*x[(2,1)]*x[(3,0)] - x[(0,3)]*x[(1,1)]*x[(2,2)]*x[(3,0)] + x[(0,1)]*x[(1,3)]*x[(2,2)]*x[(3,0)] + x[(0,2)]*x[(1,1)]*x[(2,3)]*x[(3,0)] - x[(0,1)]*x[(1,2)]*x[(2,3)]*x[(3,0)] - x[(0,3)]*x[(1,2)]*x[(2,0)]*x[(3,1)] + x[(0,2)]*x[(1,3)]*x[(2,0)]*x[(3,1)] + x[(0,3)]*x[(1,0)]*x[(2,2)]*x[(3,1)] - x[(0,0)]*x[(1,3)]*x[(2,2)]*x[(3,1)] - x[(0,2)]*x[(1,0)]*x[(2,3)]*x[(3,1)] + x[(0,0)]*x[(1,2)]*x[(2,3)]*x[(3,1)] + x[(0,3)]*x[(1,1)]*x[(2,0)]*x[(3,2)] - x[(0,1)]*x[(1,3)]*x[(2,0)]*x[(3,2)] - x[(0,3)]*x[(1,0)]*x[(2,1)]*x[(3,2)] + x[(0,0)]*x[(1,3)]*x[(2,1)]*x[(3,2)] + x[(0,1)]*x[(1,0)]*x[(2,3)]*x[(3,2)] - x[(0,0)]*x[(1,1)]*x[(2,3)]*x[(3,2)] - x[(0,2)]*x[(1,1)]*x[(2,0)]*x[(3,3)] + x[(0,1)]*x[(1,2)]*x[(2,0)]*x[(3,3)] + x[(0,2)]*x[(1,0)]*x[(2,1)]*x[(3,3)] - x[(0,0)]*x[(1,2)]*x[(2,1)]*x[(3,3)] - x[(0,1)]*x[(1,0)]*x[(2,2)]*x[(3,3)] + x[(0,0)]*x[(1,1)]*x[(2,2)]*x[(3,3)];
    if det!=1 && det!=-1
    {
        panic!("Inverting an integer matrix that has no integer inverse");
    }


    let mut y = Mat4Int::zeros();

    y[( 0 , 0 )]= -x[(0,2)]*x[(1,1)]*x[(2,0)] + x[(0,1)]*x[(1,2)]*x[(2,0)] + x[(0,2)]*x[(1,0)]*x[(2,1)] - x[(0,0)]*x[(1,2)]*x[(2,1)] - x[(0,1)]*x[(1,0)]*x[(2,2)] + x[(0,0)]*x[(1,1)]*x[(2,2)] + x[(0,3)]*x[(1,0)]*x[(3,1)] - x[(0,0)]*x[(1,3)]*x[(3,1)] - x[(1,3)]*x[(2,2)]*x[(3,1)] + x[(1,2)]*x[(2,3)]*x[(3,1)] + x[(0,3)]*x[(2,0)]*x[(3,2)] + x[(1,3)]*x[(2,1)]*x[(3,2)] - x[(0,0)]*x[(2,3)]*x[(3,2)] - x[(1,1)]*x[(2,3)]*x[(3,2)] - x[(0,1)]*x[(1,0)]*x[(3,3)] + x[(0,0)]*x[(1,1)]*x[(3,3)] - x[(0,2)]*x[(2,0)]*x[(3,3)] - x[(1,2)]*x[(2,1)]*x[(3,3)] + x[(0,0)]*x[(2,2)]*x[(3,3)] + x[(1,1)]*x[(2,2)]*x[(3,3)] - x[(0,3)]*x[(1,1)]*x[(3,0)] + x[(0,1)]*x[(1,3)]*x[(3,0)] - x[(0,3)]*x[(2,2)]*x[(3,0)] + x[(0,2)]*x[(2,3)]*x[(3,0)] + (x[(0,0)]*(x[(1,1)] + x[(2,2)] + x[(3,3)]) - x[(0,0)]*x[(1,1)] + x[(1,2)]*x[(2,1)] - x[(0,0)]*x[(2,2)] - x[(1,1)]*x[(2,2)] + x[(1,3)]*x[(3,1)] + x[(2,3)]*x[(3,2)] - x[(0,0)]*x[(3,3)] - x[(1,1)]*x[(3,3)] - x[(2,2)]*x[(3,3)])*x[(0,0)] + (x[(0,1)]*(x[(1,1)] + x[(2,2)] + x[(3,3)]) - x[(0,1)]*x[(1,1)] - x[(0,2)]*x[(2,1)] - x[(0,3)]*x[(3,1)])*x[(1,0)] + (x[(0,2)]*(x[(1,1)] + x[(2,2)] + x[(3,3)]) - x[(0,1)]*x[(1,2)] - x[(0,2)]*x[(2,2)] - x[(0,3)]*x[(3,2)])*x[(2,0)] + (x[(0,3)]*(x[(1,1)] + x[(2,2)] + x[(3,3)]) - x[(0,1)]*x[(1,3)] - x[(0,2)]*x[(2,3)] - x[(0,3)]*x[(3,3)])*x[(3,0)] ;
    y[( 0 , 1 )]= (x[(0,0)]*(x[(1,1)] + x[(2,2)] + x[(3,3)]) - x[(0,0)]*x[(1,1)] + x[(1,2)]*x[(2,1)] - x[(0,0)]*x[(2,2)] - x[(1,1)]*x[(2,2)] + x[(1,3)]*x[(3,1)] + x[(2,3)]*x[(3,2)] - x[(0,0)]*x[(3,3)] - x[(1,1)]*x[(3,3)] - x[(2,2)]*x[(3,3)])*x[(0,1)] + (x[(0,1)]*(x[(1,1)] + x[(2,2)] + x[(3,3)]) - x[(0,1)]*x[(1,1)] - x[(0,2)]*x[(2,1)] - x[(0,3)]*x[(3,1)])*x[(1,1)] + (x[(0,2)]*(x[(1,1)] + x[(2,2)] + x[(3,3)]) - x[(0,1)]*x[(1,2)] - x[(0,2)]*x[(2,2)] - x[(0,3)]*x[(3,2)])*x[(2,1)] + (x[(0,3)]*(x[(1,1)] + x[(2,2)] + x[(3,3)]) - x[(0,1)]*x[(1,3)] - x[(0,2)]*x[(2,3)] - x[(0,3)]*x[(3,3)])*x[(3,1)] ;
    y[( 0 , 2 )]= (x[(0,0)]*(x[(1,1)] + x[(2,2)] + x[(3,3)]) - x[(0,0)]*x[(1,1)] + x[(1,2)]*x[(2,1)] - x[(0,0)]*x[(2,2)] - x[(1,1)]*x[(2,2)] + x[(1,3)]*x[(3,1)] + x[(2,3)]*x[(3,2)] - x[(0,0)]*x[(3,3)] - x[(1,1)]*x[(3,3)] - x[(2,2)]*x[(3,3)])*x[(0,2)] + (x[(0,1)]*(x[(1,1)] + x[(2,2)] + x[(3,3)]) - x[(0,1)]*x[(1,1)] - x[(0,2)]*x[(2,1)] - x[(0,3)]*x[(3,1)])*x[(1,2)] + (x[(0,2)]*(x[(1,1)] + x[(2,2)] + x[(3,3)]) - x[(0,1)]*x[(1,2)] - x[(0,2)]*x[(2,2)] - x[(0,3)]*x[(3,2)])*x[(2,2)] + (x[(0,3)]*(x[(1,1)] + x[(2,2)] + x[(3,3)]) - x[(0,1)]*x[(1,3)] - x[(0,2)]*x[(2,3)] - x[(0,3)]*x[(3,3)])*x[(3,2)] ;
    y[( 0 , 3 )]= (x[(0,0)]*(x[(1,1)] + x[(2,2)] + x[(3,3)]) - x[(0,0)]*x[(1,1)] + x[(1,2)]*x[(2,1)] - x[(0,0)]*x[(2,2)] - x[(1,1)]*x[(2,2)] + x[(1,3)]*x[(3,1)] + x[(2,3)]*x[(3,2)] - x[(0,0)]*x[(3,3)] - x[(1,1)]*x[(3,3)] - x[(2,2)]*x[(3,3)])*x[(0,3)] + (x[(0,1)]*(x[(1,1)] + x[(2,2)] + x[(3,3)]) - x[(0,1)]*x[(1,1)] - x[(0,2)]*x[(2,1)] - x[(0,3)]*x[(3,1)])*x[(1,3)] + (x[(0,2)]*(x[(1,1)] + x[(2,2)] + x[(3,3)]) - x[(0,1)]*x[(1,2)] - x[(0,2)]*x[(2,2)] - x[(0,3)]*x[(3,2)])*x[(2,3)] + (x[(0,3)]*(x[(1,1)] + x[(2,2)] + x[(3,3)]) - x[(0,1)]*x[(1,3)] - x[(0,2)]*x[(2,3)] - x[(0,3)]*x[(3,3)])*x[(3,3)] ;
    y[( 1 , 0 )]= ((x[(0,0)] + x[(2,2)] + x[(3,3)])*x[(1,0)] - x[(0,0)]*x[(1,0)] - x[(1,2)]*x[(2,0)] - x[(1,3)]*x[(3,0)])*x[(0,0)] + ((x[(0,0)] + x[(2,2)] + x[(3,3)])*x[(1,1)] - x[(0,0)]*x[(1,1)] + x[(0,2)]*x[(2,0)] - x[(0,0)]*x[(2,2)] - x[(1,1)]*x[(2,2)] + x[(2,3)]*x[(3,2)] - x[(0,0)]*x[(3,3)] - x[(1,1)]*x[(3,3)] - x[(2,2)]*x[(3,3)] + x[(0,3)]*x[(3,0)])*x[(1,0)] - (x[(0,2)]*x[(1,0)] - (x[(0,0)] + x[(2,2)] + x[(3,3)])*x[(1,2)] + x[(1,2)]*x[(2,2)] + x[(1,3)]*x[(3,2)])*x[(2,0)] - (x[(0,3)]*x[(1,0)] - (x[(0,0)] + x[(2,2)] + x[(3,3)])*x[(1,3)] + x[(1,2)]*x[(2,3)] + x[(1,3)]*x[(3,3)])*x[(3,0)] ;
    y[( 1 , 1 )]= -x[(0,2)]*x[(1,1)]*x[(2,0)] + x[(0,1)]*x[(1,2)]*x[(2,0)] + x[(0,2)]*x[(1,0)]*x[(2,1)] - x[(0,0)]*x[(1,2)]*x[(2,1)] - x[(0,1)]*x[(1,0)]*x[(2,2)] + x[(0,0)]*x[(1,1)]*x[(2,2)] + x[(0,3)]*x[(1,0)]*x[(3,1)] - x[(0,0)]*x[(1,3)]*x[(3,1)] - x[(1,3)]*x[(2,2)]*x[(3,1)] + x[(1,2)]*x[(2,3)]*x[(3,1)] + x[(0,3)]*x[(2,0)]*x[(3,2)] + x[(1,3)]*x[(2,1)]*x[(3,2)] - x[(0,0)]*x[(2,3)]*x[(3,2)] - x[(1,1)]*x[(2,3)]*x[(3,2)] - x[(0,1)]*x[(1,0)]*x[(3,3)] + x[(0,0)]*x[(1,1)]*x[(3,3)] - x[(0,2)]*x[(2,0)]*x[(3,3)] - x[(1,2)]*x[(2,1)]*x[(3,3)] + x[(0,0)]*x[(2,2)]*x[(3,3)] + x[(1,1)]*x[(2,2)]*x[(3,3)] - x[(0,3)]*x[(1,1)]*x[(3,0)] + x[(0,1)]*x[(1,3)]*x[(3,0)] - x[(0,3)]*x[(2,2)]*x[(3,0)] + x[(0,2)]*x[(2,3)]*x[(3,0)] + ((x[(0,0)] + x[(2,2)] + x[(3,3)])*x[(1,0)] - x[(0,0)]*x[(1,0)] - x[(1,2)]*x[(2,0)] - x[(1,3)]*x[(3,0)])*x[(0,1)] + ((x[(0,0)] + x[(2,2)] + x[(3,3)])*x[(1,1)] - x[(0,0)]*x[(1,1)] + x[(0,2)]*x[(2,0)] - x[(0,0)]*x[(2,2)] - x[(1,1)]*x[(2,2)] + x[(2,3)]*x[(3,2)] - x[(0,0)]*x[(3,3)] - x[(1,1)]*x[(3,3)] - x[(2,2)]*x[(3,3)] + x[(0,3)]*x[(3,0)])*x[(1,1)] - (x[(0,2)]*x[(1,0)] - (x[(0,0)] + x[(2,2)] + x[(3,3)])*x[(1,2)] + x[(1,2)]*x[(2,2)] + x[(1,3)]*x[(3,2)])*x[(2,1)] - (x[(0,3)]*x[(1,0)] - (x[(0,0)] + x[(2,2)] + x[(3,3)])*x[(1,3)] + x[(1,2)]*x[(2,3)] + x[(1,3)]*x[(3,3)])*x[(3,1)] ;
    y[( 1 , 2 )]= ((x[(0,0)] + x[(2,2)] + x[(3,3)])*x[(1,0)] - x[(0,0)]*x[(1,0)] - x[(1,2)]*x[(2,0)] - x[(1,3)]*x[(3,0)])*x[(0,2)] + ((x[(0,0)] + x[(2,2)] + x[(3,3)])*x[(1,1)] - x[(0,0)]*x[(1,1)] + x[(0,2)]*x[(2,0)] - x[(0,0)]*x[(2,2)] - x[(1,1)]*x[(2,2)] + x[(2,3)]*x[(3,2)] - x[(0,0)]*x[(3,3)] - x[(1,1)]*x[(3,3)] - x[(2,2)]*x[(3,3)] + x[(0,3)]*x[(3,0)])*x[(1,2)] - (x[(0,2)]*x[(1,0)] - (x[(0,0)] + x[(2,2)] + x[(3,3)])*x[(1,2)] + x[(1,2)]*x[(2,2)] + x[(1,3)]*x[(3,2)])*x[(2,2)] - (x[(0,3)]*x[(1,0)] - (x[(0,0)] + x[(2,2)] + x[(3,3)])*x[(1,3)] + x[(1,2)]*x[(2,3)] + x[(1,3)]*x[(3,3)])*x[(3,2)] ;
    y[( 1 , 3 )]= ((x[(0,0)] + x[(2,2)] + x[(3,3)])*x[(1,0)] - x[(0,0)]*x[(1,0)] - x[(1,2)]*x[(2,0)] - x[(1,3)]*x[(3,0)])*x[(0,3)] + ((x[(0,0)] + x[(2,2)] + x[(3,3)])*x[(1,1)] - x[(0,0)]*x[(1,1)] + x[(0,2)]*x[(2,0)] - x[(0,0)]*x[(2,2)] - x[(1,1)]*x[(2,2)] + x[(2,3)]*x[(3,2)] - x[(0,0)]*x[(3,3)] - x[(1,1)]*x[(3,3)] - x[(2,2)]*x[(3,3)] + x[(0,3)]*x[(3,0)])*x[(1,3)] - (x[(0,2)]*x[(1,0)] - (x[(0,0)] + x[(2,2)] + x[(3,3)])*x[(1,2)] + x[(1,2)]*x[(2,2)] + x[(1,3)]*x[(3,2)])*x[(2,3)] - (x[(0,3)]*x[(1,0)] - (x[(0,0)] + x[(2,2)] + x[(3,3)])*x[(1,3)] + x[(1,2)]*x[(2,3)] + x[(1,3)]*x[(3,3)])*x[(3,3)] ;
    y[( 2 , 0 )]= ((x[(0,0)] + x[(1,1)] + x[(3,3)])*x[(2,0)] - x[(0,0)]*x[(2,0)] - x[(1,0)]*x[(2,1)] - x[(2,3)]*x[(3,0)])*x[(0,0)] - (x[(0,1)]*x[(2,0)] - (x[(0,0)] + x[(1,1)] + x[(3,3)])*x[(2,1)] + x[(1,1)]*x[(2,1)] + x[(2,3)]*x[(3,1)])*x[(1,0)] + (x[(0,1)]*x[(1,0)] - x[(0,0)]*x[(1,1)] + (x[(0,0)] + x[(1,1)] + x[(3,3)])*x[(2,2)] - x[(0,0)]*x[(2,2)] - x[(1,1)]*x[(2,2)] + x[(1,3)]*x[(3,1)] - x[(0,0)]*x[(3,3)] - x[(1,1)]*x[(3,3)] - x[(2,2)]*x[(3,3)] + x[(0,3)]*x[(3,0)])*x[(2,0)] - (x[(0,3)]*x[(2,0)] + x[(1,3)]*x[(2,1)] - (x[(0,0)] + x[(1,1)] + x[(3,3)])*x[(2,3)] + x[(2,3)]*x[(3,3)])*x[(3,0)] ;
    y[( 2 , 1 )]= ((x[(0,0)] + x[(1,1)] + x[(3,3)])*x[(2,0)] - x[(0,0)]*x[(2,0)] - x[(1,0)]*x[(2,1)] - x[(2,3)]*x[(3,0)])*x[(0,1)] - (x[(0,1)]*x[(2,0)] - (x[(0,0)] + x[(1,1)] + x[(3,3)])*x[(2,1)] + x[(1,1)]*x[(2,1)] + x[(2,3)]*x[(3,1)])*x[(1,1)] + (x[(0,1)]*x[(1,0)] - x[(0,0)]*x[(1,1)] + (x[(0,0)] + x[(1,1)] + x[(3,3)])*x[(2,2)] - x[(0,0)]*x[(2,2)] - x[(1,1)]*x[(2,2)] + x[(1,3)]*x[(3,1)] - x[(0,0)]*x[(3,3)] - x[(1,1)]*x[(3,3)] - x[(2,2)]*x[(3,3)] + x[(0,3)]*x[(3,0)])*x[(2,1)] - (x[(0,3)]*x[(2,0)] + x[(1,3)]*x[(2,1)] - (x[(0,0)] + x[(1,1)] + x[(3,3)])*x[(2,3)] + x[(2,3)]*x[(3,3)])*x[(3,1)] ;
    y[( 2 , 2 )]= -x[(0,2)]*x[(1,1)]*x[(2,0)] + x[(0,1)]*x[(1,2)]*x[(2,0)] + x[(0,2)]*x[(1,0)]*x[(2,1)] - x[(0,0)]*x[(1,2)]*x[(2,1)] - x[(0,1)]*x[(1,0)]*x[(2,2)] + x[(0,0)]*x[(1,1)]*x[(2,2)] + x[(0,3)]*x[(1,0)]*x[(3,1)] - x[(0,0)]*x[(1,3)]*x[(3,1)] - x[(1,3)]*x[(2,2)]*x[(3,1)] + x[(1,2)]*x[(2,3)]*x[(3,1)] + x[(0,3)]*x[(2,0)]*x[(3,2)] + x[(1,3)]*x[(2,1)]*x[(3,2)] - x[(0,0)]*x[(2,3)]*x[(3,2)] - x[(1,1)]*x[(2,3)]*x[(3,2)] - x[(0,1)]*x[(1,0)]*x[(3,3)] + x[(0,0)]*x[(1,1)]*x[(3,3)] - x[(0,2)]*x[(2,0)]*x[(3,3)] - x[(1,2)]*x[(2,1)]*x[(3,3)] + x[(0,0)]*x[(2,2)]*x[(3,3)] + x[(1,1)]*x[(2,2)]*x[(3,3)] - x[(0,3)]*x[(1,1)]*x[(3,0)] + x[(0,1)]*x[(1,3)]*x[(3,0)] - x[(0,3)]*x[(2,2)]*x[(3,0)] + x[(0,2)]*x[(2,3)]*x[(3,0)] + ((x[(0,0)] + x[(1,1)] + x[(3,3)])*x[(2,0)] - x[(0,0)]*x[(2,0)] - x[(1,0)]*x[(2,1)] - x[(2,3)]*x[(3,0)])*x[(0,2)] - (x[(0,1)]*x[(2,0)] - (x[(0,0)] + x[(1,1)] + x[(3,3)])*x[(2,1)] + x[(1,1)]*x[(2,1)] + x[(2,3)]*x[(3,1)])*x[(1,2)] + (x[(0,1)]*x[(1,0)] - x[(0,0)]*x[(1,1)] + (x[(0,0)] + x[(1,1)] + x[(3,3)])*x[(2,2)] - x[(0,0)]*x[(2,2)] - x[(1,1)]*x[(2,2)] + x[(1,3)]*x[(3,1)] - x[(0,0)]*x[(3,3)] - x[(1,1)]*x[(3,3)] - x[(2,2)]*x[(3,3)] + x[(0,3)]*x[(3,0)])*x[(2,2)] - (x[(0,3)]*x[(2,0)] + x[(1,3)]*x[(2,1)] - (x[(0,0)] + x[(1,1)] + x[(3,3)])*x[(2,3)] + x[(2,3)]*x[(3,3)])*x[(3,2)] ;
    y[( 2 , 3 )]= ((x[(0,0)] + x[(1,1)] + x[(3,3)])*x[(2,0)] - x[(0,0)]*x[(2,0)] - x[(1,0)]*x[(2,1)] - x[(2,3)]*x[(3,0)])*x[(0,3)] - (x[(0,1)]*x[(2,0)] - (x[(0,0)] + x[(1,1)] + x[(3,3)])*x[(2,1)] + x[(1,1)]*x[(2,1)] + x[(2,3)]*x[(3,1)])*x[(1,3)] + (x[(0,1)]*x[(1,0)] - x[(0,0)]*x[(1,1)] + (x[(0,0)] + x[(1,1)] + x[(3,3)])*x[(2,2)] - x[(0,0)]*x[(2,2)] - x[(1,1)]*x[(2,2)] + x[(1,3)]*x[(3,1)] - x[(0,0)]*x[(3,3)] - x[(1,1)]*x[(3,3)] - x[(2,2)]*x[(3,3)] + x[(0,3)]*x[(3,0)])*x[(2,3)] - (x[(0,3)]*x[(2,0)] + x[(1,3)]*x[(2,1)] - (x[(0,0)] + x[(1,1)] + x[(3,3)])*x[(2,3)] + x[(2,3)]*x[(3,3)])*x[(3,3)] ;
    y[( 3 , 0 )]= -(x[(1,0)]*x[(3,1)] + x[(2,0)]*x[(3,2)] - (x[(0,0)] + x[(1,1)] + x[(2,2)])*x[(3,0)] + x[(0,0)]*x[(3,0)])*x[(0,0)] + ((x[(0,0)] + x[(1,1)] + x[(2,2)])*x[(3,1)] - x[(1,1)]*x[(3,1)] - x[(2,1)]*x[(3,2)] - x[(0,1)]*x[(3,0)])*x[(1,0)] - (x[(1,2)]*x[(3,1)] - (x[(0,0)] + x[(1,1)] + x[(2,2)])*x[(3,2)] + x[(2,2)]*x[(3,2)] + x[(0,2)]*x[(3,0)])*x[(2,0)] + (x[(0,1)]*x[(1,0)] - x[(0,0)]*x[(1,1)] + x[(0,2)]*x[(2,0)] + x[(1,2)]*x[(2,1)] - x[(0,0)]*x[(2,2)] - x[(1,1)]*x[(2,2)] + (x[(0,0)] + x[(1,1)] + x[(2,2)])*x[(3,3)] - x[(0,0)]*x[(3,3)] - x[(1,1)]*x[(3,3)] - x[(2,2)]*x[(3,3)])*x[(3,0)] ;
    y[( 3 , 1 )]= -(x[(1,0)]*x[(3,1)] + x[(2,0)]*x[(3,2)] - (x[(0,0)] + x[(1,1)] + x[(2,2)])*x[(3,0)] + x[(0,0)]*x[(3,0)])*x[(0,1)] + ((x[(0,0)] + x[(1,1)] + x[(2,2)])*x[(3,1)] - x[(1,1)]*x[(3,1)] - x[(2,1)]*x[(3,2)] - x[(0,1)]*x[(3,0)])*x[(1,1)] - (x[(1,2)]*x[(3,1)] - (x[(0,0)] + x[(1,1)] + x[(2,2)])*x[(3,2)] + x[(2,2)]*x[(3,2)] + x[(0,2)]*x[(3,0)])*x[(2,1)] + (x[(0,1)]*x[(1,0)] - x[(0,0)]*x[(1,1)] + x[(0,2)]*x[(2,0)] + x[(1,2)]*x[(2,1)] - x[(0,0)]*x[(2,2)] - x[(1,1)]*x[(2,2)] + (x[(0,0)] + x[(1,1)] + x[(2,2)])*x[(3,3)] - x[(0,0)]*x[(3,3)] - x[(1,1)]*x[(3,3)] - x[(2,2)]*x[(3,3)])*x[(3,1)] ;
    y[( 3 , 2 )]= -(x[(1,0)]*x[(3,1)] + x[(2,0)]*x[(3,2)] - (x[(0,0)] + x[(1,1)] + x[(2,2)])*x[(3,0)] + x[(0,0)]*x[(3,0)])*x[(0,2)] + ((x[(0,0)] + x[(1,1)] + x[(2,2)])*x[(3,1)] - x[(1,1)]*x[(3,1)] - x[(2,1)]*x[(3,2)] - x[(0,1)]*x[(3,0)])*x[(1,2)] - (x[(1,2)]*x[(3,1)] - (x[(0,0)] + x[(1,1)] + x[(2,2)])*x[(3,2)] + x[(2,2)]*x[(3,2)] + x[(0,2)]*x[(3,0)])*x[(2,2)] + (x[(0,1)]*x[(1,0)] - x[(0,0)]*x[(1,1)] + x[(0,2)]*x[(2,0)] + x[(1,2)]*x[(2,1)] - x[(0,0)]*x[(2,2)] - x[(1,1)]*x[(2,2)] + (x[(0,0)] + x[(1,1)] + x[(2,2)])*x[(3,3)] - x[(0,0)]*x[(3,3)] - x[(1,1)]*x[(3,3)] - x[(2,2)]*x[(3,3)])*x[(3,2)] ;
    y[( 3 , 3 )]= -x[(0,2)]*x[(1,1)]*x[(2,0)] + x[(0,1)]*x[(1,2)]*x[(2,0)] + x[(0,2)]*x[(1,0)]*x[(2,1)] - x[(0,0)]*x[(1,2)]*x[(2,1)] - x[(0,1)]*x[(1,0)]*x[(2,2)] + x[(0,0)]*x[(1,1)]*x[(2,2)] + x[(0,3)]*x[(1,0)]*x[(3,1)] - x[(0,0)]*x[(1,3)]*x[(3,1)] - x[(1,3)]*x[(2,2)]*x[(3,1)] + x[(1,2)]*x[(2,3)]*x[(3,1)] + x[(0,3)]*x[(2,0)]*x[(3,2)] + x[(1,3)]*x[(2,1)]*x[(3,2)] - x[(0,0)]*x[(2,3)]*x[(3,2)] - x[(1,1)]*x[(2,3)]*x[(3,2)] - x[(0,1)]*x[(1,0)]*x[(3,3)] + x[(0,0)]*x[(1,1)]*x[(3,3)] - x[(0,2)]*x[(2,0)]*x[(3,3)] - x[(1,2)]*x[(2,1)]*x[(3,3)] + x[(0,0)]*x[(2,2)]*x[(3,3)] + x[(1,1)]*x[(2,2)]*x[(3,3)] - x[(0,3)]*x[(1,1)]*x[(3,0)] + x[(0,1)]*x[(1,3)]*x[(3,0)] - x[(0,3)]*x[(2,2)]*x[(3,0)] + x[(0,2)]*x[(2,3)]*x[(3,0)] - (x[(1,0)]*x[(3,1)] + x[(2,0)]*x[(3,2)] - (x[(0,0)] + x[(1,1)] + x[(2,2)])*x[(3,0)] + x[(0,0)]*x[(3,0)])*x[(0,3)] + ((x[(0,0)] + x[(1,1)] + x[(2,2)])*x[(3,1)] - x[(1,1)]*x[(3,1)] - x[(2,1)]*x[(3,2)] - x[(0,1)]*x[(3,0)])*x[(1,3)] - (x[(1,2)]*x[(3,1)] - (x[(0,0)] + x[(1,1)] + x[(2,2)])*x[(3,2)] + x[(2,2)]*x[(3,2)] + x[(0,2)]*x[(3,0)])*x[(2,3)] + (x[(0,1)]*x[(1,0)] - x[(0,0)]*x[(1,1)] + x[(0,2)]*x[(2,0)] + x[(1,2)]*x[(2,1)] - x[(0,0)]*x[(2,2)] - x[(1,1)]*x[(2,2)] + (x[(0,0)] + x[(1,1)] + x[(2,2)])*x[(3,3)] - x[(0,0)]*x[(3,3)] - x[(1,1)]*x[(3,3)] - x[(2,2)]*x[(3,3)])*x[(3,3)] ;


    return y*det;
}

