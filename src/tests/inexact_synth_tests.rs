use crate::algorithms::inexact_synth::grid_problem;
use crate::algorithms::inexact_synth::grid_problem_given_depth;
use crate::algorithms::inexact_synth::get_comp_point_from_integer_coord;
use crate::algorithms::inexact_synth::ellipse_parameters_for_region_a;
use crate::algorithms::inexact_synth::test_this_complex_pair_of_points;
use crate::algorithms::inexact_synth::generate_coordinates_and_center;
use crate::algorithms::exact_synth::exact_synth;





use num_traits::One;


use crate::structs::rings::Float;
use crate::structs::rings::local_ring::Local;
use crate::structs::rings::zroot2::Zroot2;
use crate::structs::rings::Int;
use crate::structs::rings::LogDepInt;
use crate::structs::rings::Localizable;
use crate::structs::unimat::ExactUniMat;
type Comp = num_complex::Complex<Float>;
type Loc = Local::<Zroot2>;
type CompLoc = num_complex::Complex<Loc>;
type Mat2 = nalgebra::Matrix2<Float>;
type Mat4 = nalgebra::Matrix4<Float>;
type Vec4 = nalgebra::Matrix4x1<Float>;
type Vec4Int = nalgebra::Matrix4x1<Int>;
type GridParams = (Comp, Float);

// Value of SQRT2
const SQRT2:Float = 1.414213562373095048801688724209698078569671875376948073176679737990732478462;

// Random number generators
use rand::thread_rng;
use rand::Rng;

// error threshold for testing Float equality
const threshold: Float = 0.000000000001;


// function to test floating point equality of complex numbers
pub fn two_comp_equal_up_to_theshold(left: Comp, right: Comp) -> bool
{
    return (left-right).norm_sqr() < threshold;

}

// function to test floating point equality 
pub fn two_float_equal_up_to_threshold(left: Float, right: Float) -> bool
{
    return (left-right).abs() < threshold;
}



#[test]
pub fn comp_from_basis_and_vector_test()
{


    let this_point = Vec4Int::new(1,0,0,0);


    let (left,right) = get_comp_point_from_integer_coord(this_point, 0);

    assert!(two_comp_equal_up_to_theshold(left, Comp::one()) );
    assert!(two_comp_equal_up_to_theshold(right, Comp::one()) );

    let (left,right) = get_comp_point_from_integer_coord(this_point, 2);

    assert!(two_comp_equal_up_to_theshold(left, 0.5* Comp::one()) );
    assert!(two_comp_equal_up_to_theshold(right, 0.5 * Comp::one()) );

    let (left,right) = get_comp_point_from_integer_coord(this_point, -2);

    assert!(two_comp_equal_up_to_theshold(left, 2.0* Comp::one()) );
    assert!(two_comp_equal_up_to_theshold(right, 2.0 * Comp::one()) );
}


pub fn produce_random_point_inside_miniscule(direction: Comp, epsilon: Float) -> Comp
{

    let mut rng = thread_rng();
    let a1: Float = rng.gen_range((1.0-epsilon)..1.0);
    let pythogoras_vert = (1.0 - a1*a1).sqrt();
    assert!( two_float_equal_up_to_threshold( a1*a1 + pythogoras_vert*pythogoras_vert, 1.0));

    let perp_comp = direction*Comp{re:0.0,im:1.0};
    assert!( two_float_equal_up_to_threshold( perp_comp.re*direction.re + perp_comp.im*direction.im , 0.0));
    let a2: Float = rng.gen_range(-pythogoras_vert..pythogoras_vert);

    let test_point = direction*a1 + perp_comp*a2;

    return test_point;

}

// Not really uniform distribution, but its okay
pub fn random_points_on_2d_circle() -> (Float,Float)
{
    let mut rng = thread_rng();

    let x1 :Float = rng.gen_range(-1.0..1.0);
    let y1 :Float = rng.gen_range(-1.0..1.0);

    let norm = (x1*x1 + y1*y1).sqrt();

    let (x,y) = (x1/norm, y1/norm);

    return (x,y);
}

pub fn produce_random_grid_paramters() -> GridParams
{

    let mut rng = thread_rng();
    let (x,y) = random_points_on_2d_circle();

    let rand_comp = Comp{re:x,im:y};
    let epsilon: Float = rng.gen_range( 0.0000001..0.02 );
    return (rand_comp, epsilon);


}



pub fn testing_the_ellipse_paramters_randomly() 
{
    let (rand_comp,epsilon) = produce_random_grid_paramters();
    let (center,mat,radius) =    ellipse_parameters_for_region_a( rand_comp, epsilon);

    assert!(two_float_equal_up_to_threshold( center.re*rand_comp.re + center.im*rand_comp.im , 1.0-epsilon) );

    let test_point = produce_random_point_inside_miniscule(rand_comp, epsilon);

    // this random test_point should be in the ellipse
    let offset = test_point - center;
    let offset_vector  = nalgebra::matrix![ offset.re  , offset.im ];
    let offset_vector_trans = offset_vector.transpose();
    let offset_skewed = mat*offset_vector_trans;

    let skewed_distance = offset_skewed.transpose() * offset_skewed;

    assert!( test_this_complex_pair_of_points(test_point, test_point, ( rand_comp, epsilon) )  );
}


#[test]
pub fn rapidly_testing_ellipse_parameters() 
{
    for i in 1..1000
    {
        println!("Test number: {}", i);
        testing_the_ellipse_paramters_randomly();
    }

    println!("It worked 1000 times!");

}




pub fn testing_4d_ellipsoid_parameters_randomly()
{
    let (direction,epsilon) = produce_random_grid_paramters();

    let region_a_point = produce_random_point_inside_miniscule(direction,epsilon);

    let mut rng = thread_rng();
    let random_length : Float  = rng.gen_range(0.0..1.0);
    let (x,y) = random_points_on_2d_circle();
    let region_b_point = Comp{re:x*random_length,im:y*random_length};


    let (center, matrix,  radius) = generate_coordinates_and_center(direction,epsilon);


    // let reducable = reduced*lattice_autmorphism;

    assert!( test_this_complex_pair_of_points( region_a_point, region_b_point, (direction, epsilon) ));

    let testvec =  Vec4::new(region_a_point.re,region_a_point.im,region_b_point.re,region_b_point.im);
    let offset = testvec - center;
    let offset_skewed = matrix*offset;
    let skewed_distance = ( offset_skewed.transpose() * offset_skewed )[(0,0)];

    println!("Point's distance is {} % of radius",100.0*skewed_distance/radius );
    assert! ( skewed_distance <=  radius)

}

#[test]
pub fn rapidly_testing_4d_ellipsoid_parameters()
{

    for i in 0..5000
    {
        testing_4d_ellipsoid_parameters_randomly();
    }
    println!("Test worked 5000 times!");
}

#[test]
pub fn matrix_norm_is_what_i_expect()
{
    let test = Vec4::new(2.0,0.0,0.0,0.0);
    assert!( two_float_equal_up_to_threshold( test.norm() , 2.0 )  );
}






pub fn generate_random_4by4_matrix() -> Mat4
{

    let mut rng = thread_rng();

    let mut out = Mat4::zeros();
    for i in 0..4
    {
        for j in 0..4
        {
            out[(i,j)] = rng.gen_range(-10.0..10.0);
        }
    }


    return out;
}

#[test]
pub fn random_inexact_synth_test()
{
    let (direction,epsilon) = produce_random_grid_paramters();

    println!("(direction ,epsilon) = ({},{})",direction,epsilon );

    let answer = grid_problem(direction,epsilon);
    println!("answer =\n {}",answer );
    println!("as a float gate = \n{}",answer.to_float_gate_upto_t_count() );

    println!("Gate sequence {}", exact_synth(answer) );

}



#[test]
pub fn inexact_synth_testing() 
{
    // grid_problem_given_depth(0, (Comp::one(), 0.19));

    let answer = grid_problem(Comp::one(), 0.020);

    println!("{}", answer.mat);
    println!("{}", answer.omega_exp);

    assert_eq!( ExactUniMat::one() , answer);

}






