use crate::algorithms::inexact_synth::grid_problem;
use crate::algorithms::inexact_synth::get_comp_point_from_basis_and_vector;
use crate::algorithms::inexact_synth::ellipse_parameters_for_region_a;
use crate::algorithms::inexact_synth::test_this_complex_pair_of_points;


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
    let coordinate_basis = Mat4::new(1.0,0.0,0.0,0.0,
                                     0.0,1.0,0.0,0.0,
                                     0.0,0.0,1.0,0.0,
                                     0.0,0.0,0.0,1.0);

    
    let (left,right) = get_comp_point_from_basis_and_vector(this_point, coordinate_basis, 0);
    
    assert!(two_comp_equal_up_to_theshold(left, Comp::one()) );
    assert!(two_comp_equal_up_to_theshold(right, Comp::one()) );
             
    let (left,right) = get_comp_point_from_basis_and_vector(this_point, coordinate_basis, 2);
    
    assert!(two_comp_equal_up_to_theshold(left, 0.5* Comp::one()) );
    assert!(two_comp_equal_up_to_theshold(right, 0.5 * Comp::one()) );

    let (left,right) = get_comp_point_from_basis_and_vector(this_point, coordinate_basis, -2);
    
    assert!(two_comp_equal_up_to_theshold(left, 2.0* Comp::one()) );
    assert!(two_comp_equal_up_to_theshold(right, 2.0 * Comp::one()) );
}



pub fn testing_the_ellipse_paramters_randomly() 
{
    let mut rng = thread_rng();
    
    let x1 :Float = rng.gen_range(-1.0..1.0);
    let y1 :Float = rng.gen_range(-1.0..1.0);
    let epsilon: Float = rng.gen_range( 0.0000000000001..0.2 );

    let norm = (x1*x1 + y1*y1).sqrt();

    let (x,y) = (x1/norm, y1/norm);

    let rand_comp = Comp{re:x,im:y};

    println!(" ----------------------   ");
    println!("Random direction {}", rand_comp);
    println!("Random epsilon {}", epsilon);

    let (center,mat,radius) =    ellipse_parameters_for_region_a( rand_comp, epsilon);

    assert!(two_float_equal_up_to_threshold( center.re*rand_comp.re + center.im*rand_comp.im , 1.0-epsilon) );
    
    let a1: Float = rng.gen_range((1.0-epsilon)..1.0);
    let pythogoras_vert = (1.0 - a1*a1).sqrt();
    assert!( two_float_equal_up_to_threshold( a1*a1 + pythogoras_vert*pythogoras_vert, 1.0));
    
    let perp_comp = rand_comp*Comp{re:0.0,im:1.0};
    assert!( two_float_equal_up_to_threshold( perp_comp.re*rand_comp.re + perp_comp.im*rand_comp.im , 0.0));
    let a2: Float = rng.gen_range(-pythogoras_vert..pythogoras_vert);
    
    let test_point = rand_comp*a1 + perp_comp*a2;
    println!("Random test_point {}", test_point);

    // this random test_point should be in the ellipse
    let offset = test_point - center;
    let offset_vector  = nalgebra::matrix![ offset.re  , offset.im ];
    let offset_vector_trans = offset_vector.transpose();
    let offset_skewed = mat*offset_vector_trans;


    let skewed_distance = offset_skewed.transpose() * offset_skewed;

    println!("Testing {} < {}", skewed_distance[0] , radius);
    assert!( test_this_complex_pair_of_points(test_point, test_point, ( rand_comp, epsilon) )  );
    // assert!(skewed_distance[0] <= radius);
    println!(" ----------------------  \n\n\n ");
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


#[test]
pub fn inexact_synth_testing() 
{
    let answer = grid_problem(Comp::one(), 0.19);

    println!("{}", answer); }
