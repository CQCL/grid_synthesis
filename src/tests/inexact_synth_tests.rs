use crate::algorithms::inexact_synth::grid_problem;
use crate::algorithms::inexact_synth::test_this_integer_point;
use crate::algorithms::inexact_synth::get_comp_point_from_basis_and_vector;


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

#[test]
pub fn basic_vicinity_test()
{

    let this_point = Vec4Int::new(1,0,0,0);
    let coordinate_basis = Mat4::new(1.0,0.0,0.0,0.0,
                                     0.0,1.0,0.0,0.0,
                                     0.0,0.0,1.0,0.0,
                                     0.0,0.0,0.0,1.0);

    let result = test_this_integer_point( this_point , coordinate_basis , 0 , ( Comp::one(), 0.1) ) ;

    println!("{}",this_point);
    assert!(result);

}

#[test]
pub fn comp_from_basis_and_vector_test()
{
    let this_point = Vec4Int::new(1,0,0,0);
    let coordinate_basis = Mat4::new(1.0,0.0,0.0,0.0,
                                     0.0,1.0,0.0,0.0,
                                     0.0,0.0,1.0,0.0,
                                     0.0,0.0,0.0,1.0);

    
    let (left,right) = get_comp_point_from_basis_and_vector(this_point, coordinate_basis, 2);
    println!("{}", left);
    println!("{}", right);
             

}


#[test]
pub fn inexact_synth_testing() 
{
    let answer = grid_problem(Comp::one(), 0.19);

    println!("{}", answer);
}
