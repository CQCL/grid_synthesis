use crate::algorithms::lll::gram_schmidt_orthogonalization;
use crate::algorithms::lll::swap_columns;
use crate::algorithms::lll::nearest_plane;
use crate::algorithms::lll::lll_reduce;
use crate::algorithms::lll::size_reduce_given_gram_schmidt;
use crate::algorithms::lll::LLL_DELTA;
use crate::algorithms::inexact_synth::vec4int_to_vec4float;
use crate::algorithms::inexact_synth::mat4int_to_mat4;
use crate::algorithms::inexact_synth::mat4int_inverse;




use crate::structs::rings::Float;
use crate::structs::rings::Int;
type Mat4 = nalgebra::Matrix4<Float>;
type Mat4Int = nalgebra::Matrix4<Int>;
type Vec4 = nalgebra::Matrix4x1<Float>;
type Vec4Int = nalgebra::Matrix4x1<Int>;


use num_traits::One;

// Random number generators
use rand::thread_rng;
use rand::Rng;

// some auxilary functions from inexact synth tests are also useful here
use crate::tests::inexact_synth_tests::two_float_equal_up_to_threshold;
use crate::tests::inexact_synth_tests::generate_random_4by4_matrix;

// function to test floating point equality 
pub fn two_float_equal_up_to_given_threshold(left: Float, right: Float, epsilon :Float) -> bool
{
    return (left-right).abs() < epsilon;

}
#[test]
pub fn nearest_plane_test_randomly()
{

    let random_lattice = generate_random_4by4_matrix();
    let bstar = gram_schmidt_orthogonalization(random_lattice);

    let mut rng = thread_rng();
    let r = 1.0; //range
    let a1: Float = rng.gen_range(-r..r);
    let a2: Float = rng.gen_range(-r..r);
    let a3: Float = rng.gen_range(-r..r);
    let a4: Float = rng.gen_range(-r..r);

    let random_vector = Vec4::new(a1,a2,a3,a4);

    let (output,intvec) = nearest_plane(random_lattice, bstar,  random_vector);

    let diff = random_vector - output;


    // println!("input is {}",random_vector );
    // println!("diff is {}",diff );

    // take dot product
    // for some reason, rust won't let me use diff.dot( ?? )
    for j in 0..4
    {
        let mut dot = 0.0;
        for k in 0..4
        {
            dot = dot + diff[k] * bstar[(k,j)];
        }

        let norm = bstar.column(j).dot(&bstar.column(j));

        // println!("dot/norm is {} for j = {}", dot/norm, j );
        assert!( (dot/norm).abs() <= 0.5 );
    }

    let test_vec =  random_lattice * vec4int_to_vec4float(intvec);

    for k in 0..4
    {
        // println!("testing equality of {} {}",test_vec[k], output[k] );
        assert!( two_float_equal_up_to_threshold( test_vec[k], output[k] ) );
    }
}


#[test]
pub fn nearest_plane_rapid_testing()
{
    let n = 1000;
    for i in 1..n
    {
        nearest_plane_test_randomly();
    }

    println!("Nearest plane algorithm works {} times!", n );

}

#[test]
pub fn gram_schmidt_test()
{
    let n  = 1000;
    for test_this_a_lot_of_times in 0..n
    {

        let input = generate_random_4by4_matrix();
        let mat = gram_schmidt_orthogonalization(input);

        for i in 0..4
        {
            for j in 0..i
            {
                let dot = mat.column(i).dot(&mat.column(j));
                assert!(two_float_equal_up_to_threshold(dot,0.0) );
            }
        }
    }

    println!("Gram- Schmidt worked {} times! ",n );

}

pub fn is_size_reduced( b : Mat4) -> bool
{

    let bstar = gram_schmidt_orthogonalization(b);
    for i in 0..4 
    {
        for j in 0..i 
        {
            let vec = Vec4::new( b[(0,i)],  b[(1,i)],  b[(2,i)],  b[(3,i)] );

            let mut dot = 0.0;
            for k in 0..4
            {
                dot = dot + vec[k] * bstar[(k,j)];
            }

            let norm = bstar.column(j).dot(&bstar.column(j));


            if ( (dot/norm).abs() > 0.5 )
            {
                // println!("dot/norm is {}", dot/norm );
                // println!("dot is {}", dot);
                // println!("norm is {}", norm );
                // println!("(i,j) is ({},{})", i,j );
                println!("Size reduced failed " );
                return false;
            }
        }
    }

    return true
}



#[test]
pub fn random_size_reduced_test()
{


    let n = 1000;
    for i in 0..n
    {
        let mat = generate_random_4by4_matrix();
        let matstar = gram_schmidt_orthogonalization(mat);
        let (reduced,x) = size_reduce_given_gram_schmidt(mat,matstar);
        assert!(is_size_reduced(reduced));
        // println!(" x {}", x);
        // println!("mat {}", mat);
        let recover = reduced * mat4int_to_mat4(x);


        for j in 0..4
        {
            for k in 0..4
            {
                assert! ( two_float_equal_up_to_given_threshold( recover[(j,k)],mat[(j,k)], 0.000001) );
            }
        }
    }

    println!("size_reduce worked {} times!", n );

}

#[test]
pub fn test_lll_reduction_rapidly() 
{
    let n = 1000;
    for i in 0..n
    {
        let mat = generate_random_4by4_matrix();
        let (reduced,x) = lll_reduce(mat);


        let test_mat = reduced * mat4int_to_mat4(x) ;

        // The output matrices multiply to return the input matrix
        for j in 0..4
        {
            for k in 0..4
            {
                assert! ( two_float_equal_up_to_given_threshold( test_mat[(j,k)],mat[(j,k)], 0.000001) );
            }
        }

        assert!(is_size_reduced(reduced));
        assert!(test_lovasz_condition(reduced));


        let xinv = mat4int_inverse(x);
        assert_eq!( xinv * x , Mat4Int::identity() );
    }

    println!("lll_reduce worked {} times!", n );

}

pub fn test_lovasz_condition( b : Mat4) -> bool
{
    let bstar = gram_schmidt_orthogonalization(b);
    for i in 0..3
    {
        let left = bstar.column(i);
        let uij = b.column(i+1).dot(&bstar.column(i))/bstar.column(i).dot(&bstar.column(i)) ;
        let right = bstar.column(i) * uij + bstar.column(i+1);

        if LLL_DELTA * (left.dot(&left)) >=  (right.dot(&right))
        {
            return false;
        }
    }

    return true;
}



#[test]
pub fn swapping_works()
{

    // I only want to check swapping, it doens't matter what SQRT2 is here
    let SQRT2 = 1.414;

    let  mut mat: Mat4  = 
        Mat4::new( 1.0 ,   SQRT2  , 0.0  , 0.0,
            0.0 ,      0.0 , 1.0  ,  SQRT2 ,
            1.0 ,  -SQRT2  , 0.0  , 0.0,
            0.0 ,      0.0 , 1.0  , -SQRT2 );

    let  matswapped_first2: Mat4  = 
        Mat4::new( SQRT2,   1.0  , 0.0  , 0.0,
            0.0 ,      0.0 , 1.0  ,  SQRT2 ,
            -SQRT2,      1.0 , 0.0  , 0.0,
            0.0 ,      0.0 , 1.0  , -SQRT2 );

    swap_columns(0,1,&mut mat);

    assert_eq!( mat , matswapped_first2);
}
