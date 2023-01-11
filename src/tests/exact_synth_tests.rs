use crate::structs::rings::Conj;
use crate::structs::rings::local_ring::Local; 
use crate::structs::rings::zroot2::Zroot2;
// use crate::structs::rings::quaternion::Quaternion;
use crate::structs::rings::special_values::{mu_8, onebyroot2comp, sqrt2, sqrtminus1};
use num_complex::Complex;
use crate::structs::sunimat::SUniMat; 
use crate::structs::rings::Int;
use crate::structs::rings::Float;
use crate::algorithms::exact_synth::apply_gate_string_to_state;
use crate::algorithms::exact_synth::multiply_H_times_T_to_n;
use crate::algorithms::exact_synth::apply_t_gate;
use crate::algorithms::exact_synth::apply_tinv_gate;
use crate::algorithms::exact_synth::apply_h_gate;
use crate::algorithms::exact_synth::exact_synth;

use num_traits::{One, Zero};
use num_traits::Pow;

type Loc = Local<Zroot2>;
type Comp = Complex<Loc>;
type Mat = SUniMat<Comp>;
// type Quat = Quaternion<Loc>;
//


use crate::structs::unimat::ExactUniMat;

// Want to check if my gate_string gives output when applied to gamma
#[test]
pub fn apply_gate_string_to_states_and_check_output() 
{
    let g = Mat::one();

    // println!("Test 1:");
    let mut g1 = apply_gate_string_to_state("H".to_string(), g);
    let mut g2 = apply_gate_string_to_state("H".to_string(), g1);

    let mut g2prime = apply_gate_string_to_state("HH".to_string(), g);

    assert_eq!(g2,g2prime);
    assert_eq!(g,g2prime);
    assert_eq!(g,g2);

    // println!("Test 2:");
    g1 = apply_gate_string_to_state("H".to_string(), g);
    g2 = apply_gate_string_to_state("TTT".to_string(), g1);
    g2prime = apply_gate_string_to_state("TTTH".to_string(), g);
    // println!("{}", g2);
    // println!("{}", g2prime);
    assert_eq!(g2,g2prime);


    // println!("Test 3:");
    g1 = apply_gate_string_to_state("TTTTTTHTTTHTH".to_string(), g);
    g2 = apply_gate_string_to_state("TTH".to_string(), g1);

    g2prime = apply_gate_string_to_state("TTHTTTTTTHTTTHTH".to_string(), g);
    // println!("{}", g2);
    // println!("{}", g2prime);
    assert_eq!(g2,g2prime);

    // println!("Test 4:");
    g1 = apply_gate_string_to_state("TTTTT".to_string(), g);
    g2 = apply_gate_string_to_state("TTT".to_string(), g1);

    g2prime = apply_gate_string_to_state("TTTTTTTT".to_string(), g);
    // println!("{}", g2);
    // println!("{}", g2prime);
    assert_eq!(g2,g2prime);
    assert_eq!(g,g2);
    assert_eq!(g,g2prime);

}

#[test]
pub fn comp_pow_test() {
    let one = Comp::one();
    for exponent in -5..5 {
        assert_eq!(one.pow(exponent), one);
    }

    let root2 = sqrt2();
    assert_eq!(root2.pow( 3), root2 * root2 * root2);

    assert_eq!(root2.pow( -2), one / root2 / root2);
}

#[test]
pub fn multiply_H_times_T_to_n_test() {
    let mut gamma = Mat::one();
    let n = 5;

    let mat_2 = multiply_H_times_T_to_n(gamma, n);

    for i in 0..n {
        gamma = apply_t_gate(gamma);
    }

    gamma = apply_h_gate(gamma);
    assert_eq!(gamma.u, mat_2.u);
    assert_eq!(gamma.t, mat_2.t);
}

#[test]
pub fn apply_tinv_test() {
    let mut gamma = Mat::one();
    gamma = apply_t_gate(gamma);
    gamma = apply_tinv_gate(gamma);
    assert_eq!(gamma, Mat::one());
}



#[test]
fn exact_synth_given_norm_1_test_old() 
{
    let hadamard = Mat{
        u: sqrtminus1() / sqrt2(),
        t: sqrtminus1() / sqrt2(),
    };

    let t_gate = Mat{
        u: mu_8(),
        t: Comp::zero(),
    };

    // Note: As of Dec 12, 2022, I think matrix multiplication might be broken.
    let hththt = hadamard * t_gate * hadamard * t_gate * hadamard * t_gate;
    // let seq  = exact_synth_given_norm_1(hththt);
    println!(" ------ HADAMARD ------- \n  {} \n --------------------", hadamard);

    println!("----- ------");

    println!(" ------ T_Gate ------- \n  {} \n --------------------", t_gate);
    // assert_eq!(mat, hththt);
    // assert_eq!(seq, "HTHTHT"); // fails
    println!("{}", hththt);
    assert_eq!(Mat::one(), hththt);
    assert_eq!(Mat::one(), t_gate *t_gate *t_gate *t_gate *t_gate *t_gate * t_gate * t_gate);
    println!("{}", t_gate * t_gate * t_gate * t_gate);
    // assert_eq!(Mat::one(),t_gate *t_gate * t_gate * t_gate); // fails

    // let htht = hadamard * t_gate * hadamard * t_gate;
    // let seq = exact_synth_given_norm_1(htht);
    // let test_mat = apply_gate_string_to_state("HTHT".to_string(),Mat::one())* htht.inv() ;
    // // assert_eq!(seq, "HTHT"); // fails
}





pub fn exact_synth_tests_with_short_sequence( inputseq : String) {

    let inputgate  = ExactUniMat::from_string(&inputseq);

    // println!("This is the input \n : {}", inputgate);

    let outputseq = exact_synth(inputgate);

    let output = ExactUniMat::from_string(&outputseq);

    assert_eq!(outputseq, inputseq); // short sequences should be recovered optimally

    assert_eq!(output, inputgate);

}

pub fn exact_synth_tests_with_longer_sequence( inputseq : String) {

    let inputgate  = ExactUniMat::from_string(&inputseq);


    let outputseq = exact_synth(inputgate);

    let output = ExactUniMat::from_string(&outputseq);

    // assert_eq!(outputseq, inputseq); // this may or may not fail for long sequences
    println!("\n \n -------------------------- \n This is the input sequence: \n {}", inputseq);
    println!(" \n This is the gate I found!  : \n {}", outputseq);
    println!("-------------------------------- ");

    assert_eq!(output, inputgate);

    println!("And it checks out after multiplication");
    println!("-------------------------------- ");
}


#[test]
pub fn testing_exact_synth_rapidly_with_short_sequences() 
{

    exact_synth_tests_with_short_sequence( "H".to_string() ) ;
    exact_synth_tests_with_short_sequence( "T".to_string() ) ;
    exact_synth_tests_with_short_sequence( "HT".to_string() ) ;
    exact_synth_tests_with_short_sequence( "HTHT".to_string() ) ;
    exact_synth_tests_with_short_sequence( "TTTHTT".to_string() ) ;
    exact_synth_tests_with_short_sequence( "TTTHTT".to_string() ) ;
    exact_synth_tests_with_short_sequence( "TTTHTT".to_string() ) ;
    exact_synth_tests_with_short_sequence( "TTTHTT".to_string() ) ;

}


#[test]
pub fn testing_exact_synth_rapidly_with_long_sequences() 
{


    // Lookup is failing for these
    exact_synth_tests_with_short_sequence( "HTTTHTTTHTT".to_string() ) ;
    exact_synth_tests_with_short_sequence( "HTTTHTTTHTTTH".to_string() ) ;
    exact_synth_tests_with_longer_sequence( "THTTTHTHTTTTTTTHTHT".to_string() ) ;
    exact_synth_tests_with_longer_sequence( "TTTHTTHTTTHTTTHTTTH".to_string() ) ;
    exact_synth_tests_with_longer_sequence( "TTTHTTHTTTHTTHTTTHTHTTTH".to_string() ) ;
    exact_synth_tests_with_longer_sequence( "TTTHTHTTTHTTTHTTTHTHTTTHTTHTTTHTHTTTH".to_string() ) ;
    exact_synth_tests_with_longer_sequence( "TTTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTH".to_string() ) ;
    exact_synth_tests_with_longer_sequence( "THTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTHTTTTTTTHTHT".to_string() ) ;

    exact_synth_tests_with_longer_sequence( "THTTHTTTHTTHTTTHTTHTTTHTHTTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTH".to_string() ) ;

    exact_synth_tests_with_longer_sequence( "TTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTH".to_string() ) ;

    exact_synth_tests_with_longer_sequence( "THTTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTHTTTTTTTHTHTHTTTHTTHTTTHTTHTTTHTHTTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTH".to_string() ) ;

    exact_synth_tests_with_longer_sequence( "TTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTHTTTTTTTHTHTTTHTTHTTHTTTHTTTHTTH".to_string() ) ;

    exact_synth_tests_with_longer_sequence( "THTTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTHTTTTTTTHTHTHTTTHTTHTTTHTTHTTTHTHTTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTH".to_string() ) ;

    exact_synth_tests_with_longer_sequence( "THTTTHTTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTHTTTTTTTHTHTHTTTHTTHTTTHTTHTTTHTHTTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTHTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTHTTTTTTTHTHTHTTTHTTHTTTHTTHTTTHTHTTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTH".to_string() ) ;

    exact_synth_tests_with_longer_sequence( "THTTHTTHTTTHTTTHTTTHTTHTTTHTTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTHTTTTTTTHTHTHTTTHTTHTTTHTTHTTTHTHTTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTHHTTTHTTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTTHTTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTHTTTTTTTHTHTHTTTHTTHTTTHTTHTTTHTHTTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTHTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTHTTTTTTTHTHTHTTTHTTHTTTHTTHTTTHTHTTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTH".to_string() ) ;

    exact_synth_tests_with_longer_sequence( "THTTHTTHTTTHTTTHTTTHTTHTTTHTTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTTHTHTTTHTTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTHTTTTTTTHTHTHTTTHTTHTTTHTTHTTTHTHTTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTHTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTHTTTTTTTHTHTHTTTHTTHTTTHTTHTTTHTHTTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTHTTTHTTTHTTHTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTHTTTTTTTHTHTHTTTHTTHTTTHTTHTTTHTHTTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTHHTTTHTTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTTHTTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTHTTTTTTTHTHTHTTTHTTHTTTHTTHTTTHTHTTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTHTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTHTTTTHTHTHTTTTTHTTHTTHTHTTTTTTTHTHTHTTTHTTHTTTHTTHTTTHTHTTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTH".to_string() );

    exact_synth_tests_with_longer_sequence( "THTTHTTHTTTHTTTHTTTHTTHTTTHTTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTTTHHTTTHTTHTHTHTTTTTHTTTHTTTTTHHHHTHTHTHTTTHTTHTTTHTTTHTTTHTHTTTHTTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTHTTTTTTTHTHTHTTTHTTHTTTHTTHTTTHTHTTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTHTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTHTTTTTTTHTHTHTTTHTTHTTTHTTHTTTHTHTTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTHTTTHTTTHTTHTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTHTTTTTTTHTHTHTTTHTTHTTTHTTHTTTHTHTTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTHHTTTHTTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTTHTTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTHTTTTTTTHTHTHTTTHTTHTTTHTTHTTTHTHTTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTHTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTTTHTTHTTHTTHTTHTTTHTTTHTTTHTTHTTHTTTHTHTTTTHTHTHTTTTTHTTHTTHTHTTTTTTTHTHTHTTTHTTHTTTHTTHTTTHTHTTTHTTHTTTHTTTHTTHTTHTTTHTTTHTTH".to_string() );

}
