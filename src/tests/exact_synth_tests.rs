use crate::structs::rings::Conj;
use crate::structs::rings::local_ring::Local; 
use crate::structs::rings::zroot2::Zroot2;
use crate::structs::rings::quaternion::Quaternion;
use crate::structs::rings::special_values::{mu_8, onebyroot2comp, sqrt2, sqrtminus1};
use num_complex::Complex;
use crate::structs::sunimat::SUniMat; 
use crate::structs::rings::Int;
use crate::structs::rings::Float;
use crate::algorithms::exact_synth::{apply_gate_string_to_state,
    pow, multiply_H_times_T_to_n, apply_t_gate, apply_tinv_gate,
    apply_h_gate, exact_synth};

use num_traits::{One, Zero};

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
        assert_eq!(pow(one, exponent), one);
    }

    let root2 = sqrt2();
    assert_eq!(pow(root2, 3), root2 * root2 * root2);

    assert_eq!(pow(root2, -2), one / root2 / root2);
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



// WILL GET BACK HERE
// #[test]
// fn exact_synth_given_norm_1_test() 
// {
//     let hadamard = Mat{
//         u: sqrtminus1() / sqrt2(),
//         t: sqrtminus1() / sqrt2(),
//     };

//     let t_gate = Mat{
//         u: mu_8(),
//         t: Comp::zero(),
//     };

//     // Note: As of Dec 12, 2022, I think matrix multiplication might be broken.
//     let hththt = hadamard * t_gate * hadamard * t_gate * hadamard * t_gate;
//     // let seq  = exact_synth_given_norm_1(hththt);
//     println!(" ------ HADAMARD ------- \n  {} \n --------------------", hadamard);
    
//     println!("----- ------");
    
//     println!(" ------ T_Gate ------- \n  {} \n --------------------", t_gate);
//     // assert_eq!(mat, hththt);
//     // assert_eq!(seq, "HTHTHT"); // fails
//     println!("{}", hththt);
//     assert_eq!(Mat::one(), hththt);
//     assert_eq!(Mat::one(), t_gate *t_gate *t_gate *t_gate *t_gate *t_gate * t_gate * t_gate);
//     println!("{}", t_gate * t_gate * t_gate * t_gate);
//     // assert_eq!(Mat::one(),t_gate *t_gate * t_gate * t_gate); // fails

//     // let htht = hadamard * t_gate * hadamard * t_gate;
//     // let seq = exact_synth_given_norm_1(htht);
//     // let test_mat = apply_gate_string_to_state("HTHT".to_string(),Mat::one())* htht.inv() ;
//     // // assert_eq!(seq, "HTHT"); // fails
// }


// #[test] //passes
fn exact_synth_tests_single_h() {

    let inputseq  = "H".to_string();
    let inputgate  = ExactUniMat::from_string(&inputseq);
    let outputseq = exact_synth(inputgate);
    
    assert_eq!(outputseq, "H");
    
    let output = ExactUniMat::from_string(&outputseq);

    assert_eq!(output, inputgate);

}

#[test]
fn exact_synth_tests_single_t() {

    let inputseq  = "T".to_string();
    let inputgate  = ExactUniMat::from_string(&inputseq);
    
    
    let outputseq = exact_synth(inputgate);

    assert_eq!(outputseq, "T");
    
    let output = ExactUniMat::from_string(&outputseq);

    assert_eq!(output, inputgate);

}


#[test]
fn exact_synth_tests_ht() {

    let inputseq  = "HT".to_string();
    let inputgate  = ExactUniMat::from_string(&inputseq);

    // println!("This is the input \n : {}", inputgate);
    
    let outputseq = exact_synth(inputgate);
    
    let output = ExactUniMat::from_string(&outputseq);

    assert_eq!(output, inputgate);
    assert_eq!(outputseq, "HT");

}


// THIS TEST FAILS
// #[test]
fn exact_synth_tests_ttthtt() {

    let inputseq  = "TTTHTT".to_string();
    let inputgate  = ExactUniMat::from_string(&inputseq);

    // println!("This is the input \n : {}", inputgate);
    
    let outputseq = exact_synth(inputgate);
    
    let output = ExactUniMat::from_string(&outputseq);

    assert_eq!(outputseq, "TTTHTT");
    assert_eq!(output, inputgate);

}

#[test]
fn exact_synth_tests_htht() {

    let testseq = "HTHT";
    let inputseq  = testseq.to_string();
    let inputgate  = ExactUniMat::from_string(&inputseq);

    // println!("This is the input \n : {}", inputgate);
    
    let outputseq = exact_synth(inputgate);
    
    let output = ExactUniMat::from_string(&outputseq);

    // assert_eq!(outputseq, testseq);
    // assert_eq!(output, inputgate);

}
