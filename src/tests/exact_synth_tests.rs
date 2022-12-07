
use crate::structs::rings::Conj;
use crate::structs::rings::local_ring::Local; 
use crate::structs::rings::zroot2::Zroot2;
use crate::structs::rings::quaternion::Quaternion;
use crate::structs::rings::special_values::mu_8;
use crate::structs::rings::special_values::onebyroot2comp;
use crate::structs::rings::special_values::sqrt2;

use num_complex::Complex;
use crate::structs::sunimat::UniMat; 
use crate::structs::rings::Int;
use crate::structs::rings::Float;
use crate::algorithms::exact_synth::apply_gate_string_to_state;

type Loc = Local<Zroot2>;
type Comp = Complex<Loc>;
type Mat = UniMat<Comp>;
// type Quat = Quaternion<Loc>;

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



