use crate::structs::unimat::ExactUniMat;
use crate::structs::rings::local_ring::Local; 
use crate::structs::rings::zroot2::Zroot2;
use crate::structs::sunimat::SUniMat; 
use crate::algorithms::exact_synth::exact_synth;
use crate::algorithms::exact_synth::apply_gate_string_to_state;
use num_complex::Complex;
type Loc = Local<Zroot2>;
type Comp = Complex<Loc>;
type Mat = SUniMat<Comp>;


#[test]
pub fn basic_identies() 
{
    let one = ExactUniMat::one();

    assert_eq!(one.inv(), one);
    assert_eq!(one*one, one);
    assert_eq!(one*one*one, one);

    let h_gate = ExactUniMat::h_gate();
    let t_gate = ExactUniMat::t_gate();

    assert_eq!(h_gate*h_gate.inv() , one);
    assert_eq!(h_gate*h_gate , one);

    assert_eq!(t_gate*t_gate.inv() , one);
    assert_eq!(t_gate.inv() *t_gate , one);

    assert_eq!(t_gate * t_gate * t_gate * t_gate * t_gate * t_gate * t_gate * t_gate , one);


    
    let temp = t_gate * h_gate * t_gate * h_gate * t_gate ;
    assert_eq!( (temp)  *  (t_gate.inv() * h_gate * t_gate.inv() * h_gate * t_gate.inv() ), one);
}

#[test]
pub fn more_identities()
{
    let h_gate = ExactUniMat::h_gate();
    let t_gate = ExactUniMat::t_gate();

    let temp = t_gate * t_gate;

    assert_eq!( (temp).inv()  *  (temp), ExactUniMat::one());

    
    
    let temp = t_gate * h_gate * t_gate * h_gate;
    assert_eq!( (temp).inv()  *  (temp), ExactUniMat::one());

}

pub fn exact_test_with_string(test_seq: String) 
{

    let prod_gate = ExactUniMat::from_string(&test_seq );
    let state_expected = apply_gate_string_to_state(test_seq , Mat::one());

    println!("{}", prod_gate.mat);
    println!("{}", prod_gate);
    println!("{}", state_expected);

    // let guess = 5;
    // assert_eq!(prod_gate.omega_exp, guess);

    assert_eq!(prod_gate.mat, state_expected);
}

#[test]
pub fn exact_test_with_some_strings()
{
    exact_test_with_string("TH".to_string() );
    exact_test_with_string("TT".to_string() );
    exact_test_with_string("HH".to_string() );
    exact_test_with_string("HT".to_string() );
    exact_test_with_string("TTHHTTTTTH".to_string() );
    exact_test_with_string("TTHHTTHTTTHTHTTTH".to_string() );
}



#[test]
pub fn multiply_h_and_t() 
{
    let h_gate = ExactUniMat::h_gate();
    let t_gate = ExactUniMat::t_gate();
    let prod = h_gate * t_gate;


    
    println!("h_gate \n = {}", h_gate);
    println!("t_gate \n = {}", t_gate);
    println!("h_gate * t_gate \n = {}", prod );

    assert_eq!(prod.omega_exp , 5) ;

    let ht_state=  apply_gate_string_to_state("HT".to_string(),Mat::one());

    assert_eq!(ht_state, prod.mat);

    assert_eq!(prod, ExactUniMat::from_string(&"HT".to_string()) );

}
