use crate::algorithms::exact_synth_hashtable_lookup::generate_gate_table;
use crate::algorithms::exact_synth_hashtable_lookup::read_hash_table;
use crate::algorithms::exact_synth_hashtable_lookup::GATE_STRING_LENGTH;

// Rings and matrices
use crate::structs::rings::zroot2::Zroot2;
type Loc = crate::structs::rings::local_ring::Local::<Zroot2>;
use num_complex::Complex;
type ExactState = crate::structs::sunimat::SUniMat<Complex<Loc>>;

// Matrix multiplication
use crate::algorithms::exact_synth::apply_gate_string_to_state;



// This test will also generate the gate table
// #[test]
pub fn print_hashtable()
{
    generate_gate_table();

    let filename = "data/gates_with_small_t_count.dat";
    let hashtable_from_file = read_hash_table(filename).unwrap();

    let mut count = 0;

    for (key,value) in hashtable_from_file 
    {
        assert_eq!(apply_gate_string_to_state(value,ExactState::one()), key);
        count+=1;
    }

}


#[test]
pub fn checking_hashtable_values_easy_cases() 
{
    // Assuming that hashtable exists
    // we will check that some values are what we expect

    let filename = "data/gates_with_small_t_count.dat";
    let hashtable = read_hash_table(filename).unwrap();

    let state = apply_gate_string_to_state("HT".to_string(), ExactState::one());
    let seq = hashtable.get(&state).unwrap();
    assert_eq!(seq, "H");
    
    let state = apply_gate_string_to_state("HTTT".to_string(), ExactState::one());
    let seq = hashtable.get(&state).unwrap();
    assert_eq!(seq, "H");
    
    let state = apply_gate_string_to_state("T".to_string(), ExactState::one());
    assert_eq!(state, ExactState::one());
    let seq = hashtable.get(&state).unwrap();
    assert_eq!(seq, "I");
    
    let state = apply_gate_string_to_state("TT".to_string(), ExactState::one());
    let seq = hashtable.get(&state).unwrap();
    assert_eq!(seq, "I");
    
    let state = apply_gate_string_to_state("TTH".to_string(), ExactState::one());
    let seq = hashtable.get(&state).unwrap();
    assert_eq!(seq, "TTH");
    
    let state = apply_gate_string_to_state("TH".to_string(), ExactState::one());
    let seq = hashtable.get(&state).unwrap();
    assert_eq!(seq, "TH");
    

}

#[test]
pub fn checking_hashtable_values_more_difficult_cases()
{
    let filename = "data/gates_with_small_t_count.dat";
    let hashtable = read_hash_table(filename).unwrap();



    let input_string = "TTHT";
    let expected_answer = "TTH";
    let state = apply_gate_string_to_state(input_string.to_string(), ExactState::one());
    let seq = hashtable.get(&state).unwrap();
    assert_eq!(seq, expected_answer);


    let input_string = "TTTHTT";
    let expected_answer = "TTTH";
    let state = apply_gate_string_to_state(input_string.to_string(), ExactState::one());
    let seq = hashtable.get(&state).unwrap();
    assert_eq!(seq, expected_answer);

}




