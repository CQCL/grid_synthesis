use crate::algorithms::exact_synth_hashtable_lookup::generate_gate_table;
use crate::algorithms::exact_synth_hashtable_lookup::read_hash_table;
use crate::algorithms::exact_synth_hashtable_lookup::GATE_STRING_LENGTH;

// Rings and matrices
use crate::structs::rings::zroot2::Zroot2;
type Loc = crate::structs::rings::local_ring::Local::<Zroot2>;
use num_complex::Complex;
type ExactGate = crate::structs::sunimat::SUniMat<Complex<Loc>>;

// Matrix multiplication
use crate::algorithms::exact_synth::apply_gate_string_to_state;



// This test will also generate the gate table
#[test]
pub fn print_hashtable()
{
    generate_gate_table();

    let filename = "data/gates_with_small_t_count.dat";
    let hashtable_from_file = read_hash_table(filename).unwrap();

    let mut count = 0;

    for (key,value) in hashtable_from_file 
    {
        assert_eq!(apply_gate_string_to_state(value,ExactGate::one()), key);
        count+=1;
    }


}





