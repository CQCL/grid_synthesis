use crate::algorithms::exact_synth_hashtable_lookup::has_repeated_zeroes;
use crate::algorithms::exact_synth_hashtable_lookup::generate_gate_table;
use crate::algorithms::exact_synth_hashtable_lookup::read_hash_table;

// #[test]
pub fn print_hashtable()
{
    generate_gate_table();

    let filename = "data/gates_with_small_t_count.dat";
    let hashtable_from_file = read_hash_table(filename).unwrap();


    for (key,value) in hashtable_from_file 
    {
        println!("Table Entry: {} \t {}", key,value);
    }


}










#[test]
fn has_repeated_zeroes_test() {
	// NOTE: We cut off leading zeros.
	// Ensure this behaviour is intended!
	assert_eq!(has_repeated_zeroes(0b1010101010101010), false);
	assert_eq!(has_repeated_zeroes(0b1010100110101010), true);
	assert_eq!(has_repeated_zeroes(0b0000000000000001), false);
	assert_eq!(has_repeated_zeroes(0b1111111111111111), false);
}
