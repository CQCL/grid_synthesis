use crate::algorithms::exact_synth_hashtable_lookup::generate_gate_table;
use crate::algorithms::exact_synth_hashtable_lookup::read_hash_table;

#[test]
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





