use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufRead, BufReader};


// Matrix multiplication
use crate::algorithms::exact_synth::apply_gate_string_to_state;

// Special values
use crate::structs::rings::special_values::mu_8;
use crate::structs::rings::special_values::sqrt2;


// Rings and matrices
use crate::structs::rings::zroot2::Zroot2;
type Loc = crate::structs::rings::local_ring::Local::<Zroot2>;
use num_complex::Complex;
type ExactGate = crate::structs::sunimat::SUniMat<Complex<Loc>>;

// A state is a complex linear combination of |0> and |1>
// such that the sum of norm_sq is 1
// This is exactly the same information as a Special Unitary Matrix
type State = ExactGate;


type GateString = String; // String of H and T
pub const GATE_STRING_LENGTH : usize = 15;  // This is the number of H and T gates to be multiplied


// Hash table. It will store gatestring and corresponding exact gate
pub type GateTable = HashMap< ExactGate,GateString,>;



// This reading and writing of the gate table file could be done using a serialization.
// There is a serde package in Rust to do this. 
// This would require writing a Serialize trait for our ExactGate type
// We will leave this for future work
pub fn read_hash_table(filename: &str) -> std::io::Result<GateTable> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let mut table = HashMap::new();


    for line in reader.lines() 
    {
        let line = line?;
        let parts: Vec<&str> = line.split(" ").collect();

        let value: GateString = parts[0].parse().unwrap();

        let index = 1; // If we want to move ahead later

        // Reconstruct from the data 
        // Since we assume that the data is written by save_hash_table
        // we do not need the constructors Loc::from_base 
        let key  = ExactGate
        {
            u: Complex
            {
                re: Loc
                {
                    num: Zroot2
                    {
                        0: parts[index+0].parse().unwrap(),
                        1: parts[index+1].parse().unwrap()
                    },
                    log_den: parts[index+2].parse().unwrap()
                },
                im: Loc
                {
                    num: Zroot2
                    {
                        0: parts[index+3].parse().unwrap(),
                        1: parts[index+4].parse().unwrap()
                    },
                    log_den: parts[index+5].parse().unwrap()
                }

            },
            
            t: Complex
            {
                re: Loc
                {
                    num: Zroot2
                    {
                        0: parts[index+6].parse().unwrap(),
                        1: parts[index+7].parse().unwrap()
                    },
                    log_den: parts[index+8].parse().unwrap()
                },
                im: Loc
                {
                    num: Zroot2
                    {
                        0: parts[index+9].parse().unwrap(),
                        1: parts[index+10].parse().unwrap()
                    },
                    log_den: parts[index+11].parse().unwrap()
                }

            },
        };

        table.insert(key, value);


    }


    Ok(table)
}


pub fn save_hash_table(table: &GateTable, filename: &str) -> std::io::Result<()> {
    let mut file = File::create(filename)?;

    for (key, value) in table 
    {
        // Changing how the data is written should be done hand in hand with how the data is read 
        // See read_hash_table
        writeln!(file, "{} {} {} {} {} {} {} {} {} {} {} {} {}", 
                 value, 
                 key.u.re.num.0,
                 key.u.re.num.1,
                 key.u.re.log_den,
                 key.u.im.num.0,
                 key.u.im.num.1,
                 key.u.im.log_den,
                 key.t.re.num.0,
                 key.t.re.num.1,
                 key.t.re.log_den,
                 key.t.im.num.0,
                 key.t.im.num.1,
                 key.t.im.log_den 
                 )?;
    }

    Ok(())
}

pub fn push_h_gate( sequence: &mut GateString)
{
    sequence.push('H');
}

pub fn push_t_gate( sequence: &mut GateString)
{
    sequence.push('T');
}

pub fn pop_sequence( sequence: &mut GateString)
{
    sequence.pop();
}

// WARNING: final_length given should be less than 8*size_of(GateString)
// otherwise there shall be overflow
pub fn generate_hashtable_of_gate_sequence_recursively(final_length: &usize, sequence: &mut GateString, table: &mut GateTable )
{

    let mut len = sequence.len();
    if len == *final_length
    {
        // write sequence to hashtable
        let ident = ExactGate::one();
        let temp = apply_gate_string_to_state(sequence.to_string(), ident);

        table.insert(temp, sequence.to_string());

        return;
    }

    push_h_gate(sequence);
    len = len+1;
    
    // avoiding HH in gate sequence
    if len > 1
    {
        if sequence[(len-2)..(len)] != "HH".to_string()
        {
            generate_hashtable_of_gate_sequence_recursively(final_length,  sequence, table);
        }
    }
    {
        generate_hashtable_of_gate_sequence_recursively(final_length,  sequence, table);
    }

    pop_sequence(sequence);
    
    push_t_gate(sequence);

    // Avoiding TTTTTTTT in gate sequence
    if len > 7
    {
        if sequence[(len-8)..(len)] != "TTTTTTTT".to_string()
        {
            generate_hashtable_of_gate_sequence_recursively(final_length,  sequence, table);
        }
    }
    else
    {
        generate_hashtable_of_gate_sequence_recursively(final_length,  sequence, table);
    }

    pop_sequence(sequence);
}



pub fn generate_gate_table() {
    let mut gatetable = GateTable::new();
    let file_to_be_saved_at = "data/gates_with_small_t_count.dat";


    for i in 0..GATE_STRING_LENGTH
    {
        let gatelength = GATE_STRING_LENGTH -i ; // This reverse length makes shorter gate sequences overwrite longer ones

        generate_hashtable_of_gate_sequence_recursively( &gatelength  ,  &mut "".to_string(), &mut gatetable);
    }
    
    // Updating the value of the identity state
    let option = gatetable.insert(ExactGate::one(), "I".to_string());

    
    save_hash_table(&gatetable, file_to_be_saved_at ).unwrap();
        
    // Adding the empty gate sequence

}
