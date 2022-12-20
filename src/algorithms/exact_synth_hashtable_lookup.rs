use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufRead, BufReader};

use crate::structs::rings::zroot2::Zroot2;
type Loc = crate::structs::rings::local_ring::Local::<Zroot2>;
type ExactGate = crate::structs::sunimat::UniMat<Loc>;

type GateString = u32; // Binary string of 0 and 1
// We will represent 0 as H gate and 1 as T gate

type GateTable = HashMap< GateString, ExactGate >;

const GATE_STRING_LENGTH : usize = 16; 


// This reading and writing of the gate table file could be done using a serialization.
// There is a serde package in Rust to do this. 
// This would require writing a Serialize trait for our ExactGate type
// We will leave this for future work
pub fn read_hash_table(filename: &str) -> std::io::Result<HashMap<i32, i32>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let mut table = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split(" ").collect();

        let key: i32 = parts[0].parse().unwrap();

        let index = 1; // If we want to move ahead later

        // Reconstruct from the data 
        // Since we assume that the data is written by save_hash_table
        // we do not need the constructors Loc::from_base 
        let value  = ExactGate
        {
            u: Loc
            {
                num: Zroot2
                {
                    0: parts[index+0].parse().unwrap(),
                    1: parts[index+1].parse().unwrap()
                },
                log_den: parts[index+2].parse().unwrap()
            },
            t: Loc
            {
                num: Zroot2
                {
                    0: parts[index+3].parse().unwrap(),
                    1: parts[index+4].parse().unwrap()
                },
                log_den: parts[index+5].parse().unwrap()
            }
        };
            
        // table.insert(key, value);
    }

    Ok(table)
}


pub fn save_hash_table(table: &GateTable, filename: &str) -> std::io::Result<()> {
    let mut file = File::create(filename)?;

    for (key, value) in table 
    {
        // Changing how the data is written should be done hand in hand with how the data is read 
        // See read_hash_table
        writeln!(file, "{} {} {} {} {} {} {}", key, value.u.num.0, value.u.num.1, value.u.log_den,value.t.num.0,value.t.num.1,value.t.log_den )?;
    }

    Ok(())
}

pub fn push_zero( sequence: &mut GateString, bits_written: &mut usize)
{
    *sequence = *sequence << 1;
    *bits_written +=1;
}

pub fn push_one( sequence: &mut GateString, bits_written: &mut usize)
{
    *sequence = *sequence << 1 | 1;
    *bits_written +=1;
}

pub fn pop_sequence( sequence: &mut GateString, bits_written: &mut usize)
{
    *sequence = *sequence >> 1;
    *bits_written -=1;
}

// WARNING: final_length given should be less than 8*size_of(GateString)
pub fn generate_hashtable_of_gate_sequence_recursively(final_length: &usize, bits_written: &mut usize ,sequence: &mut GateString, table: &mut GateTable )
{

    if *bits_written == *final_length
    {
        // WRITE SEQUENCE TO HASHTABLE
        let temp = crate::structs::sunimat::UniMat::one();

        table.insert(*sequence,temp);


        return;
    }

    push_zero(sequence, bits_written);

    if *bits_written < 2 || *sequence & 2  != 0 
    {
        generate_hashtable_of_gate_sequence_recursively(final_length, bits_written, sequence, table);
    }
    
    pop_sequence(sequence, bits_written);
    push_one(sequence, bits_written);
    
    generate_hashtable_of_gate_sequence_recursively(final_length, bits_written, sequence, table);
    
    pop_sequence(sequence, bits_written);
}


pub fn generate_gate_table() {
    let mut gatetable = GateTable::new();
    let file_to_be_saved_at = "data/gates_with_small_t_count.dat";
    generate_hashtable_of_gate_sequence_recursively( &GATE_STRING_LENGTH , &mut 0, &mut 0, &mut gatetable);

    save_hash_table(&gatetable, file_to_be_saved_at ).unwrap();
}
