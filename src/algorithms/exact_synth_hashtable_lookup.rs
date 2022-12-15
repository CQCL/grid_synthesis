use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufRead, BufReader};


type UnsignedBits = u16;

type SizeOfBits = u32;

pub fn read_hash_table(filename: &str) -> std::io::Result<HashMap<i32, i32>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let mut table = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split(" ").collect();
        let key: i32 = parts[0].parse().unwrap();
        let value: i32 = parts[1].parse().unwrap();
        table.insert(key, value);
    }

    Ok(table)
}


pub fn save_hash_table(table: &HashMap<i32, i32>, filename: &str) -> std::io::Result<()> {
    let mut file = File::create(filename)?;

    for (key, value) in table 
    {
        writeln!(file, "{} {}", key, value)?;
    }

    Ok(())
}



pub fn has_repeated_zeroes(n: UnsignedBits, size: SizeOfBits ) -> bool {
    // Convert the integer to a binary string

    for i in 0..(size-1)
    {
        let mask1 = 1 << i;
        let mask2 = 1 << i+1;
        if mask1 & n == 0 && mask2 & n == 0
        {
            return true;
        }
    }

    return false;
}



pub fn generate_gate_table() {
    let mut gatetable = HashMap::new();
    let file_to_be_saved_at = "data/gates_with_small_t_count.dat";

    let mut i : UnsignedBits = 0;
    while i < 65535
    {   

        gatetable.insert(i as i32 , i as i32);
        i = i+1;

    }
    save_hash_table(&gatetable, file_to_be_saved_at ).unwrap();
}
