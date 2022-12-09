use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufRead, BufReader};

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

//fn main() {
//    let table = read_hash_table("table.txt").unwrap();
//    println!("{:?}", table);
//    // Output:
//    //     // {1: 2, 3: 4, 5: 6}
//    //     }
//    //


pub fn generate_gate_table() {
    let mut gatetable = HashMap::new();
    let file_to_be_saved_at = "data/gates_with_small_t_count.dat";


    gatetable.insert(1, 2);
    gatetable.insert(3, 4);
    gatetable.insert(5, 6);

    save_hash_table(&gatetable, file_to_be_saved_at ).unwrap();
}
