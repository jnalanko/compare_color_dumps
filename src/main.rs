use std::{collections::HashMap, fs::File, io::{BufRead, BufReader}, path::Path};
use jseqio::{record::Record, seq_db::SeqDB};

fn ascii_to_int(ascii: &[u8]) -> usize {
    std::str::from_utf8(ascii)
    .expect("Unitig id is not valid utf-8").parse()
    .expect("Could not convert unitig id string to unsigned integer")
}

fn get_color_set_id(fasta_header: &[u8]) -> usize {
    // The fasta header should look like this " unitig_id=0 color_set_id=0".
    // Note the space at the start.

    let first_part = fasta_header[1..].split(|c| *c == b' ').next().unwrap();
    let mut tokens = first_part.split(|c| *c == b'=');
    assert_eq!(tokens.next().expect("Unitig id missing"), b"unitig_id");
    ascii_to_int(tokens.next().expect("Unitig id missing"))
}

fn read_color_sets(filename: impl AsRef<Path>, num_color_sets: usize) -> Vec<Vec<usize>> {
    // Lines should look like this:
    // color_set_id=9 size=7 3 4 9 12 14 15 16

    let mut color_sets = vec![vec![]; num_color_sets];

    let mut reader = BufReader::new(File::open(filename).unwrap());
    let mut line = String::new();

    while reader.read_line(&mut line).unwrap() > 0 {
        let line_bytes = line.trim_end().as_bytes();
        let mut tokens = line_bytes.split(|c| *c == b' ');

        let first_token = tokens.next().unwrap();
        assert_eq!(&first_token[0..13], b"color_set_id=");
        let color_set_id: usize = ascii_to_int(&first_token[13..]);
        assert!(color_set_id < num_color_sets);

        let second_token = tokens.next().unwrap();
        assert_eq!(&second_token[0..5], b"size=");
        let list_len = ascii_to_int(&second_token[5..]);

        color_sets[color_set_id] = tokens.map(ascii_to_int).collect();
        assert_eq!(color_sets[color_set_id].len(), list_len);

        line.clear();
    }

    color_sets
}

// Returns unitig id if the unitig is cyclic and a cyclic rotation is found in the hashmap
#[allow(dead_code)]
fn check_cyclic_rotations<'a>(unitig: &[u8], k: usize, unitig_to_id: &HashMap<&[u8], usize>) -> Option<usize> {
    if &unitig[0..k-1] == &unitig[unitig.len()-k+1..unitig.len()] {
        // Cyclic unitig -> check rotations
        let mut extended = unitig.to_vec();
        for i in 0..unitig.len() {
            extended.push(extended[k-1+i]);
            if let Some(id) = unitig_to_id.get(&extended[1+i..]) {
                return Some(*id);
            }
        }

    }
    None
}

fn main() {
    let mut args = std::env::args();
    args.next().unwrap(); // Program name
    let dump_A_file_prefix = args.next().unwrap();
    let dump_B_file_prefix = args.next().unwrap();


}