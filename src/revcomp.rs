use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    let bytes = args[1].as_bytes();

    let revcomp = bytes.iter().map(|x| match *x {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => {
            println!("Warning: non-ACGT character: '{}'", *x as char); 
            *x
        }
    }).rev().collect::<Vec<u8>>();

    println!("{}", std::str::from_utf8(revcomp.as_slice()).expect(""));
}
