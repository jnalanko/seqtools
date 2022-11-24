use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();

    // Turn the string into a byte vector. The clone makes a new copy
    // so that we don't modify directly in args, which is immutable.
    // into_bytes transforms the clone into a byte array, consuming the clone,
    // so that we can use the same memory and there is no redundant copy.
    let mut bytes = args[1].clone().into_bytes();

    let n = bytes.len();

    // Complement
    for i in 0..n {
        bytes[i] = match bytes[i]{
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => {
                println!("Warning: non-ACGT character: '{}'", bytes[i] as char); 
                bytes[i]
            }
        };
    }

    // Reverse
    for i in 0..(n/2){
        let temp = bytes[i];
        bytes[i] = bytes[n-1-i];
        bytes[n-1-i] = temp;
    }

    println!("{}", std::str::from_utf8(bytes.as_slice()).expect(""));
}