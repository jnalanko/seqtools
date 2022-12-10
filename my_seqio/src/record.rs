
use std::fmt;
use std::str;

pub trait Record{
    fn head(&self) -> &[u8];
    fn seq(&self) -> &[u8];
    fn qual(&self) -> Option<&[u8]>;
}

#[derive(Debug)]
pub struct SeqRecord<'a>{
    pub head: &'a [u8],    
    pub seq: &'a [u8],
    pub qual: Option<&'a [u8]>, // If FASTA, this is None
}

#[derive(Debug)]
pub struct OwnedSeqRecord{
    pub head: Vec<u8>,    
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>, // If FASTA, this is None
}

impl<'a> Record for SeqRecord<'a>{
    fn head(&self) -> &[u8]{self.head}
    fn seq(&self) -> &[u8]{self.seq}
    fn qual(&self) -> Option<&[u8]>{self.qual}
}

impl<'a> Record for OwnedSeqRecord{
    fn head(&self) -> &[u8]{self.head.as_slice()}
    fn seq(&self) -> &[u8]{self.seq.as_slice()}
    fn qual(&self) -> Option<&[u8]>{
        match &self.qual{
            Some(q) => return Some(q.as_slice()),
            None => None,
        }
    }
}

impl<'a> SeqRecord<'a>{
    pub fn to_owned(&self) -> OwnedSeqRecord{
        OwnedSeqRecord { 
            head: self.head.to_vec(), 
            seq: self.seq.to_vec(), 
            qual: match self.qual {
                Some(q) => Some(q.to_vec()), 
                None => None
            }
        }
    }
}

impl<'a> fmt::Display for SeqRecord<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f,
               "SeqRecord{{ \n  Head: {}\n  Seq:  {}\n  Qual: {}\n}}", 
               str::from_utf8(self.head).unwrap(),
               str::from_utf8(self.seq).unwrap(),
               match self.qual{
                   Some(q) => str::from_utf8(q).unwrap(),
                   None => "", // No quality values
               }
               
        )
    }
}

// Does the FASTQ interpretation of the quality bytes.
#[inline(always)]
pub fn interpret_quality_value(q: u8) -> u8{
    q - 0x21
}

#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn test_interpret_quality_values(){
        let mut data: Vec<u8> = vec![b'!', b'#', b'%'];
        data.iter_mut().for_each(|x| *x = interpret_quality_value(*x));
        assert_eq!(data, vec![0,2,4])
    }
}