//! This module contains utility function 
use std::borrow::Cow;

pub struct BitReader {
    pub acc: u64,
    pub bits: u8
}

impl BitReader {
    pub fn push(&mut self, byte: u8) {
        assert!(self.bits <= 56);
        self.acc |= (byte as u64) << (56 - self.bits);
        self.bits += 8;
    }
    
    pub fn read(&mut self, n: u8)  -> Result<u16, Cow<'static, str>> {
        assert!(n <= 16);
        if n <= self.bits {
            if n == 0 { return Ok(0) }
            let res = self.acc >> (64 - n);
            self.acc <<= n;
            self.bits -= n;
            Ok(res as u16)
        } else {
            Err("Not enough bits".into())
        }
    }
}

impl Default for BitReader {
    fn default() -> BitReader {
        BitReader {
            acc: 0,
            bits: 0
        }
    }
}