//! This module contains everything related to 
//! the entropy coding.
use std::borrow::Cow;
use util::BitReader;

const FAST_BITS: usize = 9;
const BIT_MASK: [u32; 17] = [
    0,1,3,7,15,31,63,127,255,511,1023,2047,4095,8191,16383,32767,65535
];

pub struct HuffmanTable {
    fast: [u8; 1 << FAST_BITS],
    /// HUFFCODE in spec
    code: [u16; 256],
    /// HUFFVAL in spec
    pub val: [u8; 256],
    /// HUFFSIZE in spec
    pub size: [u8; 257],
    /// MAXCODE in spec
    max: [u32; 17],
    delta: [i32; 18]
}

impl HuffmanTable {
    // Rebuilds the decoding tables from the table definitions
    // found in the JPEG stream.
    //
    // For this purpose this method assumes that the sizes have bee
    // put into the first 16 elements of `Self::size` and the Huffman
    // values into `Self::val`. In the order they occured in the JPEG
    // stream.
    //
    // The accelleration code was copied from stb_image
    pub fn build(&mut self) {
        let mut bits: [u8; 16] = [0; 16];
        // copy sizes out of the size table
        for i in 0..16 { bits[i] = self.size[i] }
        let bits = bits;
        // Generate_size_table
        // Figure C.1 from ISO/IEC 10918-1
        let mut k = 0;
        for i in 0..16 {
            for _ in 0..bits[i] {
                self.size[k] = i as u8 + 1;
                k += 1;
            }
        }
        self.size[k] = 0;
        
        // generate symbols
        let mut code = 0u32;
        k = 0;
        for j in 1..17 {
            // compute delta to add to code to compute symbol id
            self.delta[j] = k as i32 - code as i32;
            if self.size[k] == j as u8 {
                while self.size[k] == j as u8 {
                    self.code[k] = code as u16;
                    k += 1;
                    code += 1;
                }
            }
            // compute largest code + 1 for this size, preshifted as needed later
            self.max[j] = code << (16-j);
            code <<= 1
        }
        self.max[16] = 0xFFFF_FFFF;
        
        // build non-spec acceleration table; 255 is flag for not-accelerated
        for i in 0..k {
            let s = self.size[i];
            if s <= FAST_BITS as u8 {
                let c = self.code[i] << (FAST_BITS - s as usize);
                let m = 1 << (FAST_BITS - s as usize);
                for j in 0..m {
                    self.fast[c as usize + j] = i as u8;
                }
            }
        }   
    }
    
    pub fn decode(&self, bits: &mut BitReader) -> Result<u8, Cow<'static, str>> {
        let c = (bits.acc >> (64 - FAST_BITS)) & ((1 << FAST_BITS)-1);
        let k = self.fast[c as usize];
        Ok(if k < 0xFF {
            let s = self.size[k as usize];
            if s > bits.bits {
                return Err("Not enough bits to decode symbol.".into())
            }
            bits.acc <<= s;
            bits.bits -= s;
            self.val[k as usize]
        } else {
            let tmp = (bits.acc >> 48) as u32;
            let mut k = FAST_BITS + 1;
            loop {
                if tmp < self.max[k] {
                    break
                }
                k += 1;
            }
            if k == 17 {
                return Err("Huffman code not found.".into())
            }
            if k > bits.bits as usize {
                return Err("Not enough bits to decode symbol.".into())
            }
            let c = ((bits.acc >> (64 - k)) as u32 & BIT_MASK[k]) as i32 + self.delta[k];
            bits.acc <<= k;
            bits.bits -= k as u8;
            self.val[c as usize]
        })
    }
}

impl Default for HuffmanTable {
    fn default() -> HuffmanTable {
        HuffmanTable {
            fast: [0xFF; 1 << FAST_BITS],
            code: [0; 256],
            val: [0; 256],
            size: [0; 257],
            max: [0; 17],
            delta: [0; 18]
        }
    }
}