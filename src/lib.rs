//! JPEG decoder
//!

// Credits go to stb_image and jpeg.go for being a reference implementation
// and to stb_image for lending some code.

mod dct;
mod decoder;
mod huffman;
mod util;

pub use decoder::Decoder;

#[derive(Debug)]
/// Available JPEG marker
enum Marker {
// Start Of Frame markers
    /// Start Of Frame markers
    ///  - SOF(0): Baseline DCT
    SOF(u8),
    /// Reserved for JPEG extensions
    JPG,
// Huffman table specification
    /// Define Huffman table(s)
    DHT,
// Arithmetic coding conditioning specification
    /// Define arithmetic coding conditioning(s)
    DAC,
// Restart interval termination
    /// Restart with modulo 8 count `m`
    RSTm(u8),
// Other markers
    /// Start of image
    SOI,
    /// End of image
    EOI,
    /// Start of scan
    SOS,
    /// Define quantization table(s)
    DQT,
    /// Define number of lines
    DNL,
    /// Define restart interval
    DRI,
    /// Define hierarchical progression
    DHP,
    /// Expand reference component(s)
    EXP,
    /// Reserved for application segments
    APPn(u8),
    /// Reserved for JPEG extensions
    JPGn(u8),
    /// Comment
    COM,
// Reserved markers
    /// For temporary private use in arithmetic coding
    TEM,
    /// Reserved
    RES,
    
}

impl Marker {
    pub fn from_u8(n: u8) -> Option<Marker> {
        use self::Marker::*;
        Some(match n {
            // Start Of Frame markers
            0xC8 => JPG,
            0xC0 ... 0xC3 | 0xC5 ... 0xCB | 0xCD ... 0xCF  => SOF(n - 0xC0),
            // Arithmetic coding conditioning specification
            0xCC => DAC,
            // Restart interval termination
            0xD0 ... 0xD7 => RSTm(n - 0xD0),
            // Huffman table specification
            0xC4 => DHT,
            // Other markers
            0xD8 => SOI,
            0xD9 => EOI,
            0xDA => SOS,
            0xDB => DQT,
            0xDC => DNL,
            0xDD => DRI,
            0xDE => DHP,
            0xDF => EXP,
            0xE0 ... 0xEF => APPn(n - 0xE0),
            0xF0 ... 0xFD => JPGn(n - 0xF0),
            0xFE => COM,
            0x01 => TEM,
            0x02 ... 0xBF => RES,
            _ => return None
        })
    }
}