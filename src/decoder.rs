//! The JPEG-decoder

use std::borrow::Cow;
use std::cmp::min;
use std::fmt;
use std::ops::Range;

use dct;
use huffman::HuffmanTable;
use util::BitReader;
use Marker;

enum State {
    FindMarker,
    AtMarker,
    Byte(u16, ByteValue),
    U16Byte1(U16Value, u16),
    U16Byte0(U16Value),
    Bits(u8, u8, u8, Option<(i16, i16)>),
    Skip(u16),
}

enum ByteValue {
    Tuple(TupleValue),
    Qk(u8, u8),
    Ns,
    Csj(u8),
    Ss,
    Se,
    Li(u8, u8),
    Vij(u8, u8, u8),
    P,
    Nf,
    Ci,
    Tqi(u8)
}

enum TupleValue {
    PqTq,
    TcTh,
    TdjTaj(u8, u8),
    HiVi(u8),
    AhAl,
}

enum U16Value {
    Lq,
    Lh,
    Ls,
    Lr,
    Lf,
    Y(u16),
    X(u16),
    Ri,
    LSkip
}

macro_rules! unpack_4bit {
    ($byte:expr) => {
        ($byte >> 4, $byte & 0x0F)
    }
}

/// Data unit iterator
#[derive(Debug)]
struct DataUnits {
    c: u8,
    h: u8,
    v: u8
}

impl DataUnits {
    fn next(&mut self, scan_components: &[u8], components: &[Component]) -> Option<(u8, u8, u8)> {
        let c = &components[scan_components[self.c as usize] as usize];
        self.h = self.h.wrapping_add(1);
        if self.h >= c.h_scale {
            self.h = 0;
            self.v += 1;
            if self.v >= c.v_scale {
                self.v = 0;
                self.c += 1;
                if self.c >= scan_components.len() as u8 {
                    return None
                }
            }
        }
        Some((scan_components[self.c as usize], self.v,  self.h))
    }
}

impl Default for DataUnits {
    fn default() -> DataUnits {
        DataUnits {
            c: 0,
            v: 0,
            h: 255,
        }
    }
}


pub struct Component {
    // identifier
    id: u8,
    // horizontal scaling
    h_scale: u8,
    // vertical scaling
    v_scale: u8,
    // quantisation table
    q_table: u8,
    // dc dct table
    dc_table: u8,
    // ac dct table
    ac_table: u8,
    // previous dv value
    prev_dc: i16,
    // width of component
    pub width: u16,
    // height of component
    pub height: u16,
    // # x blocks in non-interleaved config
    pub blocks_x: u16,
    // # y blocks in non-interleaved config
    pub blocks_y: u16,
    // absolute width (#MCU * scale * 8)
    pub buffer_width: u16,
    // absolute height
    pub buffer_height: u16,
    pub blocks: Vec<[i16; 64]>,
    pub data: Vec<u8>
}

impl Component {
    fn new(id: u8, h_scale: u8, v_scale: u8) -> Component {
        Component {
            id: id,
            h_scale: h_scale,
            v_scale: v_scale,
            q_table: 0,
            dc_table: 0,
            ac_table: 0,
            prev_dc: 0,
            width: 0,
            height: 0,
            blocks_x: 0,
            blocks_y: 0,
            buffer_width: 0,
            buffer_height: 0,
            blocks: Vec::new(),
            data: Vec::new(),
        }
    }
}

pub struct Decoder {
    state: Option<State>,
    at_marker: Option<Marker>,
    pub output: Vec<u8>,
    bits: BitReader,
    unit: DataUnits,
    units: u32,
    progressive: bool,
    restart_interval: u16,
    eob_count: u16,
    pub components: Vec<Component>,
    q_tables: [[u8; 64]; 4],
    dc_tables: [HuffmanTable; 4],
    ac_tables: [HuffmanTable; 4],
    current_block: [i16; 64],
    /// Width of the image
    pub width: u16,
    /// Height of the image
    pub height: u16,
    h_scale_max: u8,
    v_scale_max: u8,
    /// # of MCUs in X direction
    mcu_x: u16, 
    /// # of MCUs in Y direction
    mcu_y: u16,
    /// scan data
    spectral_select: Range<u8>,
    succ_approx: Range<u8>,
    scan_components: Vec<u8>,
}

impl fmt::Debug for Decoder {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        try!(writeln!(fmt, "Decoder {{"));
        try!(writeln!(fmt, "    q_tables: ["));
        for table in self.q_tables.iter() {
            try!(write!(fmt, "        ["));
            for line in table.as_ref().chunks(8) {
                try!(write!(fmt, "\n         "));
                for val in line {
                    try!(write!(fmt, "{:X}, ", val));
                }
            }
            try!(writeln!(fmt, "],"));
        }
        try!(writeln!(fmt, "    ]"));
        write!(fmt, "}}")
    }
}

macro_rules! current_block {
    ($this:expr, $c:expr) => ({
        if $this.progressive {
            if $this.scan_components.len() > 1 {
                let x = ($this.units % $this.mcu_x as u32) as u16;
                let y = ($this.units / $this.mcu_x as u32) as u16;
                let x2 = (x * $c.h_scale as u16 + $this.unit.h as u16) as usize;
                let y2 = (y * $c.v_scale as u16 + $this.unit.v as u16) as usize;
                &mut $c.blocks[$c.buffer_width as usize / 8 * y2 + x2]
            } else {
                let x = ($this.units % $c.blocks_x as u32) as usize;
                let y = ($this.units / $c.blocks_x as u32) as usize;
                &mut $c.blocks[$c.buffer_width as usize / 8 * y + x]
            }
        } else {
            &mut $this.current_block
        }
    })
}

macro_rules! next_unit {
    ($this:expr) => (
        $this.unit.next(&$this.scan_components, &$this.components)
    )
}

impl Decoder {
    pub fn new() -> Decoder {
        Decoder {
            state: Some(State::FindMarker),
            at_marker: None,
            output: Vec::new(),
            bits: BitReader::default(),
            unit: DataUnits::default(),
            units: 0,
            progressive: false,
            restart_interval: 0,
            eob_count: 0,
            components: Vec::with_capacity(4), // Max 4 components
            q_tables: [[0; 64]; 4],
            width: 0,
            height: 0,
            dc_tables: [HuffmanTable::default(), HuffmanTable::default(), HuffmanTable::default(), HuffmanTable::default()],
            ac_tables: [HuffmanTable::default(), HuffmanTable::default(), HuffmanTable::default(), HuffmanTable::default()],
            current_block: [0; 64],
            h_scale_max: 0,
            v_scale_max: 0,
            mcu_x: 0,
            mcu_y: 0,
            spectral_select: Range { start: 0, end: 0 },
            succ_approx: Range { start: 0, end: 0 },
            scan_components: Vec::with_capacity(4) // Max 4 components
        }
    }
    
    pub fn update(&mut self, mut buf: &[u8]) {
        while buf.len() > 0 {
            buf = &buf[self.next_state(buf).unwrap()..]
        }
    }
    
    fn next_state(&mut self, buf: &[u8]) -> Result<usize, Cow<'static, str>> {
        use self::State::*;
        macro_rules! goto {
            ($state:expr) => ({
                self.state = Some($state);
                Ok(1)
            });
            ($n:expr, $state:expr) => ({
                self.state = Some($state);
                Ok($n)
            })
        }
        macro_rules! need_byte {
            ($val:ident, $params:expr) => {
            }
        }
        let byte = buf[0];
        match self.state.take().unwrap() {
            Bits(c_idx, mut k, eob_ss_start, s_r) => {
                let total_units = if self.interleaved() {
                    self.mcu_x as u32 * self.mcu_y as u32
                } else {
                    let c = &self.components[c_idx as usize];
                    c.blocks_x as u32 * c.blocks_y as u32
                };
                let mut n = 0;
                loop {                    
                    n += {
                        let bits = try!(self.fill_bits(&buf[n..]));
                        if bits == 0 && self.bits.bits == 0 && self.eob_count == 0 {
                            break
                            //panic!("unexpected eof")
                        }
                        bits
                    };
                    if self.bits.bits < 32 && self.at_marker.is_none() {
                        return goto!(n, Bits(c_idx, k, eob_ss_start, None))
                    }
                    if self.succ_approx.end != 0 {
                        let delta = 1i16 << self.succ_approx.start;
                        if self.eob_count > 0 {
                            println!("eob_count {}", self.eob_count);
                            let c = &mut self.components[c_idx as usize];
                            let block = &mut current_block!(self, c);
                            println!("{:?}", &block[..]);
                            for i in eob_ss_start..self.spectral_select.end {
                                let p = &mut block[ZIGZAG[i as usize] as usize];
                                if *p != 0 {
                                    if self.bits.bits == 0 {
                                        println!("refill buf, len: {}, bits: {}", buf.len(), self.bits.bits);
                                        println!("marker: {:?}", self.at_marker);
                                        return goto!(n, Bits(c_idx, k, i, None))
                                    }
                                    if try!(self.bits.read(1)) != 0 {
                                        if *p & delta == 0 {
                                            if *p > 0 {
                                                *p += delta
                                            } else {
                                                *p -= delta
                                            }
                                            println!("refined non-0-delta-symbol 0x{:x}", *p);
                                        }
                                    }
                                }
                            }
                            self.eob_count -= 1;
                            break
                        }
                        match k {
                            0 => {
                                let c = &mut self.components[c_idx as usize];
                                if try!(self.bits.read(1)) != 0 {
                                    current_block!(self, c)[0] += 1i16 << self.succ_approx.start;
                                    println!("DC offset 0x{:02X}", 1i16 << self.succ_approx.start);
                                    println!("refined DC coef. 0x{:02X}", current_block!(self, c)[0]);
                                }
                                k = 1;
                            }
                            _ if k <= self.spectral_select.end => {
                                let c = &mut self.components[c_idx as usize];
                                let (s, mut r) = if let Some((s, r)) = s_r {
                                    (s, r)
                                } else {
                                    let symbol = try!(self.ac_tables[c.ac_table as usize]
                                                     .decode(&mut self.bits)) as i16;
                                    println!("decoded symbol 0x{:x}", symbol);
                                    let mut s = symbol % 16;
                                    let mut r = symbol >> 4;
                                    match s {
                                        0 => {
                                            if r < 15 {
                                                self.eob_count = (1u16 << r) - 1;
                                                if r > 0 {
                                                    self.eob_count += try!(self.bits.read(r as u8));
                                                }
                                                r = 64;
                                            }
                                        }
                                        1 => {
                                            if try!(self.bits.read(1)) != 0 {
                                                s = delta as i16
                                            } else {
                                                s = (-delta) as i16
                                            }
                                        }
                                        _ => return Err("invalid Huffman code".into())
                                    }
                                    (s, r)
                                };
                                let block = &mut current_block!(self, c);
                                println!("{:?}", &block[..]);
                                println!("k {:?}", k);
                                println!("se {:?}", self.spectral_select.end);
                                println!("bits {:?}", self.bits.bits);
                                while k <= self.spectral_select.end {
                                    if self.bits.bits == 0 {
                                        return goto!(n, Bits(c_idx, k, self.succ_approx.start, Some((r, s))))
                                    }
                                    let p = &mut block[ZIGZAG[k as usize] as usize];
                                    k += 1;
                                    if *p != 0 {
                                        if try!(self.bits.read(1)) != 0 {
                                            if *p & delta == 0 {
                                                if *p > 0 {
                                                    *p += delta
                                                } else {
                                                    *p -= delta
                                                }
                                                println!("refined 0-coef. to {}", *p);
                                            }
                                        }
                                    } else {
                                        if r == 0 {
                                            *p = s as i16;
                                            println!("refined coef. to {}", *p);
                                            break
                                        }
                                        r -= 1;
                                    }
                                }
                                if self.eob_count > 0 {
                                    break
                                }
                            }
                            _ => {
                                break
                            }
                        }
                        continue;
                    }
                    if self.eob_count > 0 {
                        println!("eob_count {}", self.eob_count);
                        self.eob_count -= 1;
                        break // maybe
                    }
                    match k {
                        0 => {
                            let c = &mut self.components[c_idx as usize];
                            let symbol = try!(self.dc_tables[c.dc_table as usize].decode(&mut self.bits));
                            let v = (try!(self.bits.read(symbol)) as i32).extend(symbol) as i16;
                            c.prev_dc += v;
                            current_block!(self, c)[0] = c.prev_dc << self.succ_approx.start;
                            println!("Decoded DC coef. {}", current_block!(self, c)[0]);
                            k = 1;
                            continue
                        }
                        _ if k <= self.spectral_select.end => {
                            println!("bits {}, {:?}", self.bits.bits, self.at_marker);
                            let c = &mut self.components[c_idx as usize];
                            let symbol = try!(self.ac_tables[c.ac_table as usize]
                                             .decode(&mut self.bits));
                            let s = symbol % 16;
                            let r = symbol >> 4;
                            if s > 0 {
                                k += r;
                                let idx = ZIGZAG[k as usize] as usize;
                                current_block!(self, c)[idx] = 
                                    ((try!(self.bits.read(s)) as i32).extend(s) as i16)
                                    << self.succ_approx.start;
                                println!("Decoded AC coef. 0x{:02X}", current_block!(self, c)[idx]);
                                if k == self.spectral_select.end {
                                    break
                                } else if k > 63 {
                                    return Err("Invalid block.".into())
                                } else {
                                    k += 1;
                                    continue
                                }
                            } else {
                                if r == 15 { // Zero run length
                                    k += 16;
                                    continue
                                } else { // End of block/band
                                    self.eob_count = (1u16 << r) - 1;
                                    if r > 0 {
                                        self.eob_count += try!(self.bits.read(r));
                                    } 
                                    break
                                }
                            }
                        }
                        _ => {
                            break
                        }
                    }
                }
                if !self.progressive {
                    self.transform_block();
                }
                match if self.scan_components.len() == 1 { None } else { next_unit!(self) } {
                    Some((c_idx, _, _)) => {
                        if !self.progressive {
                            self.current_block = [0; 64]
                        }
                        goto!(n, Bits(c_idx, self.spectral_select.start, self.spectral_select.start, None))
                    },
                    None => {
                        println!("unit {} of {} finished", self.units, total_units);
                        println!("marker {:?} eob: {}", self.at_marker, self.eob_count);
                        if !self.progressive {
                            self.current_block = [0; 64]
                        }
                        self.units += 1;
                        self.unit = DataUnits::default();
                        n += match self.at_marker {
                            _ if self.bits.bits >= 8 => 0,
                            _ if self.eob_count > 0 => 0,
                            Some(Marker::EOI) => {
                                return goto!(n, FindMarker)
                            },
                            _ if self.units > total_units => {
                                // TODO set error flag
                                return goto!(n, FindMarker)
                            },
                            Some(Marker::RSTm(_)) => {
                                if self.units % self.restart_interval as u32 == 0 {
                                    self.bits = BitReader::default();
                                    self.at_marker = None;
                                    for c in &mut self.components {
                                        c.prev_dc = 0
                                    }
                                    2 // consume RST marker
                                } else {
                                    return Err("Unexpected restart marker.".into())
                                }
                            },
                            Some(_) if self.units == total_units => {
                                return goto!(n, FindMarker)
                            },
                            // TODO: Better handling of this
                            Some(_) if self.progressive => 0,
                            Some(ref m) => {
                            println!("k {}", k);
                            println!("bits {}", self.bits.bits);
                            println!("{} == {} {}", self.units , total_units, self.eob_count);
                            return Err(
                                format!("Unexpected marker {:?} found.", m).into()
                            )},
                            None => if self.restart_interval > 0 
                                    && self.units % self.restart_interval as u32 == 0 {
                                return Err("Restart marker missing.".into())
                            } else {
                                0 // consume nothing
                            }
                        };
                        let c_idx = if self.scan_components.len() == 1 {
                            self.scan_components[0]
                        } else {
                            next_unit!(self).unwrap().0
                        };
                        goto!(n, Bits(c_idx, self.spectral_select.start, self.spectral_select.start, None))
                    }    
                }
            }
            FindMarker => match byte {
                0xFF => goto!(AtMarker),
                _ => goto!(FindMarker),
            },
            Byte(left, byte_val) => {
                use self::ByteValue::*;
                use self::U16Value::*;
                use self::TupleValue::*;
                match byte_val {
                    Tuple(tuple_value) => {
                        let val = unpack_4bit!(byte);
                        match tuple_value {
                            PqTq => {
                                if left == 0 {
                                    return goto!(0, FindMarker)
                                }
                                if val.0 > 0 {
                                    return Err("Only baseline jpeg is supported.".into())
                                }
                                if val.1 > 3 {
                                    return Err(
                                        "Invalid quantization table destination identifier.".into()
                                    )
                                }
                                goto!(Byte(left - 1, Qk(val.1, 0)))
                            }
                            TcTh => {
                                if val.0 > 1 {
                                    return Err("Invalid Huffman table class.".into())
                                }
                                if val.1 > 3 {
                                    return Err(
                                        "Invalid Huffman table destination identifier.".into()
                                    )
                                }
                                goto!(Byte(left - 1, Li(byte, 0)))
                            }
                            TdjTaj(n, idx) => {
                                if val.0 > 3 || val.1 > 3 {
                                    return Err(
                                        "Invalid Huffman table destination identifier.".into()
                                    )
                                }
                                self.components[idx as usize].dc_table = val.0;
                                self.components[idx as usize].ac_table = val.1;
                                let n = n + 1;
                                if n < self.scan_components.len() as u8 {
                                    goto!(Byte(left - 1, Csj(n)))
                                } else {
                                    goto!(Byte(left - 1, Ss))
                                }
                            }
                            HiVi(id) => {
                                match val {
                                    (1...4, 1...4) => (),
                                    _ => return Err("Invalid scaling factors.".into())
                                }
                                self.components.push(Component::new(id, val.0, val.1));
                                goto!(Byte(left - 1, Tqi(self.components.len() as u8 - 1)))
                            },
                            AhAl => {
                                if left != 1 {
                                    return Err("Scan header too long.".into())
                                }
                                if !self.progressive && (val.0 != 0 || val.1 != 0 ) {
                                    return Err("Invalid succ. approx. bits.".into())
                                }
                                if val.0 > 13 || val.1 > 13 {
                                    return Err("Invalid succ. approx. bits.".into())
                                }
                                self.succ_approx.end = val.0;
                                self.succ_approx.start = val.1;
                                self.unit = DataUnits::default();
                                let c_idx = if self.scan_components.len() == 1 {
                                    self.scan_components[0]
                                } else {
                                    next_unit!(self).unwrap().0
                                };
                                goto!(Bits(c_idx, self.spectral_select.start, self.spectral_select.start, None))
                            }
                        }
                    }
                    Qk(idx, n) => {
                        self.q_tables[idx as usize][ZIGZAG[n as usize] as usize] = byte;
                        if n < 63 {
                            goto!(Byte(left - 1, Qk(idx, n + 1)))
                        } else {
                            goto!(Byte(left - 1, Tuple(PqTq)))
                        }
                    }
                    Ns => {
                        if byte != self.components.len() as u8 && byte != 1 {
                            return Err("Unsupported number of components in scan.".into())
                        }
                        self.scan_components.clear();
                        for _ in 0..byte {
                            self.scan_components.push(0)
                        }
                        goto!(Byte(left - 1, Csj(0)))
                    }
                    Csj(n) => {
                        let mut idx = None;
                        for (i, c) in self.components.iter().enumerate() {
                            if c.id == byte {
                                idx = Some(i);
                                break
                            }
                        }
                        match idx {
                            Some(idx) => {
                                self.scan_components[n as usize] = idx as u8;
                                goto!(Byte(left - 1, Tuple(TdjTaj(n, idx as u8))))
                            },
                            None => Err("Unknown component id.".into())
                        }
                    }
                    Ss => {
                        if byte != 0 && !self.progressive {
                            return Err("Sequential DCT only allows Ss = 0.".into())
                        } if byte > 63 {
                            return Err("Invalid Ss value.".into())
                        }
                        self.spectral_select.start = byte;
                        goto!(Byte(left - 1, Se))
                    }
                    Se => {
                        if byte != 63 && !self.progressive {
                            return Err("Sequential DCT only allows Se = 63.".into())
                        } else if byte > 63 || byte < self.spectral_select.start {
                            return Err("Invalid Se value.".into())
                        }
                        self.spectral_select.end = byte;
                        goto!(Byte(left - 1, Tuple(AhAl)))
                    }
                    Li(tcth, n) => {
                        let (tc, th) = unpack_4bit!(tcth);
                        let table = &mut match tc {
                            0 => &mut self.dc_tables,
                            1 => &mut self.ac_tables,
                            _ => unreachable!()
                        }[th as usize];
                        if n == 0 {
                            // Reset Huffman table
                            *table = HuffmanTable::default()
                        }
                        table.size[n as usize] = byte;
                        if n < 15 {
                            goto!(Byte(left - 1, Li(tcth, n + 1)))
                        } else {
                            let total = table.size
                                .iter().take(16)
                                .fold(0u16, |acc, &item| acc + item as u16);
                            if total > 256 {
                                return Err("Too many Huffman codes".into())
                            }
                            goto!(Byte(left - 1, Vij(tcth, 0, total as u8)))
                        }
                    }
                    Vij(tcth, n, total) => {
                        let (tc, th) = unpack_4bit!(tcth);
                        let table = &mut match tc {
                            0 => &mut self.dc_tables,
                            1 => &mut self.ac_tables,
                            _ => unreachable!()
                        }[th as usize];
                        table.val[n as usize] = byte;
                        let n = n + 1;
                        if n < total {
                            goto!(Byte(left - 1, Vij(tcth, n, total)))
                        } else {
                            table.build();
                            if left > 1 {
                                goto!(Byte(left - 1, Tuple(TcTh)))
                            } else {
                                goto!(FindMarker)
                            }
                        }
                    }
                    P => goto!(U16Byte0(Y(left - 1))),
                    Nf => {
                        let left = left - 1;
                        if left != 3 * byte as u16 {
                            return Err("Invalid frame header".into())
                        }
                        self.components = Vec::with_capacity(byte as usize);
                        goto!(Byte(left, Ci))
                    },
                    Ci => {
                        goto!(Byte(left - 1, Tuple(HiVi(byte))))
                    },
                    Tqi(idx) => {
                        self.components[idx as usize].q_table = byte;
                        let idx = idx + 1;
                        if idx < self.components.capacity() as u8 {
                            goto!(Byte(left - 1, Ci))
                        } else {
                            self.init_frame();
                            goto!(FindMarker)
                        }
                    }
                }
            },
            U16Byte1(u16_val, mut value) => {
                use self::U16Value::*;
                use self::ByteValue::*;
                use self::TupleValue::*;
                value |= byte as u16;
                match u16_val {
                    Lq => goto!(Byte(value - 2, Tuple(PqTq))),
                    Lh => goto!(Byte(value - 2, Tuple(TcTh))),
                    Ls => {
                        self.at_marker = None;
                        self.bits = BitReader::default();
                        self.units = 0;
                        goto!(Byte(value - 2, Ns))
                    },
                    Lf => goto!(Byte(value - 2, P)),
                    Y(left) => {
                        self.height = value;
                        goto!(U16Byte0(X(left - 2)))
                    },
                    X(left) => {
                        if value > 0 {
                            self.width = value;
                            goto!(Byte(left - 2, Nf))
                        } else {
                            Err("Line width = 0.".into())
                        }
                    },
                    Lr => goto!(U16Byte0(Ri)),
                    Ri => {
                        self.restart_interval = value;
                        goto!(FindMarker)
                    },
                    LSkip => goto!(Skip(value - 2)),
                }
            }
            U16Byte0(u16_val) => goto!(U16Byte1(u16_val, (byte as u16) << 8)),
            AtMarker => match byte {
                // Stuff byte
                0x00 => goto!(FindMarker),
                // Fill byte
                0xFF => goto!(AtMarker),
                n => {
                    match Marker::from_u8(n) {
                        Some(marker) => {
                            use ::Marker::*;
                            println!("*** marker: {:?}", marker);
                            match marker {
                                SOI => {
                                    self.restart_interval = 0;
                                    goto!(FindMarker)
                                }
                                DQT => goto!(U16Byte0(U16Value::Lq)),
                                DHT => goto!(U16Byte0(U16Value::Lh)),
                                DRI => goto!(U16Byte0(U16Value::Lr)),
                                SOS => goto!(U16Byte0(U16Value::Ls)),
                                SOF(0) => goto!(U16Byte0(U16Value::Lf)),
                                SOF(2) => {
                                    self.progressive = true;
                                    goto!(U16Byte0(U16Value::Lf))
                                },
                                SOF(_) => return Err(
                                    "Only Huffman coded baseline + progressive JPEG is supported.".into()
                                ),
                                EOI => {
                                    let mut output = vec![ 0;
                                        self.width as usize
                                        * self.height as usize
                                        * self.components.len()
                                    ];
                                    if self.progressive {
                                        self.transform_blocks()
                                    }
                                    println!("#mcus {}, bits {}", self.units, self.bits.bits);
                                    self.resample(&mut output);
                                    for ycbcr in output.chunks_mut(3) {
                                        ycbcr_to_rgb(ycbcr);
                                    }
                                    self.output = output;
                                    goto!(FindMarker)
                                }
                                _ => goto!(U16Byte0(U16Value::LSkip)),
                            }
                        }
                        None => panic!("Unknown marker {:X}", n)
                    }
                }
            },
            Skip(left) => {
                let to_skip = min(left, buf.len() as u16);
                goto!(to_skip as usize, match left - to_skip {
                    0 => FindMarker,
                    n => Skip(n)
                })
            }
        }
    }
    
    fn fill_bits(&mut self, mut buf: &[u8]) -> Result<usize, Cow<'static, str>> {
        if self.at_marker.is_some() {
            return Ok(0)
        }
        let mut n = 0;
        while self.bits.bits <= 56 && buf.len() > 1 {
            let byte = buf[0];
            if byte == 0xFF {
                match buf.get(1) {
                    Some(&byte) => if byte == 0 {
                         self.bits.push(0xFF);
                         buf = &buf[2..];
                         n += 2;
                    } else {
                        self.at_marker = match Marker::from_u8(byte) {
                            None => return Err("Unknown marker found in bit stream.".into()),
                            m => m
                        };
                        //n += 2;
                        break
                    },
                    None => break
                }
            } else {
                self.bits.push(byte);
                buf = &buf[1..];
                n += 1;
            }
        }
        Ok(n)
    }
    
    fn init_frame(&mut self) {
        for c in &self.components {
            if c.h_scale > self.h_scale_max {
                self.h_scale_max = c.h_scale;
            }
            if c.v_scale > self.v_scale_max {
                self.v_scale_max = c.v_scale;
            }
        }
        let mcu_h = self.h_scale_max as u32 * 8;
        let mcu_v = self.v_scale_max as u32 * 8;
        self.mcu_x = ((self.width as u32 + mcu_h - 1) / mcu_h) as u16;
        self.mcu_y = ((self.height as u32 + mcu_v - 1) / mcu_v) as u16;
        
        for c in &mut self.components {
            c.width = (
                (self.width as u32 * c.h_scale as u32  + self.h_scale_max as u32 - 1)
                /self.h_scale_max as u32
            ) as u16;
            c.height = (
                (self.height as u32 * c.v_scale as u32  + self.v_scale_max as u32 - 1)
                /self.v_scale_max as u32
            ) as u16;
            c.blocks_x = (c.width + 7) / 8;
            c.blocks_y = (c.height + 7) / 8;
            c.buffer_width = self.mcu_x * c.h_scale as u16 * 8;
            c.buffer_height = self.mcu_y * c.v_scale as u16 * 8;
            c.data = vec![0; c.buffer_width as usize * c.buffer_height as usize/* + c.buffer_width as usize*/];
            if self.progressive {
                //c.blocks = vec![[0; 64]; c.data.len()/16]
                c.blocks = Vec::with_capacity(c.data.len()/64);
                for _ in 0..c.data.len()/64 {
                    c.blocks.push([0; 64])
                }
            }
        }
    }
    
    // Progressive data
    fn transform_blocks(&mut self) {
        let total_units = self.mcu_x as u32 * self.mcu_y as u32;
        let mut unit = DataUnits::default();
        let components = (0..self.components.len()).map(|v| v as u8).collect::<Vec<_>>();
        for units in 0..total_units {
            loop {
                match unit.next(&*components, &*self.components) {
                    Some((_, _, _)) => self.transform_block_progressive(
                        &*components,
                        units,
                        &unit
                    ),
                    None => {
                        unit = DataUnits::default();
                        break
                    }
                }
            }
        }
    }
    
    /// Transforms progressive image data
    ///
    /// There is a significant code duplication in `transform_block_progressive`
    /// a separate function `transform_block` is needed to circumvent borrowing issues. 
    fn transform_block_progressive(
        &mut self,
        scan_components: &[u8],
        units: u32,
        unit: &DataUnits
    ) {
        let c = &mut self.components[scan_components[unit.c as usize] as usize];
        let q_table = &self.q_tables[c.q_table as usize];
        let x = (units % self.mcu_x as u32) as u16;
        let y = (units / self.mcu_x as u32) as u16;
        let x2 = (x * c.h_scale as u16 + unit.h as u16) * 8;
        let y2 = (y * c.v_scale as u16 + unit.v as u16) * 8;
        let mut block = c.blocks[c.buffer_width as usize / 8 * y2 as usize / 8 + x2 as usize / 8];
        for i in 0..64 {
            block[i] *= q_table[i] as i16;
        }
        dct::idct_block(
            &mut c.data[
                c.buffer_width as usize
                * y2 as usize
                + x2 as usize
            ..],
            c.buffer_width as usize,
            &block
        );
    }
    
    /// Transforms interleaved non-progressive data.
    ///
    /// There is a significant code duplication in `transform_block_progressive`
    /// a separate function `transform_block` is needed to circumvent borrowing issues. 
    fn transform_block(&mut self) {
        let c = &mut self.components[self.scan_components[self.unit.c as usize] as usize];
        let q_table = &self.q_tables[c.q_table as usize];
        for i in 0..64 {
            self.current_block[i] *= q_table[i] as i16;
        }
        let x = (self.units % self.mcu_x as u32) as u16;
        let y = (self.units / self.mcu_x as u32) as u16;
        let x2 = (x * c.h_scale as u16 + self.unit.h as u16) * 8;
        let y2 = (y * c.v_scale as u16 + self.unit.v as u16) * 8;
        dct::idct_block(
            &mut c.data[
                c.buffer_width as usize
                * y2 as usize
                + x2 as usize
            ..],
            c.buffer_width as usize,
            &self.current_block
        );
    }
    
    /// Returns true if data is intersleaved
    #[inline(always)]
    fn interleaved(&self) -> bool {
        self.scan_components.len() == self.components.len()
    }
    
    fn resample(&mut self, output: &mut [u8]) {
        let n = self.components.len();
        let width = n * self.width as usize;
        for (i, c) in self.components.iter().enumerate() {
            let h_scale = (self.h_scale_max / c.h_scale) as usize;
            let v_scale = (self.v_scale_max / c.v_scale) as usize;
            for y in 0..c.height as usize {
                for x in 0..c.width as usize {
                    let val = c.data[y*c.buffer_width as usize + x];
                    for y1 in 0..v_scale as usize {
                        for x1 in 0..h_scale as usize {
                            match output.get_mut(y*width*v_scale + y1*width + h_scale*x*n + x1*n + i) {
                                Some(out) => *out = val,
                                None => (), // this should fail on the edges only
                            }
                        }
                    }
                }
            }
        }
    } 
}


// take a -128..127 value and stbi__clamp it and convert to 0..255
macro_rules! clamp {
    ($x:expr) => {
        // trick to use a single test to catch both cases
        if $x < 00. { 0 }
        else if $x > 255.0 { 255 }
        else { $x as u8 }
    }
}

fn ycbcr_to_rgb(ycbcr: &mut [u8]) {
    let  y = ycbcr[0] as f32;
    let cb = ycbcr[1] as f32;
    let cr = ycbcr[2] as f32;
    let r = y as f32 + 1.402   * (cr-128.0);
    let g = y as f32 - 0.34414 * (cb-128.0) - 0.71414*(cr-128.0);
    let b = y as f32 + 1.772   * (cb-128.0);
    let rgb = ycbcr;
    rgb[0] = clamp!(r);
    rgb[1] = clamp!(g);
    rgb[2] = clamp!(b);
}


const ZIGZAG: [u8; 64] = [
    0,  1,  8, 16,  9,  2,  3, 10,
    17, 24, 32, 25, 18, 11,  4,  5,
    12, 19, 26, 33, 40, 48, 41, 34,
    27, 20, 13,  6,  7, 14, 21, 28,
    35, 42, 49, 56, 57, 50, 43, 36,
    29, 22, 15, 23, 30, 37, 44, 51,
    58, 59, 52, 45, 38, 31, 39, 46,
    53, 60, 61, 54, 47, 55, 62, 63,
];

trait ExtendDiff {
    // Figure F.12: ISO/IEC 10918-1 : 1993(E)
    fn extend(self, t: u8) -> Self;
}
impl ExtendDiff for i32 {
    fn extend(mut self, t: u8) -> Self {
        if t == 0 {
            0
        } else {
            let mut vt = 1 << (t-1);
            while self < vt {
                vt = ((-1) << t) + 1;
                self += vt
            }
            self
        }
    } 
}