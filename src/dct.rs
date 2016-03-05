//! This module contains everything related to 
//! the discrete cosine transformation.

// Code translated from stb_image to Rust. 

macro_rules! f2f {
    ($x:expr) => {($x * 4096.0 + 0.5) as i32}
}
macro_rules! fsh {
    ($x:expr) => {$x << 12}
}

macro_rules! assign2 {
    ($x:ident) => {$x = 2}
}

// take a -128..127 value and stbi__clamp it and convert to 0..255
macro_rules! clamp {
    ($x:expr) => {
        // trick to use a single test to catch both cases
        if $x as u32 > 255 {
           if $x < 0 {0}
           else if $x > 255 {255}
           else { $x as u8 }
        } else {
            $x as u8
        }
    }
}

// derived from jidctint -- DCT_ISLOW
macro_rules! idct_2d {
    (($s0:expr, $s1:expr, $s2:expr, $s3:expr, $s4:expr, $s5:expr, $s6:expr, $s7:expr)
    -> ($t0:ident, $t1:ident, $t2:ident, $t3:ident, $x0:ident, $x1:ident, $x2:ident, $x3:ident)) => {{
        let (mut p1, mut p2, mut p3, mut p4, p5);
        p2 = $s2 as i32;
        p3 = $s6 as i32;
        p1 = (p2 + p3) * f2f!(0.5411961);
        $t2 = p1 + p3 * f2f!(-1.847759065);
        $t3 = p1 + p2 * f2f!( 0.765366865);
        p2 = $s0 as i32;
        p3 = $s4 as i32;
        $t0 = fsh!(p2+p3);
        $t1 = fsh!(p2-p3);
        $x0 = $t0 + $t3;
        $x3 = $t0 - $t3;
        $x1 = $t1 + $t2;
        $x2 = $t1 - $t2;
        $t0 = $s7 as i32;
        $t1 = $s5 as i32;
        $t2 = $s3 as i32;
        $t3 = $s1 as i32;
        p3 = $t0 + $t2;
        p4 = $t1 + $t3;
        p1 = $t0 + $t3;
        p2 = $t1 + $t2;
        p5 = (p3 + p4) * f2f!( 1.175875602);
        $t0 = $t0 * f2f!( 0.298631336);
        $t1 = $t1 * f2f!( 2.053119869);
        $t2 = $t2 * f2f!( 3.072711026);
        $t3 = $t3 * f2f!( 1.501321110);
        p1 = p5 + p1 * f2f!(-0.899976223);
        p2 = p5 + p2 * f2f!(-2.562915447);
        p3 = p3 * f2f!(-1.961570560);
        p4 = p4 * f2f!(-0.390180644);
        $t3 += p1 + p4;
        $t2 += p2 + p3;
        $t1 += p2 + p4;
        $t0 += p1 + p3;
    }}
}

pub fn idct_block(dst: &mut [u8], out_stride: usize, data: &[i16; 64]) {
   let mut v = [0i32; 64];
   // columns
   for i in 0..8 {
   //for (i=0; i < 8; ++i,++d, ++v) {
       let d = &data[i..];
       let v = &mut v[i..];
       // if all zeroes, shortcut -- this avoids dequantizing 0s and IDCTing
       if d[ 8]==0 && d[16]==0 && d[24]==0 && d[32]==0 && d[40]==0 && d[48]==0 && d[56]==0 {
           //    no shortcut                 0     seconds
           //    (1|2|3|4|5|6|7)==0          0     seconds
           //    all separate               -0.047 seconds
           //    1 && 2|3 && 4|5 && 6|7:    -0.047 seconds
           let dcterm = (d[0] as i32) << 2;
           v[ 0] = dcterm;
           v[ 8] = dcterm;
           v[16] = dcterm;
           v[24] = dcterm;
           v[32] = dcterm;
           v[40] = dcterm;
           v[48] = dcterm;
           v[56] = dcterm;
       } else {
           let (mut t0, mut t1, mut t2, mut t3, mut x0, mut x1, mut x2, mut x3);
           idct_2d!((d[ 0],d[ 8],d[16],d[24],d[32],d[40],d[48],d[56])
           -> (t0, t1, t2, t3, x0, x1, x2, x3));
           // constants scaled things up by 1<<12; let's bring them back
           // down, but keep 2 extra bits of precision
           x0 += 512; x1 += 512; x2 += 512; x3 += 512;
           v[ 0] = (x0+t3) >> 10;
           v[56] = (x0-t3) >> 10;
           v[ 8] = (x1+t2) >> 10;
           v[48] = (x1-t2) >> 10;
           v[16] = (x2+t1) >> 10;
           v[40] = (x2-t1) >> 10;
           v[24] = (x3+t0) >> 10;
           v[32] = (x3-t0) >> 10;
       }
   }
   let mut o = dst;
   let mut v = &v[..];
   for i in 0..8 {
   //for (i=0, v=val, o=out; i < 8; ++i,v+=8,o+=out_stride) {
       // no fast case since the first 1D IDCT spread components out
           let (mut t0, mut t1, mut t2, mut t3, mut x0, mut x1, mut x2, mut x3);
       idct_2d!((v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7]) -> (t0, t1, t2, t3, x0, x1, x2, x3));
       // constants scaled things up by 1<<12, plus we had 1<<2 from first
       // loop, plus horizontal and vertical each scale by sqrt(8) so together
       // we've got an extra 1<<3, so 1<<17 total we need to remove.
       // so we want to round that, which means adding 0.5 * 1<<17,
       // aka 65536. Also, we'll end up with -128 to 127 that we want
       // to encode as 0..255 by adding 128, so we'll add that before the shift
       x0 += 65536 + (128<<17);
       x1 += 65536 + (128<<17);
       x2 += 65536 + (128<<17);
       x3 += 65536 + (128<<17);
       // tried computing the shifts into temps, or'ing the temps to see
       // if any were out of range, but that was slower
       o[0] = clamp!((x0+t3) >> 17);
       o[7] = clamp!((x0-t3) >> 17);
       o[1] = clamp!((x1+t2) >> 17);
       o[6] = clamp!((x1-t2) >> 17);
       o[2] = clamp!((x2+t1) >> 17);
       o[5] = clamp!((x2-t1) >> 17);
       o[3] = clamp!((x3+t0) >> 17);
       o[4] = clamp!((x3-t0) >> 17);
       if i < 7 {
           v = &v[8..];
           let tmp = o; o = &mut tmp[out_stride..];
       }
       
   }
}