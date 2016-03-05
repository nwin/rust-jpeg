extern crate jpeg;

fn main() {
    let file = include_bytes!("andreas.jpeg");
    let mut decoder = jpeg::Decoder::new();
    decoder.update(file);
}