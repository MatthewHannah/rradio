use std::fs::File;
use std::io::{BufReader, Read};
use num_complex::{c32, Complex32};

#[derive(Debug)]
pub enum SigmfError {
    BadFile(String),
}


// definitely needs to be not hardcoded to ci32_le
// and to look at metadata
pub struct SigmfStreamer {
    sample_file : BufReader<File>
}

impl SigmfStreamer {
    pub fn new(path: &str) -> Result<SigmfStreamer,SigmfError> {
        let f = File::open(path).map_err(|e| SigmfError::BadFile(e.to_string()))?;
        Ok(SigmfStreamer {
            sample_file: BufReader::new(f)
        })
    }
}

impl Iterator for SigmfStreamer {
    type Item = Complex32;

    fn next(&mut self) -> Option<Self::Item> {
        let mut real: [u8; 4] = [0; 4];
        let mut imag: [u8; 4] = [0; 4];

        self.sample_file.read_exact(&mut real).ok()?;
        self.sample_file.read_exact(&mut imag).ok()?;
        let real = (i32::from_le_bytes(real) as f32) / 2147483648.0;
        let imag = (i32::from_le_bytes(imag) as f32) / 2147483648.0;
        Some(c32(real, imag))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn try_stream_file() -> Result<(), SigmfError> {
        let file = SigmfStreamer::new("./res/fm_radio_20250920_6msps.sigmf-data")?;
        let sum: Complex32 = file.take(100).sum();

        assert_ne!(sum, c32(0.0,0.0));

        Ok(())
    }
}