use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;
use num_complex::{c32, Complex32};
use serde::Deserialize;

#[derive(Debug)]
pub enum SigmfError {
    BadFile(String),
    UnsupportedDatatype(String),
}

impl std::fmt::Display for SigmfError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SigmfError::BadFile(msg) => write!(f, "bad file: {}", msg),
            SigmfError::UnsupportedDatatype(dt) => write!(f, "unsupported datatype: {}", dt),
        }
    }
}

impl std::error::Error for SigmfError {}

#[derive(Deserialize)]
struct SigmfMetaFile {
    global: SigmfGlobal,
}

#[derive(Deserialize)]
struct SigmfGlobal {
    #[serde(rename = "core:datatype")]
    datatype: String,
    #[serde(rename = "core:sample_rate")]
    sample_rate: f64,
}

#[derive(Debug, Clone, Copy)]
enum SigmfDatatype {
    Ci8,
    Cu8,
    Ci16Le,
    Ci32Le,
    Cf32Le,
}

impl SigmfDatatype {
    fn parse(s: &str) -> Result<SigmfDatatype, SigmfError> {
        match s {
            "ci8" => Ok(SigmfDatatype::Ci8),
            "cu8" => Ok(SigmfDatatype::Cu8),
            "ci16_le" => Ok(SigmfDatatype::Ci16Le),
            "ci32_le" => Ok(SigmfDatatype::Ci32Le),
            "cf32_le" => Ok(SigmfDatatype::Cf32Le),
            _ => Err(SigmfError::UnsupportedDatatype(s.to_string())),
        }
    }
}

pub struct SigmfStreamer {
    sample_file: BufReader<File>,
    datatype: SigmfDatatype,
    sample_rate: f32,
}

impl SigmfStreamer {
    /// Open a SigMF recording. Pass the path to the `.sigmf-meta` file;
    /// the `.sigmf-data` file is discovered automatically.
    pub fn new(meta_path: &str) -> Result<SigmfStreamer, SigmfError> {
        let meta_path = Path::new(meta_path);

        let meta_file = File::open(meta_path)
            .map_err(|e| SigmfError::BadFile(format!("{}: {}", meta_path.display(), e)))?;
        let meta: SigmfMetaFile = serde_json::from_reader(BufReader::new(meta_file))
            .map_err(|e| SigmfError::BadFile(format!("invalid metadata: {}", e)))?;

        let datatype = SigmfDatatype::parse(&meta.global.datatype)?;
        let sample_rate = meta.global.sample_rate as f32;

        let data_path = meta_path.with_extension("sigmf-data");
        let data_file = File::open(&data_path)
            .map_err(|e| SigmfError::BadFile(format!("{}: {}", data_path.display(), e)))?;

        Ok(SigmfStreamer {
            sample_file: BufReader::new(data_file),
            datatype,
            sample_rate,
        })
    }

    pub fn sample_rate(&self) -> f32 {
        self.sample_rate
    }

    fn read_ci8(&mut self) -> Option<Complex32> {
        let mut buf = [0u8; 2];
        self.sample_file.read_exact(&mut buf).ok()?;
        Some(c32(buf[0] as i8 as f32 / 128.0, buf[1] as i8 as f32 / 128.0))
    }

    fn read_cu8(&mut self) -> Option<Complex32> {
        let mut buf = [0u8; 2];
        self.sample_file.read_exact(&mut buf).ok()?;
        Some(c32(
            (buf[0] as f32 - 128.0) / 128.0,
            (buf[1] as f32 - 128.0) / 128.0,
        ))
    }

    fn read_ci16_le(&mut self) -> Option<Complex32> {
        let mut buf = [0u8; 4];
        self.sample_file.read_exact(&mut buf).ok()?;
        let re = i16::from_le_bytes([buf[0], buf[1]]);
        let im = i16::from_le_bytes([buf[2], buf[3]]);
        Some(c32(re as f32 / 32768.0, im as f32 / 32768.0))
    }

    fn read_ci32_le(&mut self) -> Option<Complex32> {
        let mut buf = [0u8; 8];
        self.sample_file.read_exact(&mut buf).ok()?;
        let re = i32::from_le_bytes([buf[0], buf[1], buf[2], buf[3]]);
        let im = i32::from_le_bytes([buf[4], buf[5], buf[6], buf[7]]);
        Some(c32(re as f32 / 2147483648.0, im as f32 / 2147483648.0))
    }

    fn read_cf32_le(&mut self) -> Option<Complex32> {
        let mut buf = [0u8; 8];
        self.sample_file.read_exact(&mut buf).ok()?;
        let re = f32::from_le_bytes([buf[0], buf[1], buf[2], buf[3]]);
        let im = f32::from_le_bytes([buf[4], buf[5], buf[6], buf[7]]);
        Some(c32(re, im))
    }
}

impl Iterator for SigmfStreamer {
    type Item = Complex32;

    fn next(&mut self) -> Option<Self::Item> {
        match self.datatype {
            SigmfDatatype::Ci8 => self.read_ci8(),
            SigmfDatatype::Cu8 => self.read_cu8(),
            SigmfDatatype::Ci16Le => self.read_ci16_le(),
            SigmfDatatype::Ci32Le => self.read_ci32_le(),
            SigmfDatatype::Cf32Le => self.read_cf32_le(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn try_stream_file() -> Result<(), SigmfError> {
        let file = SigmfStreamer::new("./res/fm_radio_20250920_6msps.sigmf-meta")?;
        assert_eq!(file.sample_rate(), 6000000.0);
        let sum: Complex32 = file.take(100).sum();
        assert_ne!(sum, c32(0.0, 0.0));
        Ok(())
    }
}