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

/// SigMF recording writer. Writes cf32_le IQ data and a compliant metadata file.
pub struct SigmfWriter {
    data_file: std::io::BufWriter<File>,
    meta_path: std::path::PathBuf,
    sample_rate: f64,
    frequency: f64,
    hw: String,
    samples_written: u64,
    datetime: String,
}

impl SigmfWriter {
    /// Create a new SigMF recording. `base_path` is the path without extension;
    /// `.sigmf-data` and `.sigmf-meta` will be appended.
    pub fn new(base_path: &str, sample_rate: f64, frequency: f64, hw: &str) -> Result<SigmfWriter, SigmfError> {
        let data_path = format!("{}.sigmf-data", base_path);
        let data_file = File::create(&data_path)
            .map_err(|e| SigmfError::BadFile(format!("{}: {}", data_path, e)))?;

        // ISO 8601 UTC timestamp
        let now = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap_or_default()
            .as_secs();
        let datetime = format_utc_timestamp(now);

        Ok(SigmfWriter {
            data_file: std::io::BufWriter::new(data_file),
            meta_path: std::path::PathBuf::from(format!("{}.sigmf-meta", base_path)),
            sample_rate,
            frequency,
            hw: hw.to_string(),
            samples_written: 0,
            datetime,
        })
    }

    /// Write a batch of IQ samples.
    pub fn write_samples(&mut self, samples: &[Complex32]) -> Result<(), SigmfError> {
        use std::io::Write;
        for s in samples {
            self.data_file.write_all(&s.re.to_le_bytes())
                .map_err(|e| SigmfError::BadFile(format!("write error: {}", e)))?;
            self.data_file.write_all(&s.im.to_le_bytes())
                .map_err(|e| SigmfError::BadFile(format!("write error: {}", e)))?;
        }
        self.samples_written += samples.len() as u64;
        Ok(())
    }

    /// Finalize the recording: flush data and write the metadata file.
    pub fn finalize(mut self) -> Result<(), SigmfError> {
        use std::io::Write;
        self.data_file.flush()
            .map_err(|e| SigmfError::BadFile(format!("flush error: {}", e)))?;

        let meta = serde_json::json!({
            "global": {
                "core:datatype": "cf32_le",
                "core:sample_rate": self.sample_rate,
                "core:version": "1.2.0",
                "core:hw": self.hw,
                "core:recorder": "rradio",
                "core:description": format!("FM recording at {:.1} MHz", self.frequency / 1e6),
            },
            "captures": [{
                "core:sample_start": 0,
                "core:frequency": self.frequency,
                "core:datetime": self.datetime,
            }],
            "annotations": []
        });

        let meta_file = File::create(&self.meta_path)
            .map_err(|e| SigmfError::BadFile(format!("{}: {}", self.meta_path.display(), e)))?;
        serde_json::to_writer_pretty(std::io::BufWriter::new(meta_file), &meta)
            .map_err(|e| SigmfError::BadFile(format!("metadata write error: {}", e)))?;

        let duration = self.samples_written as f64 / self.sample_rate;
        eprintln!("SigMF: wrote {} samples ({:.1}s) to {}",
            self.samples_written, duration, self.meta_path.display());
        Ok(())
    }
}

fn format_utc_timestamp(epoch_secs: u64) -> String {
    let secs_per_day = 86400u64;
    let days = epoch_secs / secs_per_day;
    let time_of_day = epoch_secs % secs_per_day;
    let hours = time_of_day / 3600;
    let minutes = (time_of_day % 3600) / 60;
    let seconds = time_of_day % 60;

    // Simple date calculation from days since epoch (1970-01-01)
    let mut y = 1970i64;
    let mut remaining = days as i64;
    loop {
        let days_in_year = if y % 4 == 0 && (y % 100 != 0 || y % 400 == 0) { 366 } else { 365 };
        if remaining < days_in_year { break; }
        remaining -= days_in_year;
        y += 1;
    }
    let leap = y % 4 == 0 && (y % 100 != 0 || y % 400 == 0);
    let month_days = [31, if leap { 29 } else { 28 }, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    let mut m = 0;
    for &md in &month_days {
        if remaining < md { break; }
        remaining -= md;
        m += 1;
    }

    format!("{:04}-{:02}-{:02}T{:02}:{:02}:{:02}Z", y, m + 1, remaining + 1, hours, minutes, seconds)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn try_stream_file() -> Result<(), SigmfError> {
        let file = SigmfStreamer::new("../res/fm_radio_20250920_6msps.sigmf-meta")?;
        assert_eq!(file.sample_rate(), 6000000.0);
        let sum: Complex32 = file.take(100).sum();
        assert_ne!(sum, c32(0.0, 0.0));
        Ok(())
    }
}