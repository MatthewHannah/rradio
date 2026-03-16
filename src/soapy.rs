use num_complex::Complex32;
use soapysdr::{Device, Direction, RxStream};

static SOAPY_STREAM_SIZE: usize = 32 * 1024;

#[derive(Clone)]
pub struct SoapyConfig {
    pub filter: String,
    pub station: f32,
    pub bw: f32,
    pub fs: f32,
}

pub struct SoapySdr {
    device: Device,
}

pub struct SoapySdrIqStreamer {
    stream: RxStream<Complex32>,
    buf: Vec<Complex32>,
}

impl SoapySdr {
    pub fn connect(config: &SoapyConfig) -> Result<SoapySdr, soapysdr::Error> {
        let device = Device::new(&*config.filter)?;
        device.set_frequency(Direction::Rx, 0, config.station as f64, ())?;
        device.set_sample_rate(Direction::Rx, 0, config.fs as f64)?;
        device.set_bandwidth(Direction::Rx, 0, config.bw as f64)?;
        Ok(SoapySdr { device })
    }

    pub fn start_iq(&self) -> Result<SoapySdrIqStreamer, soapysdr::Error> {
        let mut stream = self.device.rx_stream::<Complex32>(&[0])?;
        stream.activate(None)?;

        let mtu = stream.mtu().unwrap_or(SOAPY_STREAM_SIZE);
        let buf_size = mtu.max(SOAPY_STREAM_SIZE);

        Ok(SoapySdrIqStreamer {
            stream,
            buf: vec![Complex32::new(0.0, 0.0); buf_size],
        })
    }
}

impl SoapySdrIqStreamer {
    pub fn collect_iq(&mut self, data: &mut Vec<Complex32>) -> Result<(), soapysdr::Error> {
        let n = self.stream.read(&mut [&mut self.buf], 1_000_000)?;
        data.clear();
        data.extend_from_slice(&self.buf[..n]);
        Ok(())
    }

    pub fn deactivate(mut self) -> Result<(), soapysdr::Error> {
        self.stream.deactivate(None)
    }
}
