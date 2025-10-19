use industrial_io as iio;
use num_complex::Complex32;
use num_traits::Zero;

static PLUTO_SDR_STREAM_SIZE: usize = 1024 * 1024;
static PLUTO_SDR_SCALE: f32 = 4096.0;

pub struct PlutoSdr {
    adc: iio::Device,
    phy: iio::Device,
    streaming: bool
}

pub struct PlutoSdrIqStreamer {
    rx_chan_i: iio::Channel,
    rx_chan_q: iio::Channel,
    rx_buf: iio::Buffer,
}

impl PlutoSdr {
    pub fn new(ctx: iio::Context) -> Option<PlutoSdr> {
        let adc = ctx.find_device("cf-ad9361-lpc")?;
        let phy = ctx.find_device("ad9361-phy")?;

        Some(PlutoSdr { adc, phy, streaming: false })
    }

    pub fn set_center(&self, center: f32) -> Result<(), iio::Error> {
        self.phy
            .find_channel("altvoltage0", iio::Direction::Output)
            .ok_or(iio::Error::General("Missing channel".to_string()))?
            .attr_write_int("frequency", center as i64)
    }

    pub fn set_rf_bandwidth(&self, bw: f32) -> Result<(), iio::Error> {
        self.phy
            .find_channel("voltage0", iio::Direction::Output)
            .ok_or(iio::Error::General("Missing channel".to_string()))?
            .attr_write_int("rf_bandwidth", bw as i64)
    }

    pub fn set_sampling_freq(&self, samp_freq: f32) -> Result<(), iio::Error> {
        self.phy
            .find_channel("voltage0", iio::Direction::Output)
            .ok_or(iio::Error::General("Missing channel".to_string()))?
            .attr_write_int("sampling_frequency", samp_freq as i64)
    }

    pub fn start_iq(&mut self) -> Result<PlutoSdrIqStreamer, iio::Error> {
        if self.streaming {
            return Err(iio::Error::General("Already streaming".to_string()));
        }

        let rx_chan_i = self
            .adc
            .find_channel("voltage0", iio::Direction::Input)
            .ok_or(iio::Error::General("Missing channel".to_string()))?;
        let rx_chan_q = self
            .adc
            .find_channel("voltage1", iio::Direction::Input)
            .ok_or(iio::Error::General("Missing channel".to_string()))?;

        rx_chan_i.enable();
        rx_chan_q.enable();

        // this has to be 1024*1024 otherwise you drop samples and distort phases
        let rx_buf = self
            .adc
            .create_buffer(PLUTO_SDR_STREAM_SIZE, false)
            .inspect_err(|e| println!("failed to create buffer {:?}", e))?;

        self.streaming = true;

        Ok(PlutoSdrIqStreamer {
            rx_chan_i,
            rx_chan_q,
            rx_buf,
        })
    }

    pub fn stop_iq(&mut self, streamer: PlutoSdrIqStreamer) {
        streamer.rx_chan_i.disable();
        streamer.rx_chan_q.disable();

        self.streaming = false;
    }
}

impl PlutoSdrIqStreamer {
    pub fn collect_iq(&mut self, data: &mut Vec<Complex32>) -> Result<(), iio::Error> {
        self.rx_buf
            .refill()
            .inspect_err(|e| println!("refill failed {:?}", e))?;

        // need to enfoce i16 here or else it will interpret bag of bytes incorrectly
        let i_it = self.rx_buf.channel_iter::<i16>(&self.rx_chan_i);
        let q_it = self.rx_buf.channel_iter::<i16>(&self.rx_chan_q);

        data.resize(PLUTO_SDR_STREAM_SIZE, Complex32::zero());

        i_it.zip(q_it)
            .zip(data.iter_mut())
            .for_each(|((&i, &q), o)| {
                *o = Complex32::new((i as f32) / PLUTO_SDR_SCALE, (q as f32) / PLUTO_SDR_SCALE)
            });

        Ok(())
    }
}
