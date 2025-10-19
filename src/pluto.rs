use industrial_io as iio;
use num_complex::Complex32;

pub struct PlutoSdr {
    adc: iio::Device,
    phy: iio::Device, 
    rx_chan_i: iio::Channel,
    rx_chan_q: iio::Channel,
    rx_buf: Option<iio::Buffer>,
}

impl PlutoSdr {
    pub fn new(ctx: iio::Context) -> Option<PlutoSdr> {
        let adc = ctx.find_device("cf-ad9361-lpc")?;
        let phy = ctx.find_device("ad9361-phy")?;

        phy.find_channel("altvoltage0", iio::Direction::Output)?.attr_write_int("frequency", 96100000).ok()?;
        phy.find_channel("voltage0", iio::Direction::Input)?.attr_write_int("sampling_frequency", 6000000).ok()?;
        phy.find_channel("voltage0", iio::Direction::Input)?.attr_write_int("rf_bandwidth", 6000000).ok()?;

        let rx_chan_i = adc.find_channel("voltage0", iio::Direction::Input)?;
        let rx_chan_q = adc.find_channel("voltage1", iio::Direction::Input)?;

        //let rx_buf = adc.create_buffer(128, false).inspect_err(|e| println!("failed to create buffer {:?}", e)).ok()?;

        Some(PlutoSdr {
            adc,
            phy,
            rx_chan_i,
            rx_chan_q,
            rx_buf: None,
        })
    }

    pub fn start(&mut self) -> Result<(), iio::Error>{
        self.rx_chan_i.enable();
        self.rx_chan_q.enable();

        // this has to be 1024*1024 otherwise you drop samples and distort phases
        self.rx_buf = Some(self.adc.create_buffer(1024*1024, false).inspect_err(|e| println!("failed to create buffer {:?}", e))?);

        Ok(())
    }

    pub fn stop(&mut self) {
        self.rx_chan_i.disable();
        self.rx_chan_q.disable();
    }

    pub fn collect_iq(&mut self) -> Option<Vec<Complex32>> {
        self.rx_buf.as_mut()?.refill().inspect_err(|e| println!("refill failed {:?}", e)).ok()?;
        let i = self.rx_chan_i.read::<i16>(self.rx_buf.as_mut()?).ok()?;
        let q = self.rx_chan_q.read::<i16>(self.rx_buf.as_mut()?).ok()?;
        
        Some(i.iter().zip(q.iter()).map(|(&i, &q)| Complex32::new(i as f32 / 4096.0, q as f32 / 4096.0)).collect())
    }

}