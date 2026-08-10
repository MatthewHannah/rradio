//! Non-blocking, latest-value telemetry for the optional FM visualizer.

use std::sync::{Arc, Mutex};

use crate::constellation_tap::ConstellationSnapshot;
use crate::rds_decoder::RdsSnapshot;
use crate::spectrum_tap::SpectrumSnapshot;

/// Single-producer, single-consumer slot that retains only the newest value.
#[derive(Clone)]
pub struct LatestWriter<T> {
    inner: Arc<Mutex<Option<T>>>,
}

pub struct LatestReader<T> {
    inner: Arc<Mutex<Option<T>>>,
}

pub fn latest_value<T>() -> (LatestWriter<T>, LatestReader<T>) {
    let inner = Arc::new(Mutex::new(None));
    (
        LatestWriter {
            inner: inner.clone(),
        },
        LatestReader { inner },
    )
}

impl<T> LatestWriter<T> {
    pub fn publish(&self, value: T) {
        *self.inner.lock().unwrap() = Some(value);
    }
}

impl<T> LatestReader<T> {
    pub fn take_latest(&self) -> Option<T> {
        self.inner.lock().unwrap().take()
    }
}

/// Publishers consumed by the FM and RDS processing threads.
#[cfg_attr(not(feature = "gui"), allow(dead_code))]
pub struct DisplayPublishers {
    pub signal: SignalDisplayPublishers,
    pub rds: RdsDisplayPublisher,
}

pub struct SignalDisplayPublishers {
    pub iq_spectrum: LatestWriter<SpectrumSnapshot>,
    pub mpx_spectrum: LatestWriter<SpectrumSnapshot>,
}

pub struct RdsDisplayPublisher {
    pub constellation: LatestWriter<ConstellationSnapshot>,
    pub snapshot: LatestWriter<RdsSnapshot>,
}

/// Receivers owned by the GUI thread.
#[cfg_attr(not(feature = "gui"), allow(dead_code))]
pub struct DisplayChannels {
    pub iq_spectrum: LatestReader<SpectrumSnapshot>,
    pub mpx_spectrum: LatestReader<SpectrumSnapshot>,
    pub rds_constellation: LatestReader<ConstellationSnapshot>,
    pub rds_display: LatestReader<RdsSnapshot>,
}

#[cfg_attr(not(feature = "gui"), allow(dead_code))]
pub fn display_channels() -> (DisplayPublishers, DisplayChannels) {
    let (iq_spectrum, iq_spectrum_rx) = latest_value();
    let (mpx_spectrum, mpx_spectrum_rx) = latest_value();
    let (constellation, constellation_rx) = latest_value();
    let (snapshot, rds_display_rx) = latest_value();

    (
        DisplayPublishers {
            signal: SignalDisplayPublishers {
                iq_spectrum,
                mpx_spectrum,
            },
            rds: RdsDisplayPublisher {
                constellation,
                snapshot,
            },
        },
        DisplayChannels {
            iq_spectrum: iq_spectrum_rx,
            mpx_spectrum: mpx_spectrum_rx,
            rds_constellation: constellation_rx,
            rds_display: rds_display_rx,
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn latest_value_overwrites_unread_data() {
        let (tx, rx) = latest_value();
        tx.publish(1);
        tx.publish(2);
        assert_eq!(rx.take_latest(), Some(2));
        assert_eq!(rx.take_latest(), None);
    }
}
