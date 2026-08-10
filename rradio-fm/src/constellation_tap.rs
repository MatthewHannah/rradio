use num_complex::Complex32;

use crate::telemetry::LatestWriter;

/// Snapshot of RDS constellation points for display.
#[derive(Clone)]
pub struct ConstellationSnapshot {
    /// Recent complex-valued chips for I/Q scatter plot.
    pub chips: Vec<Complex32>,
}

/// Iterator adapter that accumulates Complex32 RDS chips
/// and publishes periodic snapshots for a constellation plot.
/// All chips pass through unchanged.
pub struct ConstellationTap<I: Iterator> {
    inner: I,
    ring: Vec<Complex32>,
    capacity: usize,
    pos: usize,
    writer: LatestWriter<ConstellationSnapshot>,
    filled: bool,
}

impl<I> Iterator for ConstellationTap<I>
where
    I: Iterator<Item = Complex32>,
{
    type Item = Complex32;

    fn next(&mut self) -> Option<Self::Item> {
        let chip = self.inner.next()?;

        self.ring[self.pos] = chip;
        self.pos += 1;

        if self.pos >= self.capacity {
            self.pos = 0;
            self.filled = true;
        }

        // Publish a snapshot each time we wrap around,
        // or every `capacity` chips once the ring is full.
        if self.pos == 0 && self.filled {
            self.writer.publish(ConstellationSnapshot {
                chips: self.ring.clone(),
            });
        }

        Some(chip)
    }
}

pub trait ConstellationTappable: Iterator<Item = Complex32> + Sized {
    fn constellation_tap(
        self,
        capacity: usize,
        writer: LatestWriter<ConstellationSnapshot>,
    ) -> ConstellationTap<Self> {
        ConstellationTap {
            inner: self,
            ring: vec![Complex32::new(0.0, 0.0); capacity],
            capacity,
            pos: 0,
            writer,
            filled: false,
        }
    }
}

impl<I: Iterator<Item = Complex32> + Sized> ConstellationTappable for I {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::telemetry::latest_value;

    #[test]
    fn test_constellation_tap_passthrough() {
        let input: Vec<Complex32> = (0..64)
            .map(|i| Complex32::new(i as f32, -(i as f32)))
            .collect();
        let (tx, _rx) = latest_value::<ConstellationSnapshot>();

        let output: Vec<Complex32> = input
            .clone()
            .into_iter()
            .constellation_tap(32, tx)
            .collect();

        assert_eq!(input, output);
    }

    #[test]
    fn test_constellation_tap_publishes() {
        let cap = 16;
        let (tx, rx) = latest_value::<ConstellationSnapshot>();

        // Feed 2 full rings — first fill populates ring, second triggers publish
        let input: Vec<Complex32> = (0..cap * 2)
            .map(|i| Complex32::new(i as f32, 0.0))
            .collect();
        let _: Vec<Complex32> = input.into_iter().constellation_tap(cap, tx).collect();

        let snap = rx.take_latest().expect("should have published");
        assert_eq!(snap.chips.len(), cap);
    }

    #[test]
    fn test_constellation_tap_no_publish_before_full() {
        let cap = 32;
        let (tx, rx) = latest_value::<ConstellationSnapshot>();

        let input: Vec<Complex32> = (0..cap - 1)
            .map(|i| Complex32::new(i as f32, 0.0))
            .collect();
        let _: Vec<Complex32> = input.into_iter().constellation_tap(cap, tx).collect();

        assert!(rx.take_latest().is_none());
    }
}
