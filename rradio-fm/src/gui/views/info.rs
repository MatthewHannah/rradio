use eframe::egui;

use crate::spectrum_tap::SpectrumSnapshot;
use crate::rds_decoder::RdsSnapshot;

/// Info panel showing signal statistics and RDS data.
pub fn show(
    ui: &mut egui::Ui,
    frame_count: u64,
    spectrogram_rows: usize,
    iq_sample_rate: f32,
    mpx_snapshot: Option<&SpectrumSnapshot>,
    rds_snapshot: Option<&RdsSnapshot>,
) {
    ui.vertical(|ui| {
        ui.label("Signal Info");
        ui.separator();

        // RDS data (prominent, at the top)
        if let Some(rds) = rds_snapshot {
            if !rds.ps.is_empty() || rds.pi_code != 0 {
                ui.heading(if rds.ps.is_empty() { "---" } else { &rds.ps });
                if rds.pi_code != 0 {
                    ui.label(format!("PI: {:04X}", rds.pi_code));
                }
                if !rds.rt.is_empty() {
                    ui.label(&rds.rt);
                }
                ui.label(format!(
                    "Groups: {}  BLER: {:.1}%  {}",
                    rds.groups_decoded,
                    rds.bler * 100.0,
                    if rds.synced { "🔒 Synced" } else { "⏳ Searching" }
                ));
                ui.separator();
            }
        }

        ui.label(format!("GUI Frame: {}", frame_count));
        ui.label(format!("Spectrogram Rows: {}", spectrogram_rows));

        if iq_sample_rate > 0.0 {
            ui.label(format!("IQ Sample Rate: {:.0} kHz", iq_sample_rate / 1000.0));
            ui.label(format!("IQ Bandwidth: ±{:.0} kHz", iq_sample_rate / 2000.0));
        }

        if let Some(snap) = mpx_snapshot {
            ui.label(format!("MPX Sample Rate: {:.0} kHz", snap.sample_rate / 1000.0));
            ui.label(format!("MPX FFT Bins: {}", snap.magnitudes_db.len()));

            if let Some((peak_bin, &peak_db)) = snap
                .magnitudes_db
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            {
                let bin_hz = snap.sample_rate / (snap.magnitudes_db.len() as f32 * 2.0);
                let peak_freq = peak_bin as f32 * bin_hz;
                ui.label(format!("MPX Peak: {:.1} kHz ({:.1} dB)", peak_freq / 1000.0, peak_db));
            }
        }
    });
}
