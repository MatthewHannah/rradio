use std::collections::VecDeque;
use std::sync::atomic;
use std::sync::Arc;

use super::views;
use crate::constellation_tap::ConstellationSnapshot;
use crate::rds_decoder::RdsSnapshot;
use crate::spectrum_tap::SpectrumSnapshot;
use crate::telemetry::DisplayChannels;
use eframe::egui;

/// Cached display state, updated each frame from display slots.
struct DisplayState {
    iq_spectrogram: VecDeque<Vec<f32>>,
    iq_sample_rate: f32,
    iq_is_complex: bool,
    mpx_snapshot: Option<SpectrumSnapshot>,
    constellation_snapshot: Option<ConstellationSnapshot>,
    rds_snapshot: Option<RdsSnapshot>,
    frame_count: u64,
}

const SPECTROGRAM_HISTORY: usize = 200;

pub struct FmDisplayApp {
    channels: DisplayChannels,
    state: DisplayState,
    done_sig: Arc<atomic::AtomicBool>,
}

impl FmDisplayApp {
    pub fn new(
        _cc: &eframe::CreationContext<'_>,
        channels: DisplayChannels,
        done_sig: Arc<atomic::AtomicBool>,
    ) -> Self {
        Self {
            channels,
            done_sig,
            state: DisplayState {
                iq_spectrogram: VecDeque::new(),
                iq_sample_rate: 0.0,
                iq_is_complex: true,
                mpx_snapshot: None,
                constellation_snapshot: None,
                rds_snapshot: None,
                frame_count: 0,
            },
        }
    }

    fn poll_channels(&mut self) {
        // IQ spectrum → push into rolling spectrogram history
        if let Some(snap) = self.channels.iq_spectrum.take_latest() {
            self.state.iq_sample_rate = snap.sample_rate;
            self.state.iq_is_complex = snap.is_complex;
            self.state.iq_spectrogram.push_back(snap.magnitudes_db);
            if self.state.iq_spectrogram.len() > SPECTROGRAM_HISTORY {
                self.state.iq_spectrogram.pop_front();
            }
        }

        if let Some(snap) = self.channels.mpx_spectrum.take_latest() {
            self.state.mpx_snapshot = Some(snap);
        }

        if let Some(snap) = self.channels.rds_constellation.take_latest() {
            self.state.constellation_snapshot = Some(snap);
        }

        if let Some(snap) = self.channels.rds_display.take_latest() {
            self.state.rds_snapshot = Some(snap);
        }

        self.state.frame_count += 1;
    }
}

impl eframe::App for FmDisplayApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        // Exit if the pipeline signalled done (ctrl-c, end of recording, etc.)
        if self.done_sig.load(atomic::Ordering::SeqCst) {
            ui.ctx().send_viewport_cmd(egui::ViewportCommand::Close);
            return;
        }

        self.poll_channels();
        ui.ctx().request_repaint();

        ui.heading("rradio FM Visualizer");
        ui.separator();

        let available = ui.available_size();
        let half_height = available.y / 2.0 - 8.0;

        // Top row: IQ spectrogram | MPX PSD
        ui.allocate_ui(egui::vec2(available.x, half_height), |ui| {
            ui.columns(2, |cols| {
                views::spectrogram::show(
                    &mut cols[0],
                    &self.state.iq_spectrogram,
                    self.state.iq_sample_rate,
                    self.state.iq_is_complex,
                );
                views::psd::show(&mut cols[1], self.state.mpx_snapshot.as_ref());
            });
        });

        ui.separator();

        // Bottom row: RDS constellation | Info panel
        ui.allocate_ui(egui::vec2(available.x, half_height), |ui| {
            ui.columns(2, |cols| {
                views::constellation::show(
                    &mut cols[0],
                    self.state.constellation_snapshot.as_ref(),
                );
                views::info::show(
                    &mut cols[1],
                    self.state.frame_count,
                    self.state.iq_spectrogram.len(),
                    self.state.iq_sample_rate,
                    self.state.mpx_snapshot.as_ref(),
                    self.state.rds_snapshot.as_ref(),
                );
            });
        });
    }
}

/// Launch the eframe GUI on the current (main) thread.
pub fn run_gui(channels: DisplayChannels, done_sig: Arc<atomic::AtomicBool>) -> eframe::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1200.0, 800.0])
            .with_title("rradio FM Visualizer"),
        ..Default::default()
    };

    eframe::run_native(
        "rradio-fm-gui",
        options,
        Box::new(|cc| Ok(Box::new(FmDisplayApp::new(cc, channels, done_sig)))),
    )
}
