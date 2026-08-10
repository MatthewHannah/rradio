use eframe::egui;

use crate::spectrum_tap::SpectrumSnapshot;

/// MPX magnitude spectrum line chart.
pub fn show(ui: &mut egui::Ui, snapshot: Option<&SpectrumSnapshot>) {
    ui.vertical(|ui| {
        ui.label("MPX Spectrum");

        let Some(snap) = snapshot else {
            ui.label("Waiting for MPX data...");
            return;
        };

        let n_bins = snap.magnitudes_db.len();
        let bin_hz = snap.sample_rate / (n_bins as f32 * 2.0);

        let points: egui_plot::PlotPoints = snap
            .magnitudes_db
            .iter()
            .enumerate()
            .map(|(i, &mag)| [(i as f64) * bin_hz as f64 / 1000.0, mag as f64])
            .collect();

        let line = egui_plot::Line::new("MPX Spectrum", points)
            .color(egui::Color32::from_rgb(255, 160, 0));

        egui_plot::Plot::new("mpx_psd")
            .x_axis_label("Frequency (kHz)")
            .y_axis_label("Power (dB)")
            .include_y(-100.0)
            .include_y(0.0)
            .allow_boxed_zoom(true)
            .allow_drag(true)
            .show(ui, |plot_ui| {
                plot_ui.line(line);

                // Annotation lines for key FM features
                let pilot = egui_plot::VLine::new("Pilot 19 kHz", 19.0)
                    .color(egui::Color32::from_rgba_premultiplied(100, 200, 100, 80));
                plot_ui.vline(pilot);

                let rds = egui_plot::VLine::new("RDS 57 kHz", 57.0)
                    .color(egui::Color32::from_rgba_premultiplied(100, 100, 200, 80));
                plot_ui.vline(rds);
            });
    });
}
