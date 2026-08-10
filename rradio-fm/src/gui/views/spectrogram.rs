use std::collections::VecDeque;

use eframe::egui;
use eframe::egui::Color32;

const TURBO: [Color32; 10] = [
    Color32::from_rgb(48, 18, 59),
    Color32::from_rgb(35, 106, 141),
    Color32::from_rgb(30, 160, 140),
    Color32::from_rgb(88, 200, 98),
    Color32::from_rgb(164, 223, 39),
    Color32::from_rgb(228, 223, 14),
    Color32::from_rgb(250, 187, 13),
    Color32::from_rgb(246, 135, 8),
    Color32::from_rgb(213, 68, 2),
    Color32::from_rgb(122, 4, 2),
];

/// IQ spectrogram waterfall view.
/// Newest spectrum at top, scrolling downward as new data arrives.
pub fn show(ui: &mut egui::Ui, spectrogram: &VecDeque<Vec<f32>>, sample_rate: f32, is_complex: bool) {
    ui.vertical(|ui| {
        ui.label("IQ Spectrogram");

        if spectrogram.is_empty() {
            ui.label("Waiting for IQ data...");
            return;
        }

        let cols = spectrogram[0].len();

        let values: Vec<f64> = spectrogram
            .iter()
            .flat_map(|row| row.iter().map(|&v| v as f64))
            .collect();

        // Compute frequency axis positioning
        let half_bw_khz = (sample_rate / 2000.0) as f64;
        let (x_origin_khz, bin_width_khz) = if is_complex {
            (-half_bw_khz, (sample_rate / cols as f32) as f64 / 1000.0)
        } else {
            (0.0, (sample_rate / (cols as f32 * 2.0)) as f64 / 1000.0)
        };

        let heatmap = egui_plot::Heatmap::new(values, cols)
            .palette(&TURBO)
            .show_labels(false)
            .at(egui_plot::PlotPoint { x: x_origin_khz, y: 0.0 })
            .tile_size(bin_width_khz as f32, 1.0)
            .name(format!(
                "IQ Spectrogram (±{:.0} kHz)",
                half_bw_khz
            ));

        // Place ticks at round kHz intervals
        let tick_step_khz = nice_step(half_bw_khz * 2.0, 8);

        egui_plot::Plot::new("iq_spectrogram")
            .x_axis_label("Frequency (kHz)")
            .y_axis_label("Time (frames)")
            .show_axes([true, false])
            .allow_zoom(true)
            .allow_drag(true)
            .x_grid_spacer(move |input| {
                let (lo, hi) = input.bounds;
                let mut marks = Vec::new();
                let start = (lo / tick_step_khz).floor() as i64;
                let end = (hi / tick_step_khz).ceil() as i64;
                for i in start..=end {
                    let value = i as f64 * tick_step_khz;
                    marks.push(egui_plot::GridMark {
                        value,
                        step_size: tick_step_khz,
                    });
                }
                marks
            })
            .x_axis_formatter(|mark, _range| {
                format!("{:.0}", mark.value)
            })
            .show(ui, |plot_ui| {
                plot_ui.set_plot_bounds_x(-300.0..=300.0);
                plot_ui.heatmap(heatmap);
            });
    });
}

/// Pick a "nice" step size so roughly `target_ticks` ticks fit across `range`.
fn nice_step(range: f64, target_ticks: usize) -> f64 {
    let raw = range / target_ticks as f64;
    let magnitude = 10.0_f64.powf(raw.log10().floor());
    let residual = raw / magnitude;
    let nice = if residual <= 1.5 {
        1.0
    } else if residual <= 3.5 {
        2.0
    } else if residual <= 7.5 {
        5.0
    } else {
        10.0
    };
    nice * magnitude
}
