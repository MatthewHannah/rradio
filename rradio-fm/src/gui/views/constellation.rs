use eframe::egui;

use crate::constellation_tap::ConstellationSnapshot;

/// RDS constellation scatter plot (I vs Q).
pub fn show(ui: &mut egui::Ui, snapshot: Option<&ConstellationSnapshot>) {
    ui.vertical(|ui| {
        ui.label("RDS Constellation");

        let Some(snap) = snapshot else {
            ui.label("Waiting for RDS chips...");
            return;
        };

        let points: egui_plot::PlotPoints = snap
            .chips
            .iter()
            .map(|c| [c.re as f64, c.im as f64])
            .collect();

        let scatter = egui_plot::Points::new("RDS Chips", points)
            .radius(2.0)
            .color(egui::Color32::from_rgb(80, 180, 255));

        egui_plot::Plot::new("rds_constellation")
            .x_axis_label("In-Phase (I)")
            .y_axis_label("Quadrature (Q)")
            .data_aspect(1.0)
            .auto_bounds(egui::Vec2b::new(false, false))
            .allow_boxed_zoom(true)
            .allow_drag(false)
            .allow_zoom(false)
            .allow_scroll(false)
            .show(ui, |plot_ui| {
                plot_ui.set_plot_bounds(egui_plot::PlotBounds::from_min_max(
                    [-3.0, -3.0],
                    [3.0, 3.0],
                ));
                plot_ui.points(scatter);
            });
    });
}
