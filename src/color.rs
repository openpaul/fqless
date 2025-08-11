use ratatui::style::Color;
pub enum ColorScheme {
    RedGreen,
    MagentaGreen,
    NoColor,
}

impl ColorScheme {
    pub fn quality_to_color(&self, quality_score: u8) -> Color {
        let max_quality = 40;
        let clamped = quality_score.clamp(0, max_quality);
        let normalized = clamped as f32 / max_quality as f32;
        let scale = 255.0;

        match self {
            ColorScheme::RedGreen => {
                let red = ((1.0 - normalized) * scale).round() as u8;
                let green = (normalized * scale).round() as u8;
                Color::Rgb(red, green, 0)
            }

            ColorScheme::MagentaGreen => {
                let red = ((1.0 - normalized) * scale).round() as u8;
                let green = (normalized * scale).round() as u8;
                let blue = ((1.0 - normalized) * scale).round() as u8;
                Color::Rgb(red, green, blue)
            }
            ColorScheme::NoColor => {
                // No color, return a default color (e.g., a light gray)
                Color::Rgb(200, 200, 200) // Light gray
            }
        }
    }
    pub fn next(&self) -> ColorScheme {
        match self {
            ColorScheme::RedGreen => ColorScheme::MagentaGreen,
            ColorScheme::MagentaGreen => ColorScheme::NoColor,
            ColorScheme::NoColor => ColorScheme::RedGreen,
        }
    }
}
