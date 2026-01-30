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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quality_to_color_extremes() {
        let scheme = ColorScheme::RedGreen;

        // Low quality (0) should be red
        let low = scheme.quality_to_color(0);
        assert_eq!(low, Color::Rgb(255, 0, 0));

        // High quality (40) should be green
        let high = scheme.quality_to_color(40);
        assert_eq!(high, Color::Rgb(0, 255, 0));
    }

    #[test]
    fn test_no_color_scheme() {
        let scheme = ColorScheme::NoColor;

        // NoColor should always return light gray regardless of quality
        assert_eq!(scheme.quality_to_color(0), Color::Rgb(200, 200, 200));
        assert_eq!(scheme.quality_to_color(20), Color::Rgb(200, 200, 200));
        assert_eq!(scheme.quality_to_color(40), Color::Rgb(200, 200, 200));
    }

    #[test]
    fn test_color_scheme_cycling() {
        // Test that cycling through schemes works correctly
        let red_green = ColorScheme::RedGreen;
        let magenta_green = red_green.next();
        let no_color = magenta_green.next();
        let back_to_red = no_color.next();

        // Verify the cycle by checking color output
        assert_eq!(red_green.quality_to_color(0), Color::Rgb(255, 0, 0));
        assert_eq!(no_color.quality_to_color(0), Color::Rgb(200, 200, 200));
        assert_eq!(back_to_red.quality_to_color(0), Color::Rgb(255, 0, 0));
    }
}
