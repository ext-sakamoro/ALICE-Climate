//! Atmospheric model as continuous fields.
//!
//! Implements the International Standard Atmosphere (ISA) and simplified
//! field functions for temperature, pressure, humidity, wind, and density
//! at arbitrary altitude, latitude, longitude, and time.

/// Vertical layers of Earth's atmosphere.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtmosphericLayer {
    /// 0 – 12 000 m
    Troposphere,
    /// 12 000 – 50 000 m
    Stratosphere,
    /// 50 000 – 80 000 m
    Mesosphere,
    /// 80 000 – 700 000 m
    Thermosphere,
}

impl AtmosphericLayer {
    /// Altitude range (min, max) in metres for this layer.
    pub fn altitude_range_m(&self) -> (f64, f64) {
        match self {
            Self::Troposphere => (0.0, 12_000.0),
            Self::Stratosphere => (12_000.0, 50_000.0),
            Self::Mesosphere => (50_000.0, 80_000.0),
            Self::Thermosphere => (80_000.0, 700_000.0),
        }
    }

    /// Determine the atmospheric layer for a given altitude (metres).
    ///
    /// Values below 0 are clamped to `Troposphere`. Values above 700 000 m
    /// are also returned as `Thermosphere`.
    pub fn from_altitude(alt_m: f64) -> Self {
        let h = alt_m.max(0.0);
        if h < 12_000.0 {
            Self::Troposphere
        } else if h < 50_000.0 {
            Self::Stratosphere
        } else if h < 80_000.0 {
            Self::Mesosphere
        } else {
            Self::Thermosphere
        }
    }
}

/// Instantaneous state of the atmosphere at a single point.
#[derive(Debug, Clone, Copy)]
pub struct AtmosphericState {
    /// Temperature in degrees Celsius.
    pub temperature_c: f64,
    /// Pressure in hectopascals.
    pub pressure_hpa: f64,
    /// Relative humidity as a percentage (0–100).
    pub humidity_pct: f64,
    /// Wind velocity vector (u, v, w) in m/s.
    /// u = eastward component, v = northward component, w = vertical.
    pub wind_velocity_ms: [f64; 3],
    /// Air density in kg/m³.
    pub density_kg_m3: f64,
}

/// Evaluate the International Standard Atmosphere (ISA) at a given altitude.
///
/// Wind and humidity are zero in the standard atmosphere.
///
/// # Arguments
///
/// * `altitude_m` — Geometric altitude in metres above sea level.
pub fn standard_atmosphere(altitude_m: f64) -> AtmosphericState {
    let h = altitude_m.max(0.0);

    let (temperature_c, pressure_hpa) = if h < 11_000.0 {
        // Troposphere: linear temperature lapse, power-law pressure
        let t = 15.0 - 6.5e-3 * h;
        let p = 1013.25 * (1.0 - 2.2558e-5 * h).powf(5.2559);
        (t, p)
    } else if h < 20_000.0 {
        // Lower stratosphere: isothermal at -56.5 °C
        let t = -56.5_f64;
        let p = 226.32 * ((-1.5769e-4) * (h - 11_000.0)).exp();
        (t, p)
    } else if h < 50_000.0 {
        // Upper stratosphere: slight warming, simplified pressure
        let t = -56.5 + 1.0e-3 * (h - 20_000.0);
        let p_20k = 226.32 * ((-1.5769e-4) * 9_000.0_f64).exp();
        let p = p_20k * ((-1.0e-4) * (h - 20_000.0)).exp();
        (t, p)
    } else if h < 80_000.0 {
        // Mesosphere: temperature drops again
        let t = -2.5 - 2.8e-3 * (h - 50_000.0);
        let p_50k = {
            let p_20k = 226.32 * ((-1.5769e-4) * 9_000.0_f64).exp();
            p_20k * ((-1.0e-4) * 30_000.0_f64).exp()
        };
        let p = p_50k * ((-5.0e-5) * (h - 50_000.0)).exp();
        (t, p)
    } else {
        // Thermosphere: very hot, very low pressure
        let t = 200.0 + 3.0e-3 * (h - 80_000.0);
        let p = 0.001 * ((-2.0e-5) * (h - 80_000.0)).exp();
        (t, p)
    };

    // Ideal gas law: ρ = P [Pa] / (R_dry [J/(kg·K)] * T [K])
    // P [hPa] → [Pa]: multiply by 100
    // Pre-computed reciprocal: 1/287.05 ≈ 3.48368e-3
    const RCP_R_DRY: f64 = 1.0 / 287.05;
    let density_kg_m3 =
        pressure_hpa * 100.0 * RCP_R_DRY / (temperature_c + 273.15);

    AtmosphericState {
        temperature_c,
        pressure_hpa,
        humidity_pct: 0.0,
        wind_velocity_ms: [0.0; 3],
        density_kg_m3,
    }
}

/// Compute temperature (°C) at a geographic location and altitude,
/// incorporating latitude-dependent insolation and seasonal variation.
///
/// # Arguments
///
/// * `lat` — Latitude in decimal degrees (−90 to +90).
/// * `lon` — Longitude in decimal degrees (not used in this simplified model).
/// * `alt_m` — Altitude above sea level in metres.
/// * `day_of_year` — Day of year (1–365).
pub fn temperature_field(lat: f64, _lon: f64, alt_m: f64, day_of_year: u32) -> f64 {
    // ISA base temperature at this altitude
    let isa = standard_atmosphere(alt_m);
    let isa_t = isa.temperature_c;

    // Latitude effect: equator is ~30 °C warmer than poles at surface
    let lat_rad = lat.to_radians();
    let lat_correction = 30.0 * lat_rad.cos();

    // Seasonal variation: simple cosine model
    // Northern hemisphere summer peak around day 172 (June 21)
    // Amplitude: ±15 °C at poles, ±5 °C at equator
    let seasonal_amplitude = 5.0 + 10.0 * lat_rad.sin().abs();
    let day_angle = 2.0 * std::f64::consts::PI * (day_of_year as f64 - 172.0) / 365.0;
    // Flip sign for southern hemisphere
    let seasonal_correction = seasonal_amplitude * day_angle.cos() * lat.signum().max(1.0)
        * if lat >= 0.0 { 1.0 } else { -1.0 };

    isa_t + lat_correction + seasonal_correction
}

/// Compute simplified wind velocity (m/s) at a geographic location and altitude.
///
/// Returns `[u, v, w]` where u = eastward, v = northward, w = vertical.
///
/// Models three major surface wind belts:
/// - Trade winds (~0–30°): blow westward (negative u)
/// - Westerlies (~30–60°): blow eastward (positive u)
/// - Polar easterlies (~60–90°): blow westward (negative u)
///
/// Wind speed increases with altitude (simplified jet stream effect).
///
/// # Arguments
///
/// * `lat` — Latitude in decimal degrees.
/// * `lon` — Longitude in decimal degrees (not used in this model).
/// * `alt_m` — Altitude in metres.
pub fn wind_field(lat: f64, _lon: f64, alt_m: f64) -> [f64; 3] {
    let lat_rad = lat.to_radians();

    // Altitude scaling: wind increases up to ~10 km then levels off
    let alt_factor = (alt_m / 10_000.0).min(3.0) + 1.0;

    // Zonal (east-west) component based on latitude belt
    let abs_lat = lat.abs();
    let u = if abs_lat < 30.0 {
        // Trade winds: easterly (negative u), strongest near 15°
        -8.0 * (lat_rad * 3.0).cos() * alt_factor
    } else if abs_lat < 60.0 {
        // Westerlies: westerly (positive u), strongest near 45°
        let mid = (abs_lat - 45.0).to_radians();
        12.0 * mid.cos() * alt_factor * lat.signum()
    } else {
        // Polar easterlies: easterly (negative u)
        -5.0 * alt_factor * lat.signum()
    };

    // Meridional (north-south) component: small compared to zonal
    let v = 1.5 * (lat_rad * 2.0).sin() * alt_factor;

    // Vertical component: near zero at large scales
    let w = 0.0;

    [u, v, w]
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_isa_sea_level() {
        let state = standard_atmosphere(0.0);
        assert!(
            (state.temperature_c - 15.0).abs() < 0.1,
            "Expected 15°C at sea level, got {:.2}",
            state.temperature_c
        );
        assert!(
            (state.pressure_hpa - 1013.25).abs() < 0.5,
            "Expected 1013.25 hPa at sea level, got {:.2}",
            state.pressure_hpa
        );
    }

    #[test]
    fn test_isa_at_11km() {
        // At the tropopause (~11 000 m) temperature should be near -56.5 °C
        let state = standard_atmosphere(11_000.0);
        assert!(
            (state.temperature_c - (-56.5)).abs() < 1.0,
            "Expected ~-56.5°C at 11 km, got {:.2}",
            state.temperature_c
        );
    }

    #[test]
    fn test_layer_from_altitude() {
        assert_eq!(AtmosphericLayer::from_altitude(0.0), AtmosphericLayer::Troposphere);
        assert_eq!(AtmosphericLayer::from_altitude(5_000.0), AtmosphericLayer::Troposphere);
        assert_eq!(AtmosphericLayer::from_altitude(20_000.0), AtmosphericLayer::Stratosphere);
        assert_eq!(AtmosphericLayer::from_altitude(65_000.0), AtmosphericLayer::Mesosphere);
        assert_eq!(AtmosphericLayer::from_altitude(100_000.0), AtmosphericLayer::Thermosphere);
    }

    #[test]
    fn test_temperature_field_varies_with_latitude() {
        // Equator should be warmer than poles at the same altitude and day
        let equator = temperature_field(0.0, 0.0, 0.0, 180);
        let pole = temperature_field(90.0, 0.0, 0.0, 180);
        assert!(
            equator > pole,
            "Equator ({:.1}°C) should be warmer than pole ({:.1}°C)",
            equator,
            pole
        );
    }

    #[test]
    fn test_density_is_positive() {
        for alt in [0.0, 1_000.0, 5_000.0, 11_000.0, 25_000.0, 60_000.0, 100_000.0] {
            let state = standard_atmosphere(alt);
            assert!(
                state.density_kg_m3 > 0.0,
                "Density must be positive at alt={} m, got {}",
                alt,
                state.density_kg_m3
            );
        }
    }

    #[test]
    fn test_density_decreases_with_altitude() {
        let low = standard_atmosphere(0.0).density_kg_m3;
        let high = standard_atmosphere(10_000.0).density_kg_m3;
        assert!(
            high < low,
            "Density at 10 km ({:.4}) should be less than at sea level ({:.4})",
            high,
            low
        );
    }

    #[test]
    fn test_wind_field_returns_three_components() {
        let w = wind_field(35.0, 139.0, 0.0);
        assert_eq!(w.len(), 3);
    }

    #[test]
    fn test_altitude_range_coverage() {
        let (lo, hi) = AtmosphericLayer::Troposphere.altitude_range_m();
        assert_eq!(lo, 0.0);
        assert_eq!(hi, 12_000.0);

        let (lo2, hi2) = AtmosphericLayer::Thermosphere.altitude_range_m();
        assert_eq!(lo2, 80_000.0);
        assert_eq!(hi2, 700_000.0);
    }
}
