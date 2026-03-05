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
    #[must_use]
    pub const fn altitude_range_m(&self) -> (f64, f64) {
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
    #[must_use]
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
#[must_use]
pub fn standard_atmosphere(altitude_m: f64) -> AtmosphericState {
    // Pre-computed reciprocal: 1/287.05 ≈ 3.48368e-3
    const RCP_R_DRY: f64 = 1.0 / 287.05;

    let h = altitude_m.max(0.0);

    let (temperature_c, pressure_hpa) = if h < 11_000.0 {
        // Troposphere: linear temperature lapse, power-law pressure
        let t = 6.5e-3f64.mul_add(-h, 15.0);
        let p = 1013.25 * 2.2558e-5f64.mul_add(-h, 1.0).powf(5.2559);
        (t, p)
    } else if h < 20_000.0 {
        // Lower stratosphere: isothermal at -56.5 °C
        let t = -56.5_f64;
        let p = 226.32 * ((-1.5769e-4) * (h - 11_000.0)).exp();
        (t, p)
    } else if h < 50_000.0 {
        // Upper stratosphere: slight warming, simplified pressure
        let t = 1.0e-3f64.mul_add(h - 20_000.0, -56.5);
        let p_20k = 226.32 * ((-1.5769e-4) * 9_000.0_f64).exp();
        let p = p_20k * ((-1.0e-4) * (h - 20_000.0)).exp();
        (t, p)
    } else if h < 80_000.0 {
        // Mesosphere: temperature drops again
        let t = 2.8e-3f64.mul_add(-(h - 50_000.0), -2.5);
        let p_50k = {
            let p_20k = 226.32 * ((-1.5769e-4) * 9_000.0_f64).exp();
            p_20k * ((-1.0e-4) * 30_000.0_f64).exp()
        };
        let p = p_50k * ((-5.0e-5) * (h - 50_000.0)).exp();
        (t, p)
    } else {
        // Thermosphere: very hot, very low pressure
        let t = 3.0e-3f64.mul_add(h - 80_000.0, 200.0);
        let p = 0.001 * ((-2.0e-5) * (h - 80_000.0)).exp();
        (t, p)
    };

    // Ideal gas law: ρ = P [Pa] / (R_dry [J/(kg·K)] * T [K])
    // P [hPa] → [Pa]: multiply by 100
    let density_kg_m3 = pressure_hpa * 100.0 * RCP_R_DRY / (temperature_c + 273.15);

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
#[must_use]
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
    let seasonal_amplitude = 10.0f64.mul_add(lat_rad.sin().abs(), 5.0);
    let day_angle = 2.0 * std::f64::consts::PI * (f64::from(day_of_year) - 172.0) / 365.0;
    // Flip sign for southern hemisphere
    let seasonal_correction = seasonal_amplitude
        * day_angle.cos()
        * lat.signum().max(1.0)
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
#[must_use]
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
#[allow(clippy::float_cmp)]
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
        assert_eq!(
            AtmosphericLayer::from_altitude(0.0),
            AtmosphericLayer::Troposphere
        );
        assert_eq!(
            AtmosphericLayer::from_altitude(5_000.0),
            AtmosphericLayer::Troposphere
        );
        assert_eq!(
            AtmosphericLayer::from_altitude(20_000.0),
            AtmosphericLayer::Stratosphere
        );
        assert_eq!(
            AtmosphericLayer::from_altitude(65_000.0),
            AtmosphericLayer::Mesosphere
        );
        assert_eq!(
            AtmosphericLayer::from_altitude(100_000.0),
            AtmosphericLayer::Thermosphere
        );
    }

    #[test]
    fn test_temperature_field_varies_with_latitude() {
        // Equator should be warmer than poles at the same altitude and day
        let equator = temperature_field(0.0, 0.0, 0.0, 180);
        let pole = temperature_field(90.0, 0.0, 0.0, 180);
        assert!(
            equator > pole,
            "Equator ({equator:.1}°C) should be warmer than pole ({pole:.1}°C)"
        );
    }

    #[test]
    fn test_density_is_positive() {
        for alt in [
            0.0, 1_000.0, 5_000.0, 11_000.0, 25_000.0, 60_000.0, 100_000.0,
        ] {
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
            "Density at 10 km ({high:.4}) should be less than at sea level ({low:.4})"
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

    // ── Additional tests for improved coverage ──────────────────────

    #[test]
    fn test_layer_from_negative_altitude_clamps_to_troposphere() {
        // Negative altitudes should be clamped to Troposphere
        assert_eq!(
            AtmosphericLayer::from_altitude(-500.0),
            AtmosphericLayer::Troposphere
        );
        assert_eq!(
            AtmosphericLayer::from_altitude(-1.0),
            AtmosphericLayer::Troposphere
        );
        assert_eq!(
            AtmosphericLayer::from_altitude(f64::NEG_INFINITY),
            AtmosphericLayer::Troposphere
        );
    }

    #[test]
    fn test_layer_from_very_high_altitude_is_thermosphere() {
        // Extremely high altitudes remain Thermosphere
        assert_eq!(
            AtmosphericLayer::from_altitude(700_000.0),
            AtmosphericLayer::Thermosphere
        );
        assert_eq!(
            AtmosphericLayer::from_altitude(1_000_000.0),
            AtmosphericLayer::Thermosphere
        );
    }

    #[test]
    fn test_layer_boundary_exact_values() {
        // Test exact boundary values
        assert_eq!(
            AtmosphericLayer::from_altitude(11_999.9),
            AtmosphericLayer::Troposphere
        );
        assert_eq!(
            AtmosphericLayer::from_altitude(12_000.0),
            AtmosphericLayer::Stratosphere
        );
        assert_eq!(
            AtmosphericLayer::from_altitude(49_999.9),
            AtmosphericLayer::Stratosphere
        );
        assert_eq!(
            AtmosphericLayer::from_altitude(50_000.0),
            AtmosphericLayer::Mesosphere
        );
        assert_eq!(
            AtmosphericLayer::from_altitude(79_999.9),
            AtmosphericLayer::Mesosphere
        );
        assert_eq!(
            AtmosphericLayer::from_altitude(80_000.0),
            AtmosphericLayer::Thermosphere
        );
    }

    #[test]
    fn test_all_layer_altitude_ranges() {
        // Verify all four layers have correct ranges
        let strat = AtmosphericLayer::Stratosphere.altitude_range_m();
        assert_eq!(strat, (12_000.0, 50_000.0));

        let meso = AtmosphericLayer::Mesosphere.altitude_range_m();
        assert_eq!(meso, (50_000.0, 80_000.0));
    }

    #[test]
    fn test_isa_negative_altitude_clamps_to_sea_level() {
        // Negative altitude should clamp to 0 (sea level values)
        let neg = standard_atmosphere(-100.0);
        let zero = standard_atmosphere(0.0);
        assert!((neg.temperature_c - zero.temperature_c).abs() < 1e-10);
        assert!((neg.pressure_hpa - zero.pressure_hpa).abs() < 1e-10);
    }

    #[test]
    fn test_isa_lower_stratosphere_isothermal() {
        // Between 11 000 and 20 000 m, temperature should be ~-56.5 C (isothermal)
        let s15 = standard_atmosphere(15_000.0);
        assert!(
            (s15.temperature_c - (-56.5)).abs() < 1.0,
            "Lower stratosphere should be near -56.5 C, got {:.2}",
            s15.temperature_c
        );
    }

    #[test]
    fn test_isa_upper_stratosphere_warms() {
        // Above 20 km, temperature should start warming
        let s30 = standard_atmosphere(30_000.0);
        let s20 = standard_atmosphere(20_000.0);
        assert!(
            s30.temperature_c > s20.temperature_c,
            "Upper stratosphere should be warmer than lower: {:.2} vs {:.2}",
            s30.temperature_c,
            s20.temperature_c
        );
    }

    #[test]
    fn test_isa_mesosphere_cools() {
        // Mesosphere should be colder than upper stratosphere
        let s70 = standard_atmosphere(70_000.0);
        assert!(
            s70.temperature_c < 0.0,
            "Mesosphere should be very cold, got {:.2}",
            s70.temperature_c
        );
    }

    #[test]
    fn test_isa_thermosphere_very_hot() {
        // Thermosphere temperatures increase dramatically
        let s100 = standard_atmosphere(100_000.0);
        assert!(
            s100.temperature_c > 200.0,
            "Thermosphere should be very hot, got {:.2}",
            s100.temperature_c
        );
    }

    #[test]
    fn test_pressure_monotonically_decreases_with_altitude() {
        let altitudes = [
            0.0, 5_000.0, 11_000.0, 15_000.0, 25_000.0, 50_000.0, 80_000.0, 100_000.0,
        ];
        for pair in altitudes.windows(2) {
            let p_low = standard_atmosphere(pair[0]).pressure_hpa;
            let p_high = standard_atmosphere(pair[1]).pressure_hpa;
            assert!(
                p_high < p_low,
                "Pressure at {:.0} m ({:.4} hPa) should be less than at {:.0} m ({:.4} hPa)",
                pair[1],
                p_high,
                pair[0],
                p_low
            );
        }
    }

    #[test]
    fn test_isa_humidity_and_wind_are_zero() {
        // ISA standard atmosphere has zero humidity and wind
        for alt in [0.0, 5_000.0, 20_000.0, 80_000.0] {
            let state = standard_atmosphere(alt);
            assert_eq!(state.humidity_pct, 0.0);
            assert_eq!(state.wind_velocity_ms, [0.0; 3]);
        }
    }

    #[test]
    fn test_temperature_field_seasonal_variation() {
        // In the Northern hemisphere, summer (day ~172) should be warmer than winter (day ~355)
        let summer_t = temperature_field(45.0, 0.0, 0.0, 172);
        let winter_t = temperature_field(45.0, 0.0, 0.0, 355);
        assert!(
            summer_t > winter_t,
            "NH summer ({summer_t:.2}) should be warmer than winter ({winter_t:.2})"
        );
    }

    #[test]
    fn test_temperature_field_southern_hemisphere_reversed_seasons() {
        // Southern hemisphere: opposite seasonal pattern
        let sh_summer = temperature_field(-45.0, 0.0, 0.0, 355); // southern summer
        let sh_winter = temperature_field(-45.0, 0.0, 0.0, 172); // southern winter
        assert!(
            sh_summer > sh_winter,
            "SH summer ({sh_summer:.2}) should be warmer than SH winter ({sh_winter:.2})"
        );
    }

    #[test]
    fn test_wind_trade_winds_near_equator() {
        // Trade winds between 0-30 latitude should blow westward (negative u)
        let [u, _, _] = wind_field(15.0, 0.0, 0.0);
        assert!(
            u < 0.0,
            "Trade winds at 15N should be easterly (negative u), got {u:.4}"
        );
    }

    #[test]
    fn test_wind_westerlies_midlatitudes() {
        // Westerlies between 30-60 latitude should blow eastward (positive u)
        let [u, _, _] = wind_field(45.0, 0.0, 0.0);
        assert!(
            u > 0.0,
            "Westerlies at 45N should be westerly (positive u), got {u:.4}"
        );
    }

    #[test]
    fn test_wind_polar_easterlies() {
        // Polar easterlies above 60 degrees
        let [u, _, _] = wind_field(75.0, 0.0, 0.0);
        assert!(
            u < 0.0,
            "Polar easterlies at 75N should be easterly (negative u), got {u:.4}"
        );
    }

    #[test]
    fn test_wind_increases_with_altitude() {
        let [u_low, _, _] = wind_field(45.0, 0.0, 0.0);
        let [u_high, _, _] = wind_field(45.0, 0.0, 10_000.0);
        assert!(
            u_high.abs() > u_low.abs(),
            "Wind at 10 km ({:.4}) should be stronger than at surface ({:.4})",
            u_high.abs(),
            u_low.abs()
        );
    }

    #[test]
    fn test_wind_vertical_component_is_zero() {
        // The simplified model sets vertical wind to zero everywhere
        for lat in [-60.0, -30.0, 0.0, 30.0, 60.0] {
            let [_, _, w] = wind_field(lat, 0.0, 5_000.0);
            assert_eq!(w, 0.0, "Vertical wind should be zero, got {w:.6}");
        }
    }
}
