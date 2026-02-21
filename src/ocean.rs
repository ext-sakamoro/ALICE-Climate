//! Ocean model: temperature, salinity, pressure, density, and currents.
//!
//! Provides continuous field functions for ocean state at arbitrary
//! depth, latitude, and longitude. Uses simplified but physically
//! motivated parameterisations.

/// Vertical layers of the ocean.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OceanLayer {
    /// 0 – 200 m: mixed layer, sunlit, well-stirred.
    Surface,
    /// 200 – 1 000 m: rapid temperature drop.
    Thermocline,
    /// 1 000 – 4 000 m: slow, cold, stable.
    DeepOcean,
    /// 4 000 – 11 000 m: trenches, near-freezing.
    Abyssal,
}

impl OceanLayer {
    /// Depth range (min, max) in metres for this layer.
    pub fn depth_range_m(&self) -> (f64, f64) {
        match self {
            Self::Surface => (0.0, 200.0),
            Self::Thermocline => (200.0, 1_000.0),
            Self::DeepOcean => (1_000.0, 4_000.0),
            Self::Abyssal => (4_000.0, 11_000.0),
        }
    }

    /// Determine the ocean layer for a given depth (metres, positive downward).
    pub fn from_depth(depth_m: f64) -> Self {
        let d = depth_m.max(0.0);
        if d < 200.0 {
            Self::Surface
        } else if d < 1_000.0 {
            Self::Thermocline
        } else if d < 4_000.0 {
            Self::DeepOcean
        } else {
            Self::Abyssal
        }
    }
}

/// Instantaneous state of the ocean at a single point.
#[derive(Debug, Clone, Copy)]
pub struct OceanState {
    /// Temperature in degrees Celsius.
    pub temperature_c: f64,
    /// Salinity in Practical Salinity Units (PSU ≈ g/kg).
    pub salinity_psu: f64,
    /// Current velocity vector (u, v, w) in m/s.
    /// u = eastward, v = northward, w = vertical (upwelling/downwelling).
    pub current_velocity_ms: [f64; 3],
    /// Pressure in bar.
    pub pressure_bar: f64,
    /// Seawater density in kg/m³.
    pub density_kg_m3: f64,
}

/// Compute ocean temperature (°C) at a given depth and latitude.
///
/// - Surface (~0 m): 25 °C at equator, ~0 °C at poles.
/// - Thermocline (200–1 000 m): rapid exponential decrease.
/// - Deep ocean (>1 000 m): converges to 2–4 °C globally.
///
/// # Arguments
///
/// * `depth_m` — Depth in metres (positive downward).
/// * `latitude` — Latitude in decimal degrees.
pub fn ocean_temperature(depth_m: f64, latitude: f64) -> f64 {
    let d = depth_m.max(0.0);

    // Surface temperature: 25 °C at equator, 0 °C at poles
    let lat_rad = latitude.to_radians();
    let t_surface = 25.0 * lat_rad.cos().abs();

    if d < 200.0 {
        // Mixed layer: linear interpolation from surface to mixed-layer base
        let t_base = (t_surface - 10.0).max(0.0); // ~10 °C colder at 200 m
        t_surface - (t_surface - t_base) * (d * 5e-3)
    } else if d < 1_000.0 {
        // Thermocline: rapid exponential decay
        let t_200 = (t_surface - 10.0).max(0.0);
        let t_deep = 3.5; // converge to 3.5 °C
        let frac = (d - 200.0) * 1.25e-3; // 0 at 200 m, 1 at 1000 m
        t_200 + (t_deep - t_200) * (1.0 - (-3.0 * frac).exp())
    } else {
        // Deep/Abyssal: slow approach to 2–4 °C depending on location
        let t_deep = 2.0 + 2.0 * lat_rad.cos().abs();
        let t_1km = 3.5;
        let frac = ((d - 1_000.0) * 1e-4).min(1.0);
        t_1km + (t_deep - t_1km) * frac
    }
}

/// Compute ocean pressure (bar) at a given depth.
///
/// Approximation: 1 bar per 10 m of seawater, plus 1 atm (~1.01325 bar)
/// at the surface.
///
/// # Arguments
///
/// * `depth_m` — Depth in metres (positive downward).
pub fn ocean_pressure(depth_m: f64) -> f64 {
    let d = depth_m.max(0.0);
    1.01325 + d * 0.1
}

/// Compute seawater density (kg/m³) using a simplified UNESCO equation.
///
/// ρ ≈ 1025 + 0.8·(S − 35) − 0.2·(T − 10) + 0.05·P
///
/// where S is salinity in PSU, T is temperature in °C, P is pressure in bar.
///
/// # Arguments
///
/// * `temperature_c` — Temperature in degrees Celsius.
/// * `salinity_psu` — Salinity in PSU.
/// * `pressure_bar` — Pressure in bar.
pub fn ocean_density(temperature_c: f64, salinity_psu: f64, pressure_bar: f64) -> f64 {
    1025.0 + 0.8 * (salinity_psu - 35.0) - 0.2 * (temperature_c - 10.0)
        + 0.05 * pressure_bar
}

/// Compute simplified ocean current velocity [u, v, w] in m/s.
///
/// Major systems modelled:
/// - Equatorial current: westward (negative u) near equator (|lat| < 10°)
/// - Gulf Stream / Kuroshio: poleward (positive v) at 30–40°N/S
/// - Antarctic Circumpolar Current: strong eastward at |lat| > 50°
/// - Depth attenuation: exponential with depth
///
/// # Arguments
///
/// * `lat` — Latitude in decimal degrees.
/// * `lon` — Longitude in decimal degrees (not used in this simplified model).
/// * `depth_m` — Depth in metres (positive downward).
pub fn ocean_current(lat: f64, _lon: f64, depth_m: f64) -> [f64; 3] {
    let d = depth_m.max(0.0);

    // Current speed decays exponentially with depth
    let depth_factor = (-d * 2e-3).exp();

    let abs_lat = lat.abs();
    let lat_rad = lat.to_radians();

    let u = if abs_lat < 10.0 {
        // North/South equatorial currents: westward
        -0.5 * (lat_rad * 9.0).cos() * depth_factor
    } else if abs_lat < 45.0 {
        // Subtropical gyres / Western boundary currents: weak eastward drift
        0.3 * depth_factor * lat.signum()
    } else {
        // Antarctic Circumpolar Current / subpolar gyres: strong eastward
        1.0 * depth_factor
    };

    let v = if abs_lat > 25.0 && abs_lat < 45.0 {
        // Gulf Stream / Kuroshio: poleward component
        0.8 * ((abs_lat - 35.0).to_radians().cos()) * depth_factor * lat.signum()
    } else {
        0.1 * lat_rad.sin() * depth_factor
    };

    // Vertical component: small, upwelling at equator, downwelling elsewhere
    let w = if abs_lat < 5.0 {
        1e-4 * depth_factor // slight upwelling at equator
    } else {
        -5e-5 * depth_factor
    };

    [u, v, w]
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_surface_temperature_at_equator() {
        // At equator (lat=0), depth=0 → ~25 °C
        let t = ocean_temperature(0.0, 0.0);
        assert!(
            (t - 25.0).abs() < 1.0,
            "Expected ~25°C at equatorial surface, got {:.2}",
            t
        );
    }

    #[test]
    fn test_surface_temperature_at_pole() {
        // At pole (lat=90), depth=0 → ~0 °C
        let t = ocean_temperature(0.0, 90.0);
        assert!(
            t.abs() < 2.0,
            "Expected ~0°C at polar surface, got {:.2}",
            t
        );
    }

    #[test]
    fn test_pressure_at_1000m() {
        // 1000 m depth → 1.01325 + 1000/10 = 101.01325 bar
        let p = ocean_pressure(1_000.0);
        assert!(
            (p - 101.01325).abs() < 0.1,
            "Expected ~101 bar at 1000 m, got {:.4}",
            p
        );
    }

    #[test]
    fn test_pressure_at_surface() {
        let p = ocean_pressure(0.0);
        assert!((p - 1.01325).abs() < 0.01);
    }

    #[test]
    fn test_density_reasonable_range() {
        // Standard seawater: 1025 ± ~5 kg/m³ at surface
        let rho = ocean_density(15.0, 35.0, 1.01325);
        assert!(
            rho > 1020.0 && rho < 1030.0,
            "Density out of expected range: {:.2}",
            rho
        );
    }

    #[test]
    fn test_density_increases_with_salinity() {
        let rho_low = ocean_density(10.0, 30.0, 1.0);
        let rho_high = ocean_density(10.0, 40.0, 1.0);
        assert!(rho_high > rho_low);
    }

    #[test]
    fn test_density_decreases_with_temperature() {
        let rho_cold = ocean_density(5.0, 35.0, 1.0);
        let rho_warm = ocean_density(25.0, 35.0, 1.0);
        assert!(rho_cold > rho_warm);
    }

    #[test]
    fn test_current_at_equator_has_x_component() {
        let [u, _v, _w] = ocean_current(0.0, 0.0, 0.0);
        // Equatorial current should be westward (negative u)
        assert!(
            u < 0.0,
            "Equatorial current u should be negative (westward), got {:.4}",
            u
        );
    }

    #[test]
    fn test_deep_ocean_temperature() {
        // Below 1000 m, temperature should be 2–4 °C
        let t_deep = ocean_temperature(3_000.0, 30.0);
        assert!(
            t_deep >= 2.0 && t_deep <= 4.5,
            "Deep ocean temperature should be 2–4°C, got {:.2}",
            t_deep
        );
    }

    #[test]
    fn test_current_decays_with_depth() {
        let [u_surf, _, _] = ocean_current(0.0, 0.0, 0.0);
        let [u_deep, _, _] = ocean_current(0.0, 0.0, 5_000.0);
        assert!(
            u_surf.abs() > u_deep.abs(),
            "Current should be weaker at depth"
        );
    }

    #[test]
    fn test_layer_depth_ranges() {
        let (lo, hi) = OceanLayer::Surface.depth_range_m();
        assert_eq!(lo, 0.0);
        assert_eq!(hi, 200.0);

        let (lo2, hi2) = OceanLayer::Abyssal.depth_range_m();
        assert_eq!(lo2, 4_000.0);
        assert_eq!(hi2, 11_000.0);
    }
}
