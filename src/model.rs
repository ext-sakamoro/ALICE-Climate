//! Unified climate model integrating atmosphere and ocean.
//!
//! Provides the top-level [`evaluate_climate`] function and anomaly detection
//! logic. Queries can be issued at any geographic location and altitude;
//! negative altitudes are interpreted as ocean depth.

use crate::atmosphere::{
    standard_atmosphere, temperature_field, wind_field, AtmosphericLayer, AtmosphericState,
};
use crate::ocean::{ocean_current, ocean_density, ocean_pressure, ocean_temperature, OceanState};

/// A point-in-space-and-time query for climate state.
#[derive(Debug, Clone, Copy)]
pub struct ClimateQuery {
    /// Latitude in decimal degrees (−90 to +90).
    pub latitude: f64,
    /// Longitude in decimal degrees (−180 to +180).
    pub longitude: f64,
    /// Altitude above sea level in metres.
    /// Negative values are treated as ocean depth (|altitude_m|).
    pub altitude_m: f64,
    /// Timestamp in nanoseconds since Unix epoch.
    pub timestamp_ns: u64,
}

/// The climate state returned for a [`ClimateQuery`].
#[derive(Debug, Clone, Copy)]
pub struct ClimateResponse {
    /// Atmospheric state at the query point.
    pub atmosphere: AtmosphericState,
    /// Ocean state at the query point, or `None` if above sea level.
    pub ocean: Option<OceanState>,
    /// The atmospheric layer corresponding to the query altitude.
    pub layer: AtmosphericLayer,
}

/// Categories of climate anomaly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnomalyKind {
    /// Temperature significantly above baseline.
    HeatWave,
    /// Temperature significantly below baseline.
    ColdSnap,
    /// Very high wind speeds.
    Storm,
    /// Prolonged low humidity.
    Drought,
    /// Extreme precipitation or flooding conditions.
    Flood,
}

/// A detected climate anomaly at a geographic location and time.
#[derive(Debug, Clone, Copy)]
pub struct ClimateAnomaly {
    /// Latitude of the anomaly centre.
    pub location_lat: f64,
    /// Longitude of the anomaly centre.
    pub location_lon: f64,
    /// Type of anomaly detected.
    pub kind: AnomalyKind,
    /// Magnitude of the anomaly (e.g., temperature deviation in °C,
    /// or wind speed in m/s for storms).
    pub magnitude: f64,
    /// Timestamp when the anomaly was detected (nanoseconds since Unix epoch).
    pub timestamp_ns: u64,
    /// FNV-1a hash of the anomaly description for deduplication.
    pub content_hash: u64,
}

/// FNV-1a 64-bit hash used to generate `content_hash` for anomalies.
#[inline(always)]
fn fnv1a_f64_pair(a: f64, b: f64, kind: u8) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &byte in a
        .to_bits()
        .to_le_bytes()
        .iter()
        .chain(b.to_bits().to_le_bytes().iter())
        .chain(std::slice::from_ref(&kind))
    {
        h ^= byte as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

/// Evaluate the full climate state at a single query point.
///
/// If `query.altitude_m < 0`, the point is treated as being at ocean depth
/// `|altitude_m|` metres, and an [`OceanState`] is computed.
///
/// The atmospheric state is always computed (useful for sub-sea atmospheric
/// baseline comparisons and surface boundary-layer coupling).
pub fn evaluate_climate(query: &ClimateQuery) -> ClimateResponse {
    let alt = query.altitude_m;
    let lat = query.latitude;
    let lon = query.longitude;

    // Day of year from timestamp (nanoseconds → seconds → days, modulo 365)
    let seconds = query.timestamp_ns / 1_000_000_000;
    let day_of_year = ((seconds / 86_400) % 365 + 1) as u32;

    // Atmospheric state (always computed)
    let isa = standard_atmosphere(alt.max(0.0));
    let temperature_c = temperature_field(lat, lon, alt.max(0.0), day_of_year);
    let wind_velocity_ms = wind_field(lat, lon, alt.max(0.0));
    let layer = AtmosphericLayer::from_altitude(alt.max(0.0));

    let atmosphere = AtmosphericState {
        temperature_c,
        pressure_hpa: isa.pressure_hpa,
        humidity_pct: base_humidity(lat, alt.max(0.0)),
        wind_velocity_ms,
        density_kg_m3: isa.pressure_hpa * (100.0 / 287.05) / (temperature_c + 273.15),
    };

    // Ocean state: only if below sea level (negative altitude = depth)
    let ocean = if alt < 0.0 {
        let depth = (-alt).max(0.0);
        let temperature_c = ocean_temperature(depth, lat);
        let salinity_psu = 35.0; // simplified constant
        let pressure_bar = ocean_pressure(depth);
        let density_kg_m3 = ocean_density(temperature_c, salinity_psu, pressure_bar);
        let current_velocity_ms = ocean_current(lat, lon, depth);

        Some(OceanState {
            temperature_c,
            salinity_psu,
            current_velocity_ms,
            pressure_bar,
            density_kg_m3,
        })
    } else {
        None
    };

    ClimateResponse {
        atmosphere,
        ocean,
        layer,
    }
}

/// Evaluate climate at a batch of query points.
///
/// Results are returned in the same order as `queries`.
pub fn evaluate_climate_batch(queries: &[ClimateQuery]) -> Vec<ClimateResponse> {
    queries.iter().map(evaluate_climate).collect()
}

/// Detect a climate anomaly by comparing current state against a baseline.
///
/// Thresholds applied:
/// - Temperature deviation > +10 °C → [`AnomalyKind::HeatWave`]
/// - Temperature deviation < −10 °C → [`AnomalyKind::ColdSnap`]
/// - Wind speed (magnitude of u,v components) > 30 m/s → [`AnomalyKind::Storm`]
/// - Humidity < 10 % while baseline > 40 % → [`AnomalyKind::Drought`]
/// - Humidity > 90 % while baseline < 60 % → [`AnomalyKind::Flood`]
///
/// Returns `None` if no anomaly threshold is exceeded.
pub fn detect_anomaly(
    current: &AtmosphericState,
    baseline: &AtmosphericState,
    lat: f64,
    lon: f64,
    timestamp_ns: u64,
) -> Option<ClimateAnomaly> {
    let temp_dev = current.temperature_c - baseline.temperature_c;

    let u = current.wind_velocity_ms[0];
    let v = current.wind_velocity_ms[1];
    let wind_speed = (u * u + v * v).sqrt();

    let (kind, magnitude) = if temp_dev > 10.0 {
        (AnomalyKind::HeatWave, temp_dev)
    } else if temp_dev < -10.0 {
        (AnomalyKind::ColdSnap, temp_dev.abs())
    } else if wind_speed > 30.0 {
        (AnomalyKind::Storm, wind_speed)
    } else if current.humidity_pct < 10.0 && baseline.humidity_pct > 40.0 {
        (
            AnomalyKind::Drought,
            baseline.humidity_pct - current.humidity_pct,
        )
    } else if current.humidity_pct > 90.0 && baseline.humidity_pct < 60.0 {
        (
            AnomalyKind::Flood,
            current.humidity_pct - baseline.humidity_pct,
        )
    } else {
        return None;
    };

    let kind_byte = match kind {
        AnomalyKind::HeatWave => 0u8,
        AnomalyKind::ColdSnap => 1,
        AnomalyKind::Storm => 2,
        AnomalyKind::Drought => 3,
        AnomalyKind::Flood => 4,
    };

    Some(ClimateAnomaly {
        location_lat: lat,
        location_lon: lon,
        kind,
        magnitude,
        timestamp_ns,
        content_hash: fnv1a_f64_pair(lat, lon, kind_byte),
    })
}

/// Simple latitude/altitude-based humidity estimate (%).
///
/// Tropics are humid; high altitude and polar regions are dry.
fn base_humidity(lat: f64, alt_m: f64) -> f64 {
    let lat_rad = lat.to_radians();
    // Tropical belt: ~80%, polar: ~40%, altitude dries air
    let base = 40.0 + 40.0 * lat_rad.cos();
    let alt_dry = (alt_m / 5_000.0) * 20.0;
    (base - alt_dry).clamp(0.0, 100.0)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_query(lat: f64, lon: f64, alt: f64) -> ClimateQuery {
        ClimateQuery {
            latitude: lat,
            longitude: lon,
            altitude_m: alt,
            timestamp_ns: 180u64 * 24 * 3600 * 1_000_000_000,
        }
    }

    #[test]
    fn test_evaluate_at_sea_level() {
        let q = make_query(35.6, 139.7, 0.0);
        let r = evaluate_climate(&q);
        // Atmosphere should be populated
        assert!(r.atmosphere.pressure_hpa > 900.0);
        assert!(r.atmosphere.density_kg_m3 > 0.0);
        // No ocean at surface level
        assert!(r.ocean.is_none());
        assert_eq!(r.layer, AtmosphericLayer::Troposphere);
    }

    #[test]
    fn test_evaluate_below_sea_level_gives_ocean() {
        // altitude = -500 m → depth = 500 m (thermocline)
        let q = make_query(0.0, 0.0, -500.0);
        let r = evaluate_climate(&q);
        assert!(
            r.ocean.is_some(),
            "Ocean state should be present below sea level"
        );
        let ocean = r.ocean.unwrap();
        assert!(ocean.temperature_c > 0.0);
        assert!(ocean.pressure_bar > 40.0); // 500m → ~51 bar
        assert!(ocean.density_kg_m3 > 1_000.0);
    }

    #[test]
    fn test_batch_matches_individual() {
        let queries = vec![
            make_query(0.0, 0.0, 0.0),
            make_query(45.0, 90.0, 5_000.0),
            make_query(-30.0, -60.0, 0.0),
        ];
        let batch = evaluate_climate_batch(&queries);
        for (i, q) in queries.iter().enumerate() {
            let single = evaluate_climate(q);
            assert!(
                (batch[i].atmosphere.temperature_c - single.atmosphere.temperature_c).abs() < 1e-10,
                "Batch and single results differ at index {}",
                i
            );
            assert!(
                (batch[i].atmosphere.pressure_hpa - single.atmosphere.pressure_hpa).abs() < 1e-10,
                "Pressure mismatch at index {}",
                i
            );
        }
    }

    #[test]
    fn test_anomaly_detection_extreme_heat() {
        let baseline = AtmosphericState {
            temperature_c: 20.0,
            pressure_hpa: 1013.25,
            humidity_pct: 60.0,
            wind_velocity_ms: [2.0, 1.0, 0.0],
            density_kg_m3: 1.2,
        };
        let current = AtmosphericState {
            temperature_c: 45.0, // +25 °C deviation
            ..baseline
        };
        let anomaly = detect_anomaly(&current, &baseline, 35.0, 135.0, 0);
        assert!(anomaly.is_some());
        let a = anomaly.unwrap();
        assert_eq!(a.kind, AnomalyKind::HeatWave);
        assert!((a.magnitude - 25.0).abs() < 0.01);
    }

    #[test]
    fn test_anomaly_detection_cold_snap() {
        let baseline = AtmosphericState {
            temperature_c: 10.0,
            pressure_hpa: 1013.25,
            humidity_pct: 60.0,
            wind_velocity_ms: [0.0; 3],
            density_kg_m3: 1.2,
        };
        let current = AtmosphericState {
            temperature_c: -15.0, // -25 °C deviation
            ..baseline
        };
        let anomaly = detect_anomaly(&current, &baseline, 50.0, 10.0, 0);
        assert!(anomaly.is_some());
        assert_eq!(anomaly.unwrap().kind, AnomalyKind::ColdSnap);
    }

    #[test]
    fn test_anomaly_detection_storm() {
        let baseline = AtmosphericState {
            temperature_c: 20.0,
            pressure_hpa: 1013.25,
            humidity_pct: 60.0,
            wind_velocity_ms: [5.0, 2.0, 0.0],
            density_kg_m3: 1.2,
        };
        let current = AtmosphericState {
            wind_velocity_ms: [35.0, 15.0, 0.0], // ~38 m/s
            ..baseline
        };
        let anomaly = detect_anomaly(&current, &baseline, 20.0, 150.0, 0);
        assert!(anomaly.is_some());
        assert_eq!(anomaly.unwrap().kind, AnomalyKind::Storm);
    }

    #[test]
    fn test_no_anomaly_for_normal_conditions() {
        let state = AtmosphericState {
            temperature_c: 20.0,
            pressure_hpa: 1013.25,
            humidity_pct: 60.0,
            wind_velocity_ms: [3.0, 1.0, 0.0],
            density_kg_m3: 1.2,
        };
        // Same state for baseline and current → no deviation
        let anomaly = detect_anomaly(&state, &state, 35.0, 139.0, 0);
        assert!(
            anomaly.is_none(),
            "No anomaly expected for identical states"
        );
    }

    #[test]
    fn test_evaluate_high_altitude_stratosphere() {
        let q = make_query(0.0, 0.0, 25_000.0);
        let r = evaluate_climate(&q);
        assert_eq!(r.layer, AtmosphericLayer::Stratosphere);
        assert!(r.ocean.is_none());
        assert!(r.atmosphere.pressure_hpa < 100.0); // low pressure at 25 km
    }

    #[test]
    fn test_content_hash_determinism() {
        let baseline = AtmosphericState {
            temperature_c: 20.0,
            pressure_hpa: 1013.25,
            humidity_pct: 60.0,
            wind_velocity_ms: [0.0; 3],
            density_kg_m3: 1.2,
        };
        let current = AtmosphericState {
            temperature_c: 50.0,
            ..baseline
        };
        let a1 = detect_anomaly(&current, &baseline, 35.0, 139.0, 0).unwrap();
        let a2 = detect_anomaly(&current, &baseline, 35.0, 139.0, 0).unwrap();
        assert_eq!(a1.content_hash, a2.content_hash);
    }

    // ── Additional tests for improved coverage ──────────────────────

    #[test]
    fn test_anomaly_drought_detection() {
        let baseline = AtmosphericState {
            temperature_c: 20.0,
            pressure_hpa: 1013.25,
            humidity_pct: 60.0, // baseline > 40%
            wind_velocity_ms: [0.0; 3],
            density_kg_m3: 1.2,
        };
        let current = AtmosphericState {
            humidity_pct: 5.0, // current < 10%
            ..baseline
        };
        let anomaly = detect_anomaly(&current, &baseline, 30.0, 0.0, 0);
        assert!(anomaly.is_some());
        let a = anomaly.unwrap();
        assert_eq!(a.kind, AnomalyKind::Drought);
        assert!(
            (a.magnitude - 55.0).abs() < 0.01,
            "Drought magnitude should be 55.0, got {:.2}",
            a.magnitude
        );
    }

    #[test]
    fn test_anomaly_flood_detection() {
        let baseline = AtmosphericState {
            temperature_c: 20.0,
            pressure_hpa: 1013.25,
            humidity_pct: 50.0, // baseline < 60%
            wind_velocity_ms: [0.0; 3],
            density_kg_m3: 1.2,
        };
        let current = AtmosphericState {
            humidity_pct: 95.0, // current > 90%
            ..baseline
        };
        let anomaly = detect_anomaly(&current, &baseline, 10.0, 100.0, 0);
        assert!(anomaly.is_some());
        let a = anomaly.unwrap();
        assert_eq!(a.kind, AnomalyKind::Flood);
        assert!((a.magnitude - 45.0).abs() < 0.01);
    }

    #[test]
    fn test_anomaly_priority_heat_over_storm() {
        // When both temperature deviation and wind speed exceed thresholds,
        // heat wave should take priority (checked first)
        let baseline = AtmosphericState {
            temperature_c: 20.0,
            pressure_hpa: 1013.25,
            humidity_pct: 60.0,
            wind_velocity_ms: [0.0; 3],
            density_kg_m3: 1.2,
        };
        let current = AtmosphericState {
            temperature_c: 45.0,                 // +25 deviation
            wind_velocity_ms: [35.0, 15.0, 0.0], // ~38 m/s
            ..baseline
        };
        let anomaly = detect_anomaly(&current, &baseline, 0.0, 0.0, 0);
        assert!(anomaly.is_some());
        assert_eq!(anomaly.unwrap().kind, AnomalyKind::HeatWave);
    }

    #[test]
    fn test_anomaly_cold_snap_magnitude_is_absolute() {
        let baseline = AtmosphericState {
            temperature_c: 5.0,
            pressure_hpa: 1013.25,
            humidity_pct: 60.0,
            wind_velocity_ms: [0.0; 3],
            density_kg_m3: 1.2,
        };
        let current = AtmosphericState {
            temperature_c: -20.0, // -25 deviation
            ..baseline
        };
        let anomaly = detect_anomaly(&current, &baseline, 0.0, 0.0, 0).unwrap();
        assert_eq!(anomaly.kind, AnomalyKind::ColdSnap);
        assert!(
            (anomaly.magnitude - 25.0).abs() < 0.01,
            "ColdSnap magnitude should be abs value"
        );
    }

    #[test]
    fn test_anomaly_below_thresholds_returns_none() {
        let baseline = AtmosphericState {
            temperature_c: 20.0,
            pressure_hpa: 1013.25,
            humidity_pct: 50.0,
            wind_velocity_ms: [5.0, 3.0, 0.0],
            density_kg_m3: 1.2,
        };
        // Temperature deviation of +9 C (below 10 C threshold)
        let current = AtmosphericState {
            temperature_c: 29.0,
            wind_velocity_ms: [10.0, 5.0, 0.0], // ~11 m/s (below 30)
            humidity_pct: 50.0,
            ..baseline
        };
        let anomaly = detect_anomaly(&current, &baseline, 0.0, 0.0, 0);
        assert!(
            anomaly.is_none(),
            "Should not detect anomaly when below all thresholds"
        );
    }

    #[test]
    fn test_content_hash_differs_for_different_locations() {
        let baseline = AtmosphericState {
            temperature_c: 20.0,
            pressure_hpa: 1013.25,
            humidity_pct: 60.0,
            wind_velocity_ms: [0.0; 3],
            density_kg_m3: 1.2,
        };
        let current = AtmosphericState {
            temperature_c: 50.0,
            ..baseline
        };
        let a1 = detect_anomaly(&current, &baseline, 35.0, 139.0, 0).unwrap();
        let a2 = detect_anomaly(&current, &baseline, 40.0, 139.0, 0).unwrap();
        assert_ne!(
            a1.content_hash, a2.content_hash,
            "Different locations should have different hashes"
        );
    }

    #[test]
    fn test_content_hash_nonzero() {
        let baseline = AtmosphericState {
            temperature_c: 20.0,
            pressure_hpa: 1013.25,
            humidity_pct: 60.0,
            wind_velocity_ms: [0.0; 3],
            density_kg_m3: 1.2,
        };
        let current = AtmosphericState {
            temperature_c: 50.0,
            ..baseline
        };
        let a = detect_anomaly(&current, &baseline, 0.0, 0.0, 0).unwrap();
        assert_ne!(a.content_hash, 0, "Content hash should never be zero");
    }

    #[test]
    fn test_evaluate_climate_day_of_year_calculation() {
        // Day 1 of year: timestamp = 0 should map to day 1
        let q = ClimateQuery {
            latitude: 0.0,
            longitude: 0.0,
            altitude_m: 0.0,
            timestamp_ns: 0,
        };
        let r = evaluate_climate(&q);
        // Should not panic and should produce valid values
        assert!(r.atmosphere.temperature_c > -100.0 && r.atmosphere.temperature_c < 100.0);
    }

    #[test]
    fn test_evaluate_climate_extreme_latitude() {
        // North pole
        let q_north = make_query(90.0, 0.0, 0.0);
        let r_north = evaluate_climate(&q_north);
        assert!(r_north.atmosphere.pressure_hpa > 0.0);

        // South pole
        let q_south = make_query(-90.0, 0.0, 0.0);
        let r_south = evaluate_climate(&q_south);
        assert!(r_south.atmosphere.pressure_hpa > 0.0);
    }

    #[test]
    fn test_evaluate_climate_humidity_tropics_vs_poles() {
        // Tropics should be more humid than polar regions at surface
        let q_tropic = make_query(0.0, 0.0, 0.0);
        let q_polar = make_query(80.0, 0.0, 0.0);
        let r_tropic = evaluate_climate(&q_tropic);
        let r_polar = evaluate_climate(&q_polar);
        assert!(
            r_tropic.atmosphere.humidity_pct > r_polar.atmosphere.humidity_pct,
            "Tropical humidity ({:.1}%) should exceed polar ({:.1}%)",
            r_tropic.atmosphere.humidity_pct,
            r_polar.atmosphere.humidity_pct
        );
    }

    #[test]
    fn test_evaluate_climate_humidity_decreases_with_altitude() {
        let q_low = make_query(30.0, 0.0, 0.0);
        let q_high = make_query(30.0, 0.0, 8_000.0);
        let r_low = evaluate_climate(&q_low);
        let r_high = evaluate_climate(&q_high);
        assert!(
            r_low.atmosphere.humidity_pct > r_high.atmosphere.humidity_pct,
            "Humidity at surface ({:.1}%) should exceed at 8km ({:.1}%)",
            r_low.atmosphere.humidity_pct,
            r_high.atmosphere.humidity_pct
        );
    }

    #[test]
    fn test_batch_empty_input() {
        let batch = evaluate_climate_batch(&[]);
        assert!(batch.is_empty(), "Empty input should produce empty output");
    }

    #[test]
    fn test_batch_single_element() {
        let q = make_query(35.0, 139.0, 0.0);
        let batch = evaluate_climate_batch(&[q]);
        assert_eq!(batch.len(), 1);
        let single = evaluate_climate(&q);
        assert!(
            (batch[0].atmosphere.temperature_c - single.atmosphere.temperature_c).abs() < 1e-10
        );
    }

    #[test]
    fn test_evaluate_deep_ocean_query() {
        // Very deep ocean point
        let q = make_query(0.0, 0.0, -5_000.0);
        let r = evaluate_climate(&q);
        assert!(r.ocean.is_some());
        let ocean = r.ocean.unwrap();
        assert!(ocean.temperature_c >= 1.0 && ocean.temperature_c <= 5.0);
        assert!(ocean.pressure_bar > 400.0); // 5000m ~= 501 bar
        assert!(ocean.density_kg_m3 > 1020.0);
        // Atmosphere should still be valid (surface-level baseline)
        assert_eq!(r.layer, AtmosphericLayer::Troposphere);
    }

    #[test]
    fn test_base_humidity_clamp_high_altitude() {
        // At very high altitude, humidity should be clamped to >= 0
        let q = make_query(0.0, 0.0, 20_000.0);
        let r = evaluate_climate(&q);
        assert!(
            r.atmosphere.humidity_pct >= 0.0,
            "Humidity should be >= 0, got {:.2}",
            r.atmosphere.humidity_pct
        );
    }

    #[test]
    fn test_evaluate_mesosphere_layer() {
        let q = make_query(0.0, 0.0, 60_000.0);
        let r = evaluate_climate(&q);
        assert_eq!(r.layer, AtmosphericLayer::Mesosphere);
    }

    #[test]
    fn test_evaluate_thermosphere_layer() {
        let q = make_query(0.0, 0.0, 100_000.0);
        let r = evaluate_climate(&q);
        assert_eq!(r.layer, AtmosphericLayer::Thermosphere);
    }

    #[test]
    fn test_anomaly_location_stored_correctly() {
        let baseline = AtmosphericState {
            temperature_c: 20.0,
            pressure_hpa: 1013.25,
            humidity_pct: 60.0,
            wind_velocity_ms: [0.0; 3],
            density_kg_m3: 1.2,
        };
        let current = AtmosphericState {
            temperature_c: 50.0,
            ..baseline
        };
        let a = detect_anomaly(&current, &baseline, 35.5, 139.7, 12345).unwrap();
        assert!((a.location_lat - 35.5).abs() < 1e-10);
        assert!((a.location_lon - 139.7).abs() < 1e-10);
        assert_eq!(a.timestamp_ns, 12345);
    }

    #[test]
    fn test_fnv1a_deterministic_and_nonzero() {
        // Test the internal hash function through the anomaly API
        let baseline = AtmosphericState {
            temperature_c: 0.0,
            pressure_hpa: 1013.25,
            humidity_pct: 60.0,
            wind_velocity_ms: [0.0; 3],
            density_kg_m3: 1.2,
        };
        let current = AtmosphericState {
            temperature_c: 50.0,
            ..baseline
        };
        // Test at various locations
        for (lat, lon) in [(0.0, 0.0), (90.0, 180.0), (-45.0, -90.0)] {
            let a1 = detect_anomaly(&current, &baseline, lat, lon, 0).unwrap();
            let a2 = detect_anomaly(&current, &baseline, lat, lon, 0).unwrap();
            assert_eq!(a1.content_hash, a2.content_hash);
            assert_ne!(a1.content_hash, 0);
        }
    }
}
