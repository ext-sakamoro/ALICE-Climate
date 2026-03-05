//! C-ABI FFI for Unity / UE5 / C consumers
//!
//! All functions use the `al_clm_*` prefix (ALICE-Climate).
//!
//! ## Repr(C) Structs
//!
//! | Struct | Fields |
//! |--------|--------|
//! | `FfiAtmosphericState` | temperature_c, pressure_hpa, humidity_pct, wind_u/v/w, density_kg_m3 |
//! | `FfiOceanState` | temperature_c, salinity_psu, current_u/v/w, pressure_bar, density_kg_m3 |
//! | `FfiClimateQuery` | latitude, longitude, altitude_m, timestamp_ns |
//! | `FfiClimateResponse` | atmosphere, has_ocean, ocean, layer |
//! | `FfiClimateAnomaly` | location_lat/lon, kind, magnitude, timestamp_ns, content_hash |
//! | `FfiWeatherStation` | id, latitude, longitude, altitude_m, name_hash |
//! | `FfiObservation` | station_id, timestamp_ns, temperature_c, pressure_hpa, humidity_pct, wind_speed_ms, wind_direction_rad |

#![allow(clippy::missing_safety_doc)]

use core::ffi::{c_char, c_int};
use std::ffi::{CStr, CString};
use std::sync::OnceLock;

use crate::atmosphere::{standard_atmosphere, temperature_field, wind_field, AtmosphericLayer};
use crate::model::{detect_anomaly, evaluate_climate, evaluate_climate_batch, ClimateQuery};
use crate::ocean::{ocean_current, ocean_density, ocean_pressure, ocean_temperature, OceanLayer};

// ============================================================================
// Repr(C) Structs (7)
// ============================================================================

/// C-compatible atmospheric state.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct FfiAtmosphericState {
    pub temperature_c: f64,
    pub pressure_hpa: f64,
    pub humidity_pct: f64,
    pub wind_u: f64,
    pub wind_v: f64,
    pub wind_w: f64,
    pub density_kg_m3: f64,
}

/// C-compatible ocean state.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct FfiOceanState {
    pub temperature_c: f64,
    pub salinity_psu: f64,
    pub current_u: f64,
    pub current_v: f64,
    pub current_w: f64,
    pub pressure_bar: f64,
    pub density_kg_m3: f64,
}

/// C-compatible climate query.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct FfiClimateQuery {
    pub latitude: f64,
    pub longitude: f64,
    pub altitude_m: f64,
    pub timestamp_ns: u64,
}

/// C-compatible climate response.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct FfiClimateResponse {
    pub atmosphere: FfiAtmosphericState,
    /// 1 = ocean data present, 0 = above sea level
    pub has_ocean: u8,
    pub ocean: FfiOceanState,
    /// 0 = Troposphere, 1 = Stratosphere, 2 = Mesosphere, 3 = Thermosphere
    pub layer: u8,
}

/// C-compatible climate anomaly.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct FfiClimateAnomaly {
    pub location_lat: f64,
    pub location_lon: f64,
    /// 0 = HeatWave, 1 = ColdSnap, 2 = Storm, 3 = Drought, 4 = Flood
    pub kind: u8,
    pub magnitude: f64,
    pub timestamp_ns: u64,
    pub content_hash: u64,
}

/// C-compatible weather station.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct FfiWeatherStation {
    pub id: u64,
    pub latitude: f64,
    pub longitude: f64,
    pub altitude_m: f64,
    pub name_hash: u64,
}

/// C-compatible observation record.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct FfiObservation {
    pub station_id: u64,
    pub timestamp_ns: u64,
    pub temperature_c: f64,
    pub pressure_hpa: f64,
    pub humidity_pct: f64,
    pub wind_speed_ms: f64,
    pub wind_direction_rad: f64,
}

// ============================================================================
// Conversions
// ============================================================================

const ZERO_OCEAN: FfiOceanState = FfiOceanState {
    temperature_c: 0.0,
    salinity_psu: 0.0,
    current_u: 0.0,
    current_v: 0.0,
    current_w: 0.0,
    pressure_bar: 0.0,
    density_kg_m3: 0.0,
};

fn atmosphere_to_ffi(a: &crate::atmosphere::AtmosphericState) -> FfiAtmosphericState {
    FfiAtmosphericState {
        temperature_c: a.temperature_c,
        pressure_hpa: a.pressure_hpa,
        humidity_pct: a.humidity_pct,
        wind_u: a.wind_velocity_ms[0],
        wind_v: a.wind_velocity_ms[1],
        wind_w: a.wind_velocity_ms[2],
        density_kg_m3: a.density_kg_m3,
    }
}

fn ffi_to_atmospheric(f: &FfiAtmosphericState) -> crate::atmosphere::AtmosphericState {
    crate::atmosphere::AtmosphericState {
        temperature_c: f.temperature_c,
        pressure_hpa: f.pressure_hpa,
        humidity_pct: f.humidity_pct,
        wind_velocity_ms: [f.wind_u, f.wind_v, f.wind_w],
        density_kg_m3: f.density_kg_m3,
    }
}

fn layer_to_u8(l: AtmosphericLayer) -> u8 {
    match l {
        AtmosphericLayer::Troposphere => 0,
        AtmosphericLayer::Stratosphere => 1,
        AtmosphericLayer::Mesosphere => 2,
        AtmosphericLayer::Thermosphere => 3,
    }
}

fn ocean_layer_to_u8(l: OceanLayer) -> u8 {
    match l {
        OceanLayer::Surface => 0,
        OceanLayer::Thermocline => 1,
        OceanLayer::DeepOcean => 2,
        OceanLayer::Abyssal => 3,
    }
}

fn anomaly_kind_to_u8(k: crate::model::AnomalyKind) -> u8 {
    match k {
        crate::model::AnomalyKind::HeatWave => 0,
        crate::model::AnomalyKind::ColdSnap => 1,
        crate::model::AnomalyKind::Storm => 2,
        crate::model::AnomalyKind::Drought => 3,
        crate::model::AnomalyKind::Flood => 4,
    }
}

// ============================================================================
// Atmosphere (5)
// ============================================================================

/// Evaluate ISA at a given altitude. Writes result to `out`.
#[no_mangle]
pub unsafe extern "C" fn al_clm_standard_atmosphere(
    altitude_m: f64,
    out: *mut FfiAtmosphericState,
) {
    if out.is_null() {
        return;
    }
    *out = atmosphere_to_ffi(&standard_atmosphere(altitude_m));
}

/// Temperature (°C) with latitude/altitude/seasonal correction.
#[no_mangle]
pub extern "C" fn al_clm_temperature_field(
    lat: f64,
    lon: f64,
    alt_m: f64,
    day_of_year: u32,
) -> f64 {
    temperature_field(lat, lon, alt_m, day_of_year)
}

/// Wind velocity [u, v, w] in m/s. Writes 3 f64 values to `out_uvw`.
#[no_mangle]
pub unsafe extern "C" fn al_clm_wind_field(lat: f64, lon: f64, alt_m: f64, out_uvw: *mut f64) {
    if out_uvw.is_null() {
        return;
    }
    let [u, v, w] = wind_field(lat, lon, alt_m);
    *out_uvw = u;
    *out_uvw.add(1) = v;
    *out_uvw.add(2) = w;
}

/// Classify altitude to atmospheric layer.
/// Returns 0=Troposphere, 1=Stratosphere, 2=Mesosphere, 3=Thermosphere.
#[no_mangle]
pub extern "C" fn al_clm_layer_from_altitude(alt_m: f64) -> u8 {
    layer_to_u8(AtmosphericLayer::from_altitude(alt_m))
}

/// Altitude range (min, max) for a layer. Returns 1 if valid, 0 if invalid.
#[no_mangle]
pub unsafe extern "C" fn al_clm_layer_altitude_range(
    layer: u8,
    out_min: *mut f64,
    out_max: *mut f64,
) -> c_int {
    if out_min.is_null() || out_max.is_null() {
        return 0;
    }
    let l = match layer {
        0 => AtmosphericLayer::Troposphere,
        1 => AtmosphericLayer::Stratosphere,
        2 => AtmosphericLayer::Mesosphere,
        3 => AtmosphericLayer::Thermosphere,
        _ => return 0,
    };
    let (lo, hi) = l.altitude_range_m();
    *out_min = lo;
    *out_max = hi;
    1
}

// ============================================================================
// Ocean (6)
// ============================================================================

/// Ocean temperature (°C) at given depth and latitude.
#[no_mangle]
pub extern "C" fn al_clm_ocean_temperature(depth_m: f64, latitude: f64) -> f64 {
    ocean_temperature(depth_m, latitude)
}

/// Ocean pressure (bar) at given depth.
#[no_mangle]
pub extern "C" fn al_clm_ocean_pressure(depth_m: f64) -> f64 {
    ocean_pressure(depth_m)
}

/// Seawater density (kg/m³) from simplified UNESCO equation.
#[no_mangle]
pub extern "C" fn al_clm_ocean_density(
    temperature_c: f64,
    salinity_psu: f64,
    pressure_bar: f64,
) -> f64 {
    ocean_density(temperature_c, salinity_psu, pressure_bar)
}

/// Ocean current [u, v, w] in m/s. Writes 3 f64 values to `out_uvw`.
#[no_mangle]
pub unsafe extern "C" fn al_clm_ocean_current(lat: f64, lon: f64, depth_m: f64, out_uvw: *mut f64) {
    if out_uvw.is_null() {
        return;
    }
    let [u, v, w] = ocean_current(lat, lon, depth_m);
    *out_uvw = u;
    *out_uvw.add(1) = v;
    *out_uvw.add(2) = w;
}

/// Classify depth to ocean layer.
/// Returns 0=Surface, 1=Thermocline, 2=DeepOcean, 3=Abyssal.
#[no_mangle]
pub extern "C" fn al_clm_ocean_layer_from_depth(depth_m: f64) -> u8 {
    ocean_layer_to_u8(OceanLayer::from_depth(depth_m))
}

/// Depth range (min, max) for an ocean layer. Returns 1 if valid, 0 if invalid.
#[no_mangle]
pub unsafe extern "C" fn al_clm_ocean_layer_depth_range(
    layer: u8,
    out_min: *mut f64,
    out_max: *mut f64,
) -> c_int {
    if out_min.is_null() || out_max.is_null() {
        return 0;
    }
    let l = match layer {
        0 => OceanLayer::Surface,
        1 => OceanLayer::Thermocline,
        2 => OceanLayer::DeepOcean,
        3 => OceanLayer::Abyssal,
        _ => return 0,
    };
    let (lo, hi) = l.depth_range_m();
    *out_min = lo;
    *out_max = hi;
    1
}

// ============================================================================
// Model (3)
// ============================================================================

/// Evaluate climate at a single query point. Returns 1 on success, 0 on null.
#[no_mangle]
pub unsafe extern "C" fn al_clm_evaluate(
    query: *const FfiClimateQuery,
    out: *mut FfiClimateResponse,
) -> c_int {
    if query.is_null() || out.is_null() {
        return 0;
    }
    let q = ClimateQuery {
        latitude: (*query).latitude,
        longitude: (*query).longitude,
        altitude_m: (*query).altitude_m,
        timestamp_ns: (*query).timestamp_ns,
    };
    let r = evaluate_climate(&q);
    let (has_ocean, ocean_ffi) = match r.ocean {
        Some(o) => (
            1u8,
            FfiOceanState {
                temperature_c: o.temperature_c,
                salinity_psu: o.salinity_psu,
                current_u: o.current_velocity_ms[0],
                current_v: o.current_velocity_ms[1],
                current_w: o.current_velocity_ms[2],
                pressure_bar: o.pressure_bar,
                density_kg_m3: o.density_kg_m3,
            },
        ),
        None => (0u8, ZERO_OCEAN),
    };
    *out = FfiClimateResponse {
        atmosphere: atmosphere_to_ffi(&r.atmosphere),
        has_ocean,
        ocean: ocean_ffi,
        layer: layer_to_u8(r.layer),
    };
    1
}

/// Evaluate climate at a batch of query points. Caller allocates `out` array.
/// Returns number of results written.
#[no_mangle]
pub unsafe extern "C" fn al_clm_evaluate_batch(
    queries: *const FfiClimateQuery,
    count: c_int,
    out: *mut FfiClimateResponse,
) -> c_int {
    if queries.is_null() || out.is_null() || count <= 0 {
        return 0;
    }
    let n = count as usize;
    let rust_queries: Vec<ClimateQuery> = (0..n)
        .map(|i| {
            let q = &*queries.add(i);
            ClimateQuery {
                latitude: q.latitude,
                longitude: q.longitude,
                altitude_m: q.altitude_m,
                timestamp_ns: q.timestamp_ns,
            }
        })
        .collect();
    let results = evaluate_climate_batch(&rust_queries);
    for (i, r) in results.iter().enumerate() {
        let (has_ocean, ocean_ffi) = match r.ocean {
            Some(o) => (
                1u8,
                FfiOceanState {
                    temperature_c: o.temperature_c,
                    salinity_psu: o.salinity_psu,
                    current_u: o.current_velocity_ms[0],
                    current_v: o.current_velocity_ms[1],
                    current_w: o.current_velocity_ms[2],
                    pressure_bar: o.pressure_bar,
                    density_kg_m3: o.density_kg_m3,
                },
            ),
            None => (0u8, ZERO_OCEAN),
        };
        *out.add(i) = FfiClimateResponse {
            atmosphere: atmosphere_to_ffi(&r.atmosphere),
            has_ocean,
            ocean: ocean_ffi,
            layer: layer_to_u8(r.layer),
        };
    }
    count
}

/// Detect climate anomaly. Returns 1 if anomaly found, 0 otherwise.
#[no_mangle]
pub unsafe extern "C" fn al_clm_detect_anomaly(
    current: *const FfiAtmosphericState,
    baseline: *const FfiAtmosphericState,
    lat: f64,
    lon: f64,
    timestamp_ns: u64,
    out: *mut FfiClimateAnomaly,
) -> c_int {
    if current.is_null() || baseline.is_null() || out.is_null() {
        return 0;
    }
    let cur = ffi_to_atmospheric(&*current);
    let base = ffi_to_atmospheric(&*baseline);
    match detect_anomaly(&cur, &base, lat, lon, timestamp_ns) {
        Some(a) => {
            *out = FfiClimateAnomaly {
                location_lat: a.location_lat,
                location_lon: a.location_lon,
                kind: anomaly_kind_to_u8(a.kind),
                magnitude: a.magnitude,
                timestamp_ns: a.timestamp_ns,
                content_hash: a.content_hash,
            };
            1
        }
        None => 0,
    }
}

// ============================================================================
// Station (2)
// ============================================================================

/// Create a weather station from C string name. Writes to `out`.
/// Returns 1 on success, 0 on null/invalid UTF-8.
#[no_mangle]
pub unsafe extern "C" fn al_clm_station_create(
    id: u64,
    name: *const c_char,
    lat: f64,
    lon: f64,
    alt: f64,
    out: *mut FfiWeatherStation,
) -> c_int {
    if name.is_null() || out.is_null() {
        return 0;
    }
    let name_str = match CStr::from_ptr(name).to_str() {
        Ok(s) => s,
        Err(_) => return 0,
    };
    let station = crate::station::WeatherStation::new(id, name_str, lat, lon, alt);
    *out = FfiWeatherStation {
        id: station.id.0,
        latitude: station.latitude,
        longitude: station.longitude,
        altitude_m: station.altitude_m,
        name_hash: station.name_hash,
    };
    1
}

/// Haversine great-circle distance between two stations in km.
/// Returns -1.0 on null input.
#[no_mangle]
pub unsafe extern "C" fn al_clm_station_distance(
    a: *const FfiWeatherStation,
    b: *const FfiWeatherStation,
) -> f64 {
    if a.is_null() || b.is_null() {
        return -1.0;
    }
    let sa = crate::station::WeatherStation {
        id: crate::station::StationId((*a).id),
        latitude: (*a).latitude,
        longitude: (*a).longitude,
        altitude_m: (*a).altitude_m,
        name_hash: (*a).name_hash,
    };
    let sb = crate::station::WeatherStation {
        id: crate::station::StationId((*b).id),
        latitude: (*b).latitude,
        longitude: (*b).longitude,
        altitude_m: (*b).altitude_m,
        name_hash: (*b).name_hash,
    };
    sa.distance_to(&sb)
}

// ============================================================================
// Utility (1)
// ============================================================================

/// Library version string. Do NOT free the returned pointer.
#[no_mangle]
pub extern "C" fn al_clm_version() -> *const c_char {
    static VERSION: OnceLock<CString> = OnceLock::new();
    VERSION
        .get_or_init(|| {
            CString::new(env!("CARGO_PKG_VERSION"))
                .unwrap_or_else(|_| CString::new("0.0.0").unwrap())
        })
        .as_ptr()
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_atmosphere() {
        unsafe {
            let mut state = std::mem::zeroed::<FfiAtmosphericState>();
            al_clm_standard_atmosphere(0.0, &mut state);
            assert!((state.temperature_c - 15.0).abs() < 0.1);
            assert!((state.pressure_hpa - 1013.25).abs() < 0.5);
            assert!(state.density_kg_m3 > 1.0);
        }
    }

    #[test]
    fn test_temperature_field_ffi() {
        let t = al_clm_temperature_field(0.0, 0.0, 0.0, 180);
        assert!(t > -50.0 && t < 80.0);
    }

    #[test]
    fn test_wind_field_ffi() {
        unsafe {
            let mut uvw = [0.0f64; 3];
            al_clm_wind_field(45.0, 0.0, 0.0, uvw.as_mut_ptr());
            // Westerlies at 45N: positive u
            assert!(uvw[0] > 0.0);
        }
    }

    #[test]
    fn test_layer_classification() {
        assert_eq!(al_clm_layer_from_altitude(0.0), 0);
        assert_eq!(al_clm_layer_from_altitude(20_000.0), 1);
        assert_eq!(al_clm_layer_from_altitude(65_000.0), 2);
        assert_eq!(al_clm_layer_from_altitude(100_000.0), 3);
    }

    #[test]
    fn test_layer_altitude_range() {
        unsafe {
            let mut lo = 0.0f64;
            let mut hi = 0.0f64;
            assert_eq!(al_clm_layer_altitude_range(0, &mut lo, &mut hi), 1);
            assert_eq!(lo, 0.0);
            assert_eq!(hi, 12_000.0);

            // Invalid layer
            assert_eq!(al_clm_layer_altitude_range(99, &mut lo, &mut hi), 0);
        }
    }

    #[test]
    fn test_ocean_functions() {
        let t = al_clm_ocean_temperature(0.0, 0.0);
        assert!((t - 25.0).abs() < 1.0);

        let p = al_clm_ocean_pressure(1_000.0);
        assert!((p - 101.01325).abs() < 0.1);

        let d = al_clm_ocean_density(15.0, 35.0, 1.01325);
        assert!(d > 1020.0 && d < 1030.0);
    }

    #[test]
    fn test_ocean_current_ffi() {
        unsafe {
            let mut uvw = [0.0f64; 3];
            al_clm_ocean_current(0.0, 0.0, 0.0, uvw.as_mut_ptr());
            // Equatorial current: westward (negative u)
            assert!(uvw[0] < 0.0);
        }
    }

    #[test]
    fn test_ocean_layer() {
        assert_eq!(al_clm_ocean_layer_from_depth(100.0), 0);
        assert_eq!(al_clm_ocean_layer_from_depth(500.0), 1);
        assert_eq!(al_clm_ocean_layer_from_depth(2_000.0), 2);
        assert_eq!(al_clm_ocean_layer_from_depth(5_000.0), 3);

        unsafe {
            let mut lo = 0.0f64;
            let mut hi = 0.0f64;
            assert_eq!(al_clm_ocean_layer_depth_range(0, &mut lo, &mut hi), 1);
            assert_eq!(lo, 0.0);
            assert_eq!(hi, 200.0);
            assert_eq!(al_clm_ocean_layer_depth_range(99, &mut lo, &mut hi), 0);
        }
    }

    #[test]
    fn test_evaluate_surface() {
        unsafe {
            let query = FfiClimateQuery {
                latitude: 35.6,
                longitude: 139.7,
                altitude_m: 0.0,
                timestamp_ns: 180u64 * 24 * 3600 * 1_000_000_000,
            };
            let mut resp = std::mem::zeroed::<FfiClimateResponse>();
            assert_eq!(al_clm_evaluate(&query, &mut resp), 1);
            assert!(resp.atmosphere.pressure_hpa > 900.0);
            assert_eq!(resp.has_ocean, 0);
            assert_eq!(resp.layer, 0); // Troposphere
        }
    }

    #[test]
    fn test_evaluate_ocean() {
        unsafe {
            let query = FfiClimateQuery {
                latitude: 0.0,
                longitude: 0.0,
                altitude_m: -500.0,
                timestamp_ns: 0,
            };
            let mut resp = std::mem::zeroed::<FfiClimateResponse>();
            assert_eq!(al_clm_evaluate(&query, &mut resp), 1);
            assert_eq!(resp.has_ocean, 1);
            assert!(resp.ocean.pressure_bar > 40.0);
            assert!(resp.ocean.density_kg_m3 > 1000.0);
        }
    }

    #[test]
    fn test_evaluate_batch_ffi() {
        unsafe {
            let queries = [
                FfiClimateQuery {
                    latitude: 0.0,
                    longitude: 0.0,
                    altitude_m: 0.0,
                    timestamp_ns: 0,
                },
                FfiClimateQuery {
                    latitude: 45.0,
                    longitude: 90.0,
                    altitude_m: 5_000.0,
                    timestamp_ns: 0,
                },
            ];
            let mut results = [std::mem::zeroed::<FfiClimateResponse>(); 2];
            let n = al_clm_evaluate_batch(queries.as_ptr(), 2, results.as_mut_ptr());
            assert_eq!(n, 2);
            assert!(results[0].atmosphere.pressure_hpa > 900.0);
            assert!(results[1].atmosphere.pressure_hpa < results[0].atmosphere.pressure_hpa);
        }
    }

    #[test]
    fn test_detect_anomaly_heat() {
        unsafe {
            let baseline = FfiAtmosphericState {
                temperature_c: 20.0,
                pressure_hpa: 1013.25,
                humidity_pct: 60.0,
                wind_u: 0.0,
                wind_v: 0.0,
                wind_w: 0.0,
                density_kg_m3: 1.2,
            };
            let current = FfiAtmosphericState {
                temperature_c: 45.0,
                ..baseline
            };
            let mut anomaly = std::mem::zeroed::<FfiClimateAnomaly>();
            assert_eq!(
                al_clm_detect_anomaly(&current, &baseline, 35.0, 139.0, 0, &mut anomaly),
                1
            );
            assert_eq!(anomaly.kind, 0); // HeatWave
            assert!((anomaly.magnitude - 25.0).abs() < 0.01);
        }
    }

    #[test]
    fn test_detect_anomaly_none() {
        unsafe {
            let state = FfiAtmosphericState {
                temperature_c: 20.0,
                pressure_hpa: 1013.25,
                humidity_pct: 60.0,
                wind_u: 3.0,
                wind_v: 1.0,
                wind_w: 0.0,
                density_kg_m3: 1.2,
            };
            let mut anomaly = std::mem::zeroed::<FfiClimateAnomaly>();
            assert_eq!(
                al_clm_detect_anomaly(&state, &state, 0.0, 0.0, 0, &mut anomaly),
                0
            );
        }
    }

    #[test]
    fn test_station_create_and_distance() {
        unsafe {
            let name_tokyo = std::ffi::CString::new("Tokyo").unwrap();
            let name_osaka = std::ffi::CString::new("Osaka").unwrap();
            let mut tokyo = std::mem::zeroed::<FfiWeatherStation>();
            let mut osaka = std::mem::zeroed::<FfiWeatherStation>();

            assert_eq!(
                al_clm_station_create(1, name_tokyo.as_ptr(), 35.6895, 139.6917, 40.0, &mut tokyo),
                1
            );
            assert_eq!(
                al_clm_station_create(2, name_osaka.as_ptr(), 34.6937, 135.5023, 17.0, &mut osaka),
                1
            );

            let dist = al_clm_station_distance(&tokyo, &osaka);
            assert!(
                (dist - 396.0).abs() < 10.0,
                "Expected ~396 km, got {:.1}",
                dist
            );
        }
    }

    #[test]
    fn test_null_safety() {
        unsafe {
            al_clm_standard_atmosphere(0.0, std::ptr::null_mut());
            al_clm_wind_field(0.0, 0.0, 0.0, std::ptr::null_mut());
            al_clm_ocean_current(0.0, 0.0, 0.0, std::ptr::null_mut());

            assert_eq!(
                al_clm_layer_altitude_range(0, std::ptr::null_mut(), std::ptr::null_mut()),
                0
            );
            assert_eq!(
                al_clm_ocean_layer_depth_range(0, std::ptr::null_mut(), std::ptr::null_mut()),
                0
            );
            assert_eq!(al_clm_evaluate(std::ptr::null(), std::ptr::null_mut()), 0);
            assert_eq!(
                al_clm_evaluate_batch(std::ptr::null(), 0, std::ptr::null_mut()),
                0
            );
            assert_eq!(
                al_clm_detect_anomaly(
                    std::ptr::null(),
                    std::ptr::null(),
                    0.0,
                    0.0,
                    0,
                    std::ptr::null_mut()
                ),
                0
            );
            assert_eq!(
                al_clm_station_create(0, std::ptr::null(), 0.0, 0.0, 0.0, std::ptr::null_mut()),
                0
            );
            assert_eq!(
                al_clm_station_distance(std::ptr::null(), std::ptr::null()),
                -1.0
            );
        }
    }

    #[test]
    fn test_version() {
        let ptr = al_clm_version();
        assert!(!ptr.is_null());
        let v = unsafe { std::ffi::CStr::from_ptr(ptr) }.to_str().unwrap();
        assert!(v.starts_with("0."));
    }
}
