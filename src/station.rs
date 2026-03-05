//! Weather station and observation data structures.
//!
//! Provides types for representing physical weather stations and
//! their recorded observations, including haversine-based distance
//! calculations between stations.

/// Unique identifier for a weather station.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct StationId(pub u64);

/// A single weather observation recorded at a station.
#[derive(Debug, Clone, Copy)]
pub struct Observation {
    /// The station that recorded this observation.
    pub station_id: StationId,
    /// Timestamp in nanoseconds since Unix epoch.
    pub timestamp_ns: u64,
    /// Ambient temperature in degrees Celsius.
    pub temperature_c: f64,
    /// Atmospheric pressure in hectopascals.
    pub pressure_hpa: f64,
    /// Relative humidity as a percentage (0–100).
    pub humidity_pct: f64,
    /// Wind speed in metres per second.
    pub wind_speed_ms: f64,
    /// Wind direction in radians, measured clockwise from north.
    pub wind_direction_rad: f64,
}

/// A physical weather station at a fixed geographic location.
#[derive(Debug, Clone)]
pub struct WeatherStation {
    /// Unique identifier.
    pub id: StationId,
    /// Latitude in decimal degrees (−90 to +90).
    pub latitude: f64,
    /// Longitude in decimal degrees (−180 to +180).
    pub longitude: f64,
    /// Altitude above sea level in metres.
    pub altitude_m: f64,
    /// FNV-1a hash of the station name for fast equality checks.
    pub name_hash: u64,
}

/// FNV-1a 64-bit hash of a byte slice.
#[inline]
fn fnv1a(data: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in data {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

impl WeatherStation {
    /// Create a new weather station.
    ///
    /// The station name is hashed with FNV-1a and stored as `name_hash`.
    ///
    /// # Arguments
    ///
    /// * `id` — Raw numeric identifier.
    /// * `name` — Human-readable station name (hashed internally).
    /// * `lat` — Latitude in decimal degrees.
    /// * `lon` — Longitude in decimal degrees.
    /// * `alt` — Altitude above sea level in metres.
    #[must_use] 
    pub fn new(id: u64, name: &str, lat: f64, lon: f64, alt: f64) -> Self {
        Self {
            id: StationId(id),
            latitude: lat,
            longitude: lon,
            altitude_m: alt,
            name_hash: fnv1a(name.as_bytes()),
        }
    }

    /// Compute the haversine great-circle distance to another station in km.
    ///
    /// Uses Earth radius R = 6 371 km.
    #[must_use] 
    pub fn distance_to(&self, other: &WeatherStation) -> f64 {
        const R_KM: f64 = 6_371.0;
        let lat1 = self.latitude.to_radians();
        let lat2 = other.latitude.to_radians();
        let dlat = (other.latitude - self.latitude).to_radians();
        let dlon = (other.longitude - self.longitude).to_radians();

        let sin_dlat_half = (dlat * 0.5).sin();
        let sin_dlon_half = (dlon * 0.5).sin();
        let a = sin_dlat_half * sin_dlat_half
            + lat1.cos() * lat2.cos() * (sin_dlon_half * sin_dlon_half);
        let c = 2.0 * a.sqrt().atan2((1.0 - a).sqrt());
        R_KM * c
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_station() {
        let s = WeatherStation::new(1, "Tokyo", 35.6895, 139.6917, 40.0);
        assert_eq!(s.id, StationId(1));
        assert!((s.latitude - 35.6895).abs() < 1e-10);
        assert!((s.longitude - 139.6917).abs() < 1e-10);
        assert!((s.altitude_m - 40.0).abs() < 1e-10);
    }

    #[test]
    fn test_haversine_distance_known_pair() {
        // Tokyo (35.6895°N, 139.6917°E) to Osaka (34.6937°N, 135.5023°E)
        // Haversine great-circle distance ≈ 396 km for these exact coordinates
        let tokyo = WeatherStation::new(1, "Tokyo", 35.6895, 139.6917, 40.0);
        let osaka = WeatherStation::new(2, "Osaka", 34.6937, 135.5023, 17.0);
        let dist = tokyo.distance_to(&osaka);
        // Allow ±10 km tolerance around 396 km
        assert!(
            (dist - 396.0).abs() < 10.0,
            "Expected ~396 km, got {:.1} km",
            dist
        );
    }

    #[test]
    fn test_haversine_distance_same_station_is_zero() {
        let s = WeatherStation::new(1, "Tokyo", 35.6895, 139.6917, 40.0);
        assert!(s.distance_to(&s).abs() < 1e-9);
    }

    #[test]
    fn test_name_hash_determinism() {
        let s1 = WeatherStation::new(1, "TestStation", 0.0, 0.0, 0.0);
        let s2 = WeatherStation::new(2, "TestStation", 10.0, 20.0, 100.0);
        // Same name -> same hash regardless of other fields
        assert_eq!(s1.name_hash, s2.name_hash);
    }

    #[test]
    fn test_name_hash_different_names() {
        let s1 = WeatherStation::new(1, "Alpha", 0.0, 0.0, 0.0);
        let s2 = WeatherStation::new(2, "Beta", 0.0, 0.0, 0.0);
        assert_ne!(s1.name_hash, s2.name_hash);
    }

    #[test]
    fn test_observation_creation() {
        let obs = Observation {
            station_id: StationId(42),
            timestamp_ns: 1_700_000_000_000_000_000,
            temperature_c: 22.5,
            pressure_hpa: 1013.25,
            humidity_pct: 65.0,
            wind_speed_ms: 5.2,
            wind_direction_rad: 1.57,
        };
        assert_eq!(obs.station_id, StationId(42));
        assert!((obs.temperature_c - 22.5).abs() < 1e-10);
        assert!((obs.pressure_hpa - 1013.25).abs() < 1e-10);
    }

    #[test]
    fn test_haversine_antipodal_cities() {
        // London (~51.5°N, 0°) to Sydney (~33.9°S, 151.2°E)
        // Expected ~16 993 km
        let london = WeatherStation::new(1, "London", 51.5074, -0.1278, 11.0);
        let sydney = WeatherStation::new(2, "Sydney", -33.8688, 151.2093, 39.0);
        let dist = london.distance_to(&sydney);
        assert!(
            (dist - 16_993.0).abs() < 50.0,
            "Expected ~16993 km, got {:.1} km",
            dist
        );
    }

    // ── Additional tests for improved coverage ──────────────────────

    #[test]
    fn test_haversine_symmetry() {
        let a = WeatherStation::new(1, "A", 35.6895, 139.6917, 0.0);
        let b = WeatherStation::new(2, "B", 34.6937, 135.5023, 0.0);
        let d_ab = a.distance_to(&b);
        let d_ba = b.distance_to(&a);
        assert!(
            (d_ab - d_ba).abs() < 1e-6,
            "Distance should be symmetric: {:.4} vs {:.4}",
            d_ab,
            d_ba
        );
    }

    #[test]
    fn test_haversine_equatorial_points() {
        // Two points on the equator separated by 1 degree of longitude
        // At equator, 1 degree ~= 111.32 km
        let a = WeatherStation::new(1, "A", 0.0, 0.0, 0.0);
        let b = WeatherStation::new(2, "B", 0.0, 1.0, 0.0);
        let dist = a.distance_to(&b);
        assert!(
            (dist - 111.32).abs() < 1.0,
            "1 degree at equator should be ~111.32 km, got {:.2} km",
            dist
        );
    }

    #[test]
    fn test_haversine_poles_to_equator() {
        // North pole to equator should be ~10 018 km (quarter circumference)
        let pole = WeatherStation::new(1, "NP", 90.0, 0.0, 0.0);
        let equator = WeatherStation::new(2, "EQ", 0.0, 0.0, 0.0);
        let dist = pole.distance_to(&equator);
        assert!(
            (dist - 10_018.0).abs() < 30.0,
            "Pole to equator should be ~10018 km, got {:.1} km",
            dist
        );
    }

    #[test]
    fn test_station_id_equality() {
        assert_eq!(StationId(42), StationId(42));
        assert_ne!(StationId(1), StationId(2));
    }

    #[test]
    fn test_station_different_ids_same_location() {
        let s1 = WeatherStation::new(1, "StationA", 35.0, 139.0, 10.0);
        let s2 = WeatherStation::new(2, "StationA", 35.0, 139.0, 10.0);
        assert_ne!(s1.id, s2.id);
        assert_eq!(s1.name_hash, s2.name_hash);
        // Distance between same location should be ~0
        assert!(s1.distance_to(&s2).abs() < 1e-9);
    }

    #[test]
    fn test_name_hash_empty_string() {
        let s = WeatherStation::new(1, "", 0.0, 0.0, 0.0);
        // Even empty string should produce a valid non-zero hash (FNV basis)
        assert_ne!(s.name_hash, 0);
    }

    #[test]
    fn test_observation_field_access() {
        let obs = Observation {
            station_id: StationId(100),
            timestamp_ns: 1_000_000_000_000_000_000,
            temperature_c: -40.0,
            pressure_hpa: 980.0,
            humidity_pct: 0.0,
            wind_speed_ms: 50.0,
            wind_direction_rad: std::f64::consts::PI,
        };
        assert_eq!(obs.station_id, StationId(100));
        assert!((obs.temperature_c - (-40.0)).abs() < 1e-10);
        assert!((obs.humidity_pct - 0.0).abs() < 1e-10);
        assert!((obs.wind_speed_ms - 50.0).abs() < 1e-10);
        assert!((obs.wind_direction_rad - std::f64::consts::PI).abs() < 1e-10);
    }

    #[test]
    fn test_station_negative_altitude() {
        // Weather station below sea level (e.g., Dead Sea)
        let s = WeatherStation::new(1, "DeadSea", 31.5, 35.5, -430.0);
        assert!((s.altitude_m - (-430.0)).abs() < 1e-10);
    }
}
