#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss
)]

//! # ALICE-Climate
//!
//! Planetary climate as continuous SDF fields, enabling infinite-resolution
//! queries of temperature, wind, and humidity at any point and time.
//!
//! Instead of discrete grid voxels, atmospheric and oceanic states are
//! represented as continuous mathematical fields that can be sampled at
//! arbitrary resolution.
//!
//! ## Modules
//!
//! - [`atmosphere`] — Atmospheric layers and continuous field evaluation
//! - [`ocean`] — Ocean layers, temperature, salinity, and current models
//! - [`model`] — Unified climate model integrating atmosphere and ocean
//! - [`station`] — Weather station observation data structures
//!
//! ## Quick Example
//!
//! ```rust
//! use alice_climate::model::{ClimateQuery, evaluate_climate};
//!
//! // Query climate at Tokyo (35.6°N, 139.7°E), surface level, day 180
//! let query = ClimateQuery {
//!     latitude: 35.6,
//!     longitude: 139.7,
//!     altitude_m: 0.0,
//!     timestamp_ns: 180u64 * 24 * 3600 * 1_000_000_000,
//! };
//! let response = evaluate_climate(&query);
//! assert!(response.atmosphere.temperature_c > -100.0);
//! assert!(response.atmosphere.pressure_hpa > 0.0);
//! ```

pub mod atmosphere;
pub mod coupling;
#[cfg(feature = "ffi")]
pub mod ffi;
pub mod model;
pub mod ocean;
pub mod simulation;
pub mod solar;
pub mod station;

pub use atmosphere::{
    standard_atmosphere, temperature_field, wind_field, AtmosphericLayer, AtmosphericState,
};
pub use model::{
    detect_anomaly, evaluate_climate, evaluate_climate_batch, AnomalyKind, ClimateAnomaly,
    ClimateQuery, ClimateResponse,
};
pub use ocean::{
    ocean_current, ocean_density, ocean_pressure, ocean_temperature, OceanLayer, OceanState,
};
pub use station::{Observation, StationId, WeatherStation};

pub use coupling::{
    compute_coupled_state, coupled_step, evaporation_rate, latent_heat_flux, net_heat_flux,
    sensible_heat_flux, BulkCoefficients, CoupledState,
};
pub use simulation::{SimulationConfig, SimulationState, StepResult};
pub use solar::{
    cos_solar_zenith, diurnal_temperature_correction, insolation, seasonal_temperature_correction,
    solar_declination, solar_hour_angle, solar_zenith_angle, surface_insolation, SOLAR_CONSTANT,
};
