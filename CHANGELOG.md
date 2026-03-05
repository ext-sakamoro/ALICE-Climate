# Changelog

All notable changes to ALICE-Climate will be documented in this file.

## [0.1.1] - 2026-03-04

### Added
- `ffi` — 17 `al_clm_*` extern "C" functions (Atmosphere 5, Ocean 6, Model 3, Station 2, Utility 1)
- `ffi` — 7 `#[repr(C)]` structs: `FfiAtmosphericState`, `FfiOceanState`, `FfiClimateQuery`, `FfiClimateResponse`, `FfiClimateAnomaly`, `FfiWeatherStation`, `FfiObservation`
- `ffi` — 16 FFI tests (lifecycle, null safety, batch, anomaly)
- `bindings/unity/AliceClimate.cs` — 17 DllImport + 7 structs
- `bindings/ue5/AliceClimate.h` — 17 extern C + 7 structs
- `Cargo.toml` — `crate-type = ["rlib", "staticlib", "cdylib"]`, `ffi` feature
- `README.md`

### Changed
- Total tests: 95 → 111

## [0.1.0] - 2026-02-23

### Added
- `atmosphere` — continuous atmospheric SDF fields (temperature, wind, pressure, humidity)
- `AtmosphericLayer` enum and `standard_atmosphere` model
- `ocean` — ocean layer model (temperature, salinity, currents, pressure, density)
- `OceanLayer` enum with thermocline / halocline profiles
- `model` — unified `ClimateQuery` / `ClimateResponse` with `evaluate_climate` and batch evaluation
- `detect_anomaly` — climate anomaly detection (`AnomalyKind` enum)
- `station` — `WeatherStation` / `Observation` data structures with `StationId`
- 95 tests (94 unit + 1 doc-test)
