# Changelog

All notable changes to ALICE-Climate will be documented in this file.

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
