# ALICE-Climate

Planetary climate as continuous SDF fields — infinite-resolution atmospheric and oceanic simulation.

## Features

| Feature | Description |
|---------|-------------|
| **Atmosphere** | ISA model, temperature/wind fields with latitude/seasonal correction |
| **Ocean** | Temperature, pressure, density (UNESCO), current systems |
| **Model** | Unified `evaluate_climate` with Rayon batch parallelism |
| **Anomaly** | Detection: HeatWave, ColdSnap, Storm, Drought, Flood |
| **Station** | Weather station network with Haversine distance |
| **FFI** | C-ABI for Unity / UE5 / C consumers (`al_clm_*`) |

## Quick Start

```rust
use alice_climate::model::{ClimateQuery, evaluate_climate};

let query = ClimateQuery {
    latitude: 35.6,
    longitude: 139.7,
    altitude_m: 0.0,
    timestamp_ns: 180u64 * 24 * 3600 * 1_000_000_000,
};
let response = evaluate_climate(&query);
assert!(response.atmosphere.temperature_c > -100.0);
```

## Build

```bash
cargo build --release
cargo test
```

### FFI (staticlib + cdylib)

```bash
cargo build --release --features ffi
```

## Test Suite

| Category | Tests |
|----------|-------|
| atmosphere | 26 |
| ocean | 21 |
| model | 37 |
| station | 18 |
| ffi | 16 |
| doc | 1 |
| **Total** | **111** |

## FFI Functions (17)

Prefix: `al_clm_*`

| Group | Count | Functions |
|-------|-------|-----------|
| Atmosphere | 5 | standard_atmosphere, temperature_field, wind_field, layer_from_altitude, layer_altitude_range |
| Ocean | 6 | ocean_temperature, ocean_pressure, ocean_density, ocean_current, ocean_layer_from_depth, ocean_layer_depth_range |
| Model | 3 | evaluate, evaluate_batch, detect_anomaly |
| Station | 2 | station_create, station_distance |
| Utility | 1 | version |

## Bindings

- `bindings/unity/AliceClimate.cs` — 17 DllImport + 7 structs
- `bindings/ue5/AliceClimate.h` — 17 extern C + 7 structs

## License

MIT — Copyright (C) 2026 Moroya Sakamoto
