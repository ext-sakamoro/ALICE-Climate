/*
    ALICE-Climate — UE5 C++ Bindings
    Copyright (C) 2026 Moroya Sakamoto

    17 extern C + 7 structs
    Prefix: al_clm_*
*/

#pragma once

#include <cstdint>

// ============================================================================
// Repr(C) structs (7)
// ============================================================================

struct AlClmAtmosphericState {
    double temperature_c;
    double pressure_hpa;
    double humidity_pct;
    double wind_u;
    double wind_v;
    double wind_w;
    double density_kg_m3;
};

struct AlClmOceanState {
    double temperature_c;
    double salinity_psu;
    double current_u;
    double current_v;
    double current_w;
    double pressure_bar;
    double density_kg_m3;
};

struct AlClmClimateQuery {
    double latitude;
    double longitude;
    double altitude_m;
    uint64_t timestamp_ns;
};

struct AlClmClimateResponse {
    AlClmAtmosphericState atmosphere;
    uint8_t has_ocean;      // 1 = ocean data present
    AlClmOceanState ocean;
    uint8_t layer;          // 0=Troposphere, 1=Stratosphere, 2=Mesosphere, 3=Thermosphere
};

struct AlClmClimateAnomaly {
    double location_lat;
    double location_lon;
    uint8_t kind;           // 0=HeatWave, 1=ColdSnap, 2=Storm, 3=Drought, 4=Flood
    double magnitude;
    uint64_t timestamp_ns;
    uint64_t content_hash;
};

struct AlClmWeatherStation {
    uint64_t id;
    double latitude;
    double longitude;
    double altitude_m;
    uint64_t name_hash;
};

struct AlClmObservation {
    uint64_t station_id;
    uint64_t timestamp_ns;
    double temperature_c;
    double pressure_hpa;
    double humidity_pct;
    double wind_speed_ms;
    double wind_direction_rad;
};

// ============================================================================
// extern "C" (17 functions)
// ============================================================================

extern "C" {

// --- Atmosphere (5) ---
void al_clm_standard_atmosphere(double altitude_m, AlClmAtmosphericState* out);
double al_clm_temperature_field(double lat, double lon, double alt_m, uint32_t day_of_year);
void al_clm_wind_field(double lat, double lon, double alt_m, double* out_uvw);
uint8_t al_clm_layer_from_altitude(double alt_m);
int32_t al_clm_layer_altitude_range(uint8_t layer, double* out_min, double* out_max);

// --- Ocean (6) ---
double al_clm_ocean_temperature(double depth_m, double latitude);
double al_clm_ocean_pressure(double depth_m);
double al_clm_ocean_density(double temperature_c, double salinity_psu, double pressure_bar);
void al_clm_ocean_current(double lat, double lon, double depth_m, double* out_uvw);
uint8_t al_clm_ocean_layer_from_depth(double depth_m);
int32_t al_clm_ocean_layer_depth_range(uint8_t layer, double* out_min, double* out_max);

// --- Model (3) ---
int32_t al_clm_evaluate(const AlClmClimateQuery* query, AlClmClimateResponse* out);
int32_t al_clm_evaluate_batch(
    const AlClmClimateQuery* queries, int32_t count, AlClmClimateResponse* out);
int32_t al_clm_detect_anomaly(
    const AlClmAtmosphericState* current, const AlClmAtmosphericState* baseline,
    double lat, double lon, uint64_t timestamp_ns, AlClmClimateAnomaly* out);

// --- Station (2) ---
int32_t al_clm_station_create(
    uint64_t id, const char* name, double lat, double lon, double alt,
    AlClmWeatherStation* out);
double al_clm_station_distance(const AlClmWeatherStation* a, const AlClmWeatherStation* b);

// --- Utility (1) ---
const char* al_clm_version();

} // extern "C"
