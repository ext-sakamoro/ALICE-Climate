/*
    ALICE-Climate — Unity C# Bindings
    Copyright (C) 2026 Moroya Sakamoto

    17 DllImport + 7 structs
    Prefix: al_clm_*
*/

using System;
using System.Runtime.InteropServices;

namespace Alice.Climate
{
    // ========================================================================
    // Repr(C) structs (7)
    // ========================================================================

    [StructLayout(LayoutKind.Sequential)]
    public struct FfiAtmosphericState
    {
        public double temperatureC;
        public double pressureHpa;
        public double humidityPct;
        public double windU;
        public double windV;
        public double windW;
        public double densityKgM3;
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct FfiOceanState
    {
        public double temperatureC;
        public double salinityPsu;
        public double currentU;
        public double currentV;
        public double currentW;
        public double pressureBar;
        public double densityKgM3;
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct FfiClimateQuery
    {
        public double latitude;
        public double longitude;
        public double altitudeM;
        public ulong timestampNs;
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct FfiClimateResponse
    {
        public FfiAtmosphericState atmosphere;
        public byte hasOcean;
        public FfiOceanState ocean;
        public byte layer;
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct FfiClimateAnomaly
    {
        public double locationLat;
        public double locationLon;
        public byte kind;
        public double magnitude;
        public ulong timestampNs;
        public ulong contentHash;
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct FfiWeatherStation
    {
        public ulong id;
        public double latitude;
        public double longitude;
        public double altitudeM;
        public ulong nameHash;
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct FfiObservation
    {
        public ulong stationId;
        public ulong timestampNs;
        public double temperatureC;
        public double pressureHpa;
        public double humidityPct;
        public double windSpeedMs;
        public double windDirectionRad;
    }

    // ========================================================================
    // Native (17 DllImport)
    // ========================================================================

    public static class Native
    {
        const string Lib = "alice_climate";

        // --- Atmosphere (5) ---

        [DllImport(Lib)]
        public static extern void al_clm_standard_atmosphere(
            double altitudeM, out FfiAtmosphericState outState);

        [DllImport(Lib)]
        public static extern double al_clm_temperature_field(
            double lat, double lon, double altM, uint dayOfYear);

        [DllImport(Lib)]
        public static extern unsafe void al_clm_wind_field(
            double lat, double lon, double altM, double* outUvw);

        [DllImport(Lib)]
        public static extern byte al_clm_layer_from_altitude(double altM);

        [DllImport(Lib)]
        public static extern int al_clm_layer_altitude_range(
            byte layer, out double outMin, out double outMax);

        // --- Ocean (6) ---

        [DllImport(Lib)]
        public static extern double al_clm_ocean_temperature(double depthM, double latitude);

        [DllImport(Lib)]
        public static extern double al_clm_ocean_pressure(double depthM);

        [DllImport(Lib)]
        public static extern double al_clm_ocean_density(
            double temperatureC, double salinityPsu, double pressureBar);

        [DllImport(Lib)]
        public static extern unsafe void al_clm_ocean_current(
            double lat, double lon, double depthM, double* outUvw);

        [DllImport(Lib)]
        public static extern byte al_clm_ocean_layer_from_depth(double depthM);

        [DllImport(Lib)]
        public static extern int al_clm_ocean_layer_depth_range(
            byte layer, out double outMin, out double outMax);

        // --- Model (3) ---

        [DllImport(Lib)]
        public static extern int al_clm_evaluate(
            ref FfiClimateQuery query, out FfiClimateResponse outResp);

        [DllImport(Lib)]
        public static extern unsafe int al_clm_evaluate_batch(
            FfiClimateQuery* queries, int count, FfiClimateResponse* outResp);

        [DllImport(Lib)]
        public static extern int al_clm_detect_anomaly(
            ref FfiAtmosphericState current, ref FfiAtmosphericState baseline,
            double lat, double lon, ulong timestampNs, out FfiClimateAnomaly outAnomaly);

        // --- Station (2) ---

        [DllImport(Lib)]
        public static extern int al_clm_station_create(
            ulong id, [MarshalAs(UnmanagedType.LPStr)] string name,
            double lat, double lon, double alt, out FfiWeatherStation outStation);

        [DllImport(Lib)]
        public static extern double al_clm_station_distance(
            ref FfiWeatherStation a, ref FfiWeatherStation b);

        // --- Utility (1) ---

        [DllImport(Lib)] public static extern IntPtr al_clm_version();

        // --- Helpers ---

        public static string Version()
        {
            IntPtr ptr = al_clm_version();
            return ptr != IntPtr.Zero ? Marshal.PtrToStringAnsi(ptr) : "unknown";
        }
    }
}
