//! 大気-海洋結合フィードバック — 界面での熱・水蒸気交換
//!
//! 海面バルク公式に基づく熱フラックス・蒸発量の計算と、
//! 大気-海洋を結合して更新する step 関数を提供する。

use crate::atmosphere::AtmosphericState;
use crate::ocean::OceanState;

/// 大気-海洋界面の状態。
#[derive(Debug, Clone, Copy)]
pub struct CoupledState {
    /// 大気の状態。
    pub atmosphere: AtmosphericState,
    /// 海洋の状態。
    pub ocean: OceanState,
    /// 海面での正味熱フラックス (W/m², 正 = 海洋→大気)。
    pub net_heat_flux_wm2: f64,
    /// 蒸発量 (kg/m²/s)。
    pub evaporation_rate: f64,
}

/// バルク公式の係数。
#[derive(Debug, Clone, Copy)]
pub struct BulkCoefficients {
    /// 顕熱交換係数 (Stanton 数, 典型値 ~1.0e-3)。
    pub c_h: f64,
    /// 潜熱交換係数 (Dalton 数, 典型値 ~1.2e-3)。
    pub c_e: f64,
    /// 空気密度 (kg/m³, 海面付近)。
    pub rho_air: f64,
    /// 空気の比熱 (J/(kg·K))。
    pub cp_air: f64,
    /// 水の蒸発潜熱 (J/kg)。
    pub latent_heat: f64,
}

impl Default for BulkCoefficients {
    fn default() -> Self {
        Self {
            c_h: 1.0e-3,
            c_e: 1.2e-3,
            rho_air: 1.225,
            cp_air: 1005.0,
            latent_heat: 2.45e6,
        }
    }
}

/// 海面風速 (m/s) を大気状態から計算。
fn surface_wind_speed(atm: &AtmosphericState) -> f64 {
    let u = atm.wind_velocity_ms[0];
    let v = atm.wind_velocity_ms[1];
    u.hypot(v).max(0.5) // 最低 0.5 m/s (calm 条件回避)
}

/// 飽和比湿 (kg/kg) を温度から近似 (Tetens 式)。
fn saturation_specific_humidity(temp_c: f64) -> f64 {
    // 飽和蒸気圧 (hPa): e_s = 6.112 * exp(17.67 * T / (T + 243.5))
    let e_s = 6.112 * (17.67 * temp_c / (temp_c + 243.5)).exp();
    // 比湿 ≈ 0.622 * e_s / (1013.25 - 0.378 * e_s)
    0.622 * e_s / 0.378f64.mul_add(-e_s, 1013.25)
}

/// 顕熱フラックス (W/m²)。
///
/// `Q_H` = ρ × `C_H` × `C_p` × U × (`T_sea` - `T_air`)
///
/// 正 = 海洋→大気 (海が空気より暖かい)。
#[must_use]
pub fn sensible_heat_flux(
    atm: &AtmosphericState,
    ocean: &OceanState,
    coeff: &BulkCoefficients,
) -> f64 {
    let wind = surface_wind_speed(atm);
    let dt = ocean.temperature_c - atm.temperature_c;
    coeff.rho_air * coeff.c_h * coeff.cp_air * wind * dt
}

/// 潜熱フラックス (W/m²)。
///
/// `Q_E` = ρ × `C_E` × L × U × (`q_s(T_sea)` - `q_a`)
///
/// 正 = 海洋→大気 (蒸発)。
#[must_use]
pub fn latent_heat_flux(
    atm: &AtmosphericState,
    ocean: &OceanState,
    coeff: &BulkCoefficients,
) -> f64 {
    let wind = surface_wind_speed(atm);
    let q_sea = saturation_specific_humidity(ocean.temperature_c);
    let q_air = saturation_specific_humidity(atm.temperature_c) * (atm.humidity_pct / 100.0);
    let dq = q_sea - q_air;
    coeff.rho_air * coeff.c_e * coeff.latent_heat * wind * dq
}

/// 蒸発量 (kg/m²/s)。
#[must_use]
pub fn evaporation_rate(
    atm: &AtmosphericState,
    ocean: &OceanState,
    coeff: &BulkCoefficients,
) -> f64 {
    let wind = surface_wind_speed(atm);
    let q_sea = saturation_specific_humidity(ocean.temperature_c);
    let q_air = saturation_specific_humidity(atm.temperature_c) * (atm.humidity_pct / 100.0);
    let dq = (q_sea - q_air).max(0.0);
    coeff.rho_air * coeff.c_e * wind * dq
}

/// 正味海面熱フラックス (W/m²)。
///
/// 顕熱 + 潜熱。正 = 海洋→大気。
#[must_use]
pub fn net_heat_flux(atm: &AtmosphericState, ocean: &OceanState, coeff: &BulkCoefficients) -> f64 {
    sensible_heat_flux(atm, ocean, coeff) + latent_heat_flux(atm, ocean, coeff)
}

/// 結合状態を計算。
#[must_use]
pub fn compute_coupled_state(atmosphere: AtmosphericState, ocean: OceanState) -> CoupledState {
    let coeff = BulkCoefficients::default();
    let flux = net_heat_flux(&atmosphere, &ocean, &coeff);
    let evap = evaporation_rate(&atmosphere, &ocean, &coeff);

    CoupledState {
        atmosphere,
        ocean,
        net_heat_flux_wm2: flux,
        evaporation_rate: evap,
    }
}

/// 大気-海洋結合ステップ。
///
/// `dt` 秒間のフラックスに基づいて大気温度・湿度と海面温度を更新。
/// 海洋の混合層深度 `mixed_layer_depth_m` で熱容量を決定。
pub fn coupled_step(
    atm: &mut AtmosphericState,
    ocean: &mut OceanState,
    dt_s: f64,
    mixed_layer_depth_m: f64,
) {
    let coeff = BulkCoefficients::default();
    let q_h = sensible_heat_flux(atm, ocean, &coeff);
    let q_e = latent_heat_flux(atm, ocean, &coeff);
    let evap = evaporation_rate(atm, ocean, &coeff);

    // 海洋: 混合層の熱容量 (J/(m²·K))
    // ρ_water × C_p_water × depth
    let rho_water = 1025.0;
    let cp_water = 3985.0;
    let ocean_heat_capacity = rho_water * cp_water * mixed_layer_depth_m.max(1.0);

    // 海面温度変化: ΔT_ocean = -(Q_H + Q_E) × dt / heat_capacity
    let dt_ocean = -(q_h + q_e) * dt_s / ocean_heat_capacity;
    ocean.temperature_c += dt_ocean;

    // 大気温度変化 (境界層 ~1000m)
    let atm_layer_depth = 1000.0;
    let atm_heat_capacity = coeff.rho_air * coeff.cp_air * atm_layer_depth;
    let dt_atm = (q_h + q_e) * dt_s / atm_heat_capacity;
    atm.temperature_c += dt_atm;

    // 大気湿度変化: 蒸発による加湿
    let moisture_addition = evap * dt_s / (coeff.rho_air * atm_layer_depth);
    let q_current = saturation_specific_humidity(atm.temperature_c) * (atm.humidity_pct / 100.0);
    let q_new = q_current + moisture_addition;
    let q_sat = saturation_specific_humidity(atm.temperature_c);
    atm.humidity_pct = if q_sat > 0.0 {
        ((q_new / q_sat) * 100.0).clamp(0.0, 100.0)
    } else {
        atm.humidity_pct
    };
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn warm_ocean_cool_air() -> (AtmosphericState, OceanState) {
        let atm = AtmosphericState {
            temperature_c: 15.0,
            pressure_hpa: 1013.25,
            humidity_pct: 60.0,
            wind_velocity_ms: [5.0, 3.0, 0.0],
            density_kg_m3: 1.225,
        };
        let ocean = OceanState {
            temperature_c: 25.0,
            salinity_psu: 35.0,
            current_velocity_ms: [0.5, 0.2, 0.0],
            pressure_bar: 1.0,
            density_kg_m3: 1025.0,
        };
        (atm, ocean)
    }

    #[test]
    fn sensible_heat_warm_ocean() {
        let (atm, ocean) = warm_ocean_cool_air();
        let coeff = BulkCoefficients::default();
        let qh = sensible_heat_flux(&atm, &ocean, &coeff);
        // 海が暖かい → 正のフラックス (海→大気)
        assert!(qh > 0.0, "Sensible heat should be positive: {qh:.2}");
    }

    #[test]
    fn latent_heat_warm_ocean() {
        let (atm, ocean) = warm_ocean_cool_air();
        let coeff = BulkCoefficients::default();
        let qe = latent_heat_flux(&atm, &ocean, &coeff);
        // 暖かい海面からの蒸発 → 正
        assert!(qe > 0.0, "Latent heat should be positive: {qe:.2}");
    }

    #[test]
    fn evaporation_positive() {
        let (atm, ocean) = warm_ocean_cool_air();
        let coeff = BulkCoefficients::default();
        let evap = evaporation_rate(&atm, &ocean, &coeff);
        assert!(evap > 0.0, "Evaporation should be positive: {evap:.6}");
    }

    #[test]
    fn net_flux_warm_ocean_positive() {
        let (atm, ocean) = warm_ocean_cool_air();
        let coeff = BulkCoefficients::default();
        let flux = net_heat_flux(&atm, &ocean, &coeff);
        assert!(
            flux > 0.0,
            "Net flux should be positive when ocean is warmer"
        );
    }

    #[test]
    fn net_flux_equilibrium_near_zero() {
        // 大気と海洋が同温 → フラックスは小さい
        let atm = AtmosphericState {
            temperature_c: 20.0,
            pressure_hpa: 1013.25,
            humidity_pct: 100.0, // 飽和
            wind_velocity_ms: [3.0, 0.0, 0.0],
            density_kg_m3: 1.225,
        };
        let ocean = OceanState {
            temperature_c: 20.0,
            salinity_psu: 35.0,
            current_velocity_ms: [0.0; 3],
            pressure_bar: 1.0,
            density_kg_m3: 1025.0,
        };
        let coeff = BulkCoefficients::default();
        let qh = sensible_heat_flux(&atm, &ocean, &coeff);
        assert!(
            qh.abs() < 1.0,
            "Sensible heat at equilibrium should be ~0: {qh:.4}"
        );
    }

    #[test]
    fn coupled_state_fields() {
        let (atm, ocean) = warm_ocean_cool_air();
        let state = compute_coupled_state(atm, ocean);
        assert!(state.net_heat_flux_wm2 > 0.0);
        assert!(state.evaporation_rate > 0.0);
    }

    #[test]
    fn coupled_step_ocean_cools() {
        let (mut atm, mut ocean) = warm_ocean_cool_air();
        let initial_ocean_temp = ocean.temperature_c;
        // 1時間のステップ
        coupled_step(&mut atm, &mut ocean, 3600.0, 50.0);
        // 海が暖かいので海面温度は下がるはず
        assert!(
            ocean.temperature_c < initial_ocean_temp,
            "Ocean should cool: {:.4} vs {:.4}",
            ocean.temperature_c,
            initial_ocean_temp
        );
    }

    #[test]
    fn coupled_step_air_warms() {
        let (mut atm, mut ocean) = warm_ocean_cool_air();
        let initial_air_temp = atm.temperature_c;
        coupled_step(&mut atm, &mut ocean, 3600.0, 50.0);
        assert!(
            atm.temperature_c > initial_air_temp,
            "Air should warm: {:.4} vs {:.4}",
            atm.temperature_c,
            initial_air_temp
        );
    }

    #[test]
    fn coupled_step_humidity_increases() {
        let (mut atm, mut ocean) = warm_ocean_cool_air();
        let initial_humidity = atm.humidity_pct;
        coupled_step(&mut atm, &mut ocean, 3600.0, 50.0);
        assert!(
            atm.humidity_pct >= initial_humidity,
            "Humidity should increase with evaporation: {:.2} vs {:.2}",
            atm.humidity_pct,
            initial_humidity
        );
    }

    #[test]
    fn coupled_step_energy_conservation_direction() {
        // 海の冷却量と大気の加熱量は同方向
        let (mut atm, mut ocean) = warm_ocean_cool_air();
        let t_ocean_0 = ocean.temperature_c;
        let t_air_0 = atm.temperature_c;
        coupled_step(&mut atm, &mut ocean, 3600.0, 50.0);
        let dt_ocean = ocean.temperature_c - t_ocean_0;
        let dt_air = atm.temperature_c - t_air_0;
        // 海は冷え (dt < 0)、大気は温まる (dt > 0)
        assert!(dt_ocean < 0.0 && dt_air > 0.0);
    }

    #[test]
    fn default_bulk_coefficients() {
        let coeff = BulkCoefficients::default();
        assert!(coeff.c_h > 0.0);
        assert!(coeff.c_e > 0.0);
        assert!(coeff.rho_air > 0.0);
        assert!(coeff.cp_air > 0.0);
        assert!(coeff.latent_heat > 0.0);
    }

    #[test]
    fn saturation_humidity_increases_with_temp() {
        let q_cold = saturation_specific_humidity(0.0);
        let q_warm = saturation_specific_humidity(30.0);
        assert!(
            q_warm > q_cold,
            "Warmer air should hold more moisture: {q_warm:.6} vs {q_cold:.6}"
        );
    }
}
