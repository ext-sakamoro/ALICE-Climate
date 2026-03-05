//! 季節・日周補正 — 太陽位置と日射モデル
//!
//! 太陽赤緯・時角・天頂角から日射量を計算し、
//! 季節変動と日周変動を大気モデルに反映する。

use core::f64::consts::PI;

/// 太陽定数 (W/m²)。
pub const SOLAR_CONSTANT: f64 = 1361.0;

/// 太陽赤緯 (ラジアン)。
///
/// Cooper の近似式: δ = 23.45° × sin(360/365 × (284 + `day_of_year`))
#[must_use]
pub fn solar_declination(day_of_year: u32) -> f64 {
    let angle = 2.0 * PI * f64::from(284 + day_of_year) / 365.0;
    (23.45_f64).to_radians() * angle.sin()
}

/// 時角 (ラジアン)。
///
/// `hour` は地方太陽時 (0.0–24.0)。正午 = 0。
#[must_use]
pub fn solar_hour_angle(solar_hour: f64) -> f64 {
    (solar_hour - 12.0) * (PI / 12.0)
}

/// 太陽天頂角の余弦。
///
/// cos(θz) = sin(φ)sin(δ) + cos(φ)cos(δ)cos(ω)
///
/// - `latitude_rad` — 緯度 (ラジアン)
/// - `declination_rad` — 太陽赤緯 (ラジアン)
/// - `hour_angle_rad` — 時角 (ラジアン)
#[must_use]
pub fn cos_solar_zenith(latitude_rad: f64, declination_rad: f64, hour_angle_rad: f64) -> f64 {
    let sin_lat = latitude_rad.sin();
    let cos_lat = latitude_rad.cos();
    let sin_dec = declination_rad.sin();
    let cos_dec = declination_rad.cos();
    let cos_h = hour_angle_rad.cos();
    sin_lat.mul_add(sin_dec, cos_lat * cos_dec * cos_h)
}

/// 太陽天頂角 (ラジアン)。
#[must_use]
pub fn solar_zenith_angle(latitude_rad: f64, declination_rad: f64, hour_angle_rad: f64) -> f64 {
    let cos_z = cos_solar_zenith(latitude_rad, declination_rad, hour_angle_rad);
    cos_z.clamp(-1.0, 1.0).acos()
}

/// 大気上端日射量 (W/m²)。
///
/// 天頂角が 90° を超える場合 (夜間) は 0。
#[must_use]
pub fn insolation(latitude_rad: f64, declination_rad: f64, hour_angle_rad: f64) -> f64 {
    let cos_z = cos_solar_zenith(latitude_rad, declination_rad, hour_angle_rad);
    if cos_z <= 0.0 {
        return 0.0;
    }
    SOLAR_CONSTANT * cos_z
}

/// 地表到達日射量 (W/m²)。
///
/// 大気透過率 (Beer-Lambert 近似) を適用。
/// `altitude_m` が高いほど大気層が薄く透過率が高い。
#[must_use]
pub fn surface_insolation(
    latitude_rad: f64,
    declination_rad: f64,
    hour_angle_rad: f64,
    altitude_m: f64,
) -> f64 {
    let cos_z = cos_solar_zenith(latitude_rad, declination_rad, hour_angle_rad);
    if cos_z <= 0.0 {
        return 0.0;
    }

    // エアマス (簡易近似: 1/cos(θz))、最大 38 で制限
    let air_mass = (1.0 / cos_z).min(38.0);

    // 大気透過率: 高度補正 (scale height ≈ 8500m)
    let scale_height = 8500.0;
    let atm_factor = (-altitude_m / scale_height).exp();
    let transmittance = (0.7_f64).powf(air_mass * atm_factor);

    SOLAR_CONSTANT * cos_z * transmittance
}

/// 日周温度補正 (°C)。
///
/// 日中の加熱と夜間の放射冷却を近似。
/// 最高気温は太陽正午の約2時間後 (14時)。
#[must_use]
pub fn diurnal_temperature_correction(day_of_year: u32, latitude: f64, solar_hour: f64) -> f64 {
    let declination = solar_declination(day_of_year);
    let lat_rad = latitude.to_radians();
    let hour_angle = solar_hour_angle(solar_hour);

    let cos_z = cos_solar_zenith(lat_rad, declination, hour_angle);

    // 日照がある場合: 最高気温は14時頃
    // 位相シフト: (hour - 14) を使った余弦関数
    let phase_shift = (solar_hour - 14.0) * (PI / 12.0);
    let diurnal_wave = phase_shift.cos();

    // 振幅は日射量と緯度に依存 (赤道付近で大きい)
    let amplitude = 5.0 * lat_rad.cos().abs();

    if cos_z > 0.0 {
        // 日中: 正の補正
        amplitude * diurnal_wave
    } else {
        // 夜間: 放射冷却 (負の補正)
        -amplitude * 0.5 * (-diurnal_wave).max(0.0)
    }
}

/// 季節温度補正 (°C)。
///
/// 夏至/冬至を基準に緯度に応じた季節変動を返す。
#[must_use]
pub fn seasonal_temperature_correction(day_of_year: u32, latitude: f64) -> f64 {
    let declination = solar_declination(day_of_year);
    let lat_rad = latitude.to_radians();

    // 季節効果: 太陽赤緯と緯度の関係
    // 北半球夏 (declination > 0, lat > 0) → 正の補正
    let seasonal_factor = lat_rad.sin() * declination.sin();

    // 振幅: 高緯度ほど大きい (最大 ±15°C)
    15.0 * seasonal_factor / (23.45_f64).to_radians().sin()
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn declination_summer_solstice() {
        // Day ~172 (June 21): 赤緯 ≈ +23.45°
        let dec = solar_declination(172);
        let deg = dec.to_degrees();
        assert!(
            (deg - 23.45).abs() < 1.0,
            "Summer solstice declination should be ~23.45°, got {deg:.2}°"
        );
    }

    #[test]
    fn declination_winter_solstice() {
        // Day ~355 (Dec 21): 赤緯 ≈ -23.45°
        let dec = solar_declination(355);
        let deg = dec.to_degrees();
        assert!(
            (deg + 23.45).abs() < 1.5,
            "Winter solstice declination should be ~-23.45°, got {deg:.2}°"
        );
    }

    #[test]
    fn declination_equinox() {
        // Day ~80 (Mar 21): 赤緯 ≈ 0°
        let dec = solar_declination(80);
        let deg = dec.to_degrees();
        assert!(
            deg.abs() < 2.0,
            "Equinox declination should be ~0°, got {deg:.2}°"
        );
    }

    #[test]
    fn hour_angle_noon() {
        let ha = solar_hour_angle(12.0);
        assert!(ha.abs() < 1e-10, "Noon hour angle should be 0");
    }

    #[test]
    fn hour_angle_morning_evening() {
        let morning = solar_hour_angle(6.0);
        let evening = solar_hour_angle(18.0);
        assert!(morning < 0.0, "Morning hour angle should be negative");
        assert!(evening > 0.0, "Evening hour angle should be positive");
        assert!((morning + evening).abs() < 1e-10, "Should be symmetric");
    }

    #[test]
    fn zenith_noon_equator_equinox() {
        // 春分の日の赤道正午: 天頂角 ≈ 0°
        let dec = solar_declination(80);
        let ha = solar_hour_angle(12.0);
        let zenith = solar_zenith_angle(0.0, dec, ha);
        assert!(
            zenith.to_degrees() < 5.0,
            "Equator noon equinox zenith should be ~0°, got {:.2}°",
            zenith.to_degrees()
        );
    }

    #[test]
    fn insolation_noon_vs_night() {
        let dec = solar_declination(172);
        let lat = (45.0_f64).to_radians();
        let noon = insolation(lat, dec, solar_hour_angle(12.0));
        let midnight = insolation(lat, dec, solar_hour_angle(0.0));
        assert!(noon > 0.0, "Noon insolation should be positive");
        assert!(midnight.abs() < 1e-10, "Midnight insolation should be 0");
    }

    #[test]
    fn insolation_equator_max() {
        // 赤道の正午日射量は太陽定数に近い
        let dec = solar_declination(80); // equinox
        let noon = insolation(0.0, dec, solar_hour_angle(12.0));
        assert!(
            noon > SOLAR_CONSTANT * 0.95,
            "Equator equinox noon should be ~{SOLAR_CONSTANT}, got {noon:.1}"
        );
    }

    #[test]
    fn surface_insolation_less_than_toa() {
        let dec = solar_declination(172);
        let lat = (35.0_f64).to_radians();
        let ha = solar_hour_angle(12.0);
        let toa = insolation(lat, dec, ha);
        let surface = surface_insolation(lat, dec, ha, 0.0);
        assert!(
            surface < toa,
            "Surface ({surface:.1}) should be less than TOA ({toa:.1})"
        );
        assert!(surface > 0.0);
    }

    #[test]
    fn surface_insolation_altitude_increases() {
        let dec = solar_declination(172);
        let lat = (35.0_f64).to_radians();
        let ha = solar_hour_angle(12.0);
        let sea_level = surface_insolation(lat, dec, ha, 0.0);
        let mountain = surface_insolation(lat, dec, ha, 5000.0);
        assert!(
            mountain > sea_level,
            "Higher altitude ({mountain:.1}) should have more insolation than sea level ({sea_level:.1})"
        );
    }

    #[test]
    fn diurnal_afternoon_warmer() {
        // 14時付近が最も暖かい
        let correction_14 = diurnal_temperature_correction(172, 35.0, 14.0);
        let correction_06 = diurnal_temperature_correction(172, 35.0, 6.0);
        assert!(
            correction_14 > correction_06,
            "14:00 ({correction_14:.2}) should be warmer than 06:00 ({correction_06:.2})"
        );
    }

    #[test]
    fn seasonal_summer_warmer_north() {
        // 北半球夏至: 正の補正
        let summer = seasonal_temperature_correction(172, 45.0);
        let winter = seasonal_temperature_correction(355, 45.0);
        assert!(summer > 0.0, "NH summer should be positive: {summer:.2}");
        assert!(winter < 0.0, "NH winter should be negative: {winter:.2}");
    }

    #[test]
    fn seasonal_opposite_hemispheres() {
        let nh_summer = seasonal_temperature_correction(172, 45.0);
        let sh_summer = seasonal_temperature_correction(172, -45.0);
        // 反対の符号
        assert!(
            nh_summer * sh_summer < 0.0,
            "Opposite hemispheres should have opposite signs"
        );
    }

    #[test]
    fn seasonal_equator_small() {
        // 赤道では季節変動が小さい
        let correction = seasonal_temperature_correction(172, 0.0);
        assert!(
            correction.abs() < 1.0,
            "Equator seasonal correction should be small: {correction:.2}"
        );
    }
}
