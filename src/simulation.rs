//! Time-stepping シミュレーション — 気候状態の時間発展
//!
//! 大気-海洋結合モデルを時間ステップで進行させる。
//! 各ステップで日射・結合フィードバック・異常検出を実行。

use crate::atmosphere::AtmosphericState;
use crate::coupling::{compute_coupled_state, coupled_step, CoupledState};
use crate::model::{detect_anomaly, AnomalyKind, ClimateAnomaly};
use crate::ocean::OceanState;
use crate::solar::{diurnal_temperature_correction, seasonal_temperature_correction};

/// シミュレーション設定。
#[derive(Debug, Clone, Copy)]
pub struct SimulationConfig {
    /// タイムステップ (秒)。
    pub dt_s: f64,
    /// 混合層深度 (m)。
    pub mixed_layer_depth_m: f64,
    /// 緯度 (度)。
    pub latitude: f64,
    /// 経度 (度)。
    pub longitude: f64,
    /// 日射補正を有効にするか。
    pub enable_solar: bool,
    /// 異常検出を有効にするか。
    pub enable_anomaly_detection: bool,
}

impl Default for SimulationConfig {
    fn default() -> Self {
        Self {
            dt_s: 3600.0, // 1時間
            mixed_layer_depth_m: 50.0,
            latitude: 35.0,
            longitude: 139.0,
            enable_solar: true,
            enable_anomaly_detection: true,
        }
    }
}

/// シミュレーション状態。
#[derive(Debug, Clone)]
pub struct SimulationState {
    /// 現在の大気状態。
    pub atmosphere: AtmosphericState,
    /// 現在の海洋状態。
    pub ocean: OceanState,
    /// 現在時刻 (秒, シミュレーション開始からの経過)。
    pub elapsed_s: f64,
    /// 開始時のタイムスタンプ (ナノ秒, Unix epoch)。
    pub start_timestamp_ns: u64,
    /// 実行済みステップ数。
    pub step_count: u64,
    /// 検出された異常のリスト。
    pub anomalies: Vec<ClimateAnomaly>,
    /// 初期大気状態 (ベースライン)。
    baseline: AtmosphericState,
}

/// ステップ結果。
#[derive(Debug, Clone, Copy)]
pub struct StepResult {
    /// 結合状態。
    pub coupled: CoupledState,
    /// 検出された異常の種類 (あれば)。
    pub anomaly: Option<AnomalyKind>,
    /// 適用された日射温度補正 (°C)。
    pub solar_correction: f64,
}

impl SimulationState {
    /// 初期状態を作成。
    #[must_use]
    pub const fn new(
        atmosphere: AtmosphericState,
        ocean: OceanState,
        start_timestamp_ns: u64,
    ) -> Self {
        Self {
            baseline: atmosphere,
            atmosphere,
            ocean,
            elapsed_s: 0.0,
            start_timestamp_ns,
            step_count: 0,
            anomalies: Vec::new(),
        }
    }

    /// 現在のタイムスタンプ (ナノ秒)。
    #[must_use]
    pub fn current_timestamp_ns(&self) -> u64 {
        self.start_timestamp_ns + (self.elapsed_s * 1e9) as u64
    }

    /// 現在の年間通算日。
    #[must_use]
    pub fn day_of_year(&self) -> u32 {
        let total_seconds = self.current_timestamp_ns() / 1_000_000_000;
        ((total_seconds / 86_400) % 365 + 1) as u32
    }

    /// 現在の地方太陽時 (0.0–24.0)。
    #[must_use]
    pub fn solar_hour(&self, longitude: f64) -> f64 {
        let total_seconds = self.current_timestamp_ns() / 1_000_000_000;
        let utc_hour = (total_seconds % 86_400) as f64 / 3600.0;
        // 経度から地方時へ変換 (15° = 1時間)
        let local_hour = utc_hour + longitude / 15.0;
        ((local_hour % 24.0) + 24.0) % 24.0
    }

    /// 1ステップ進行。
    pub fn step(&mut self, config: &SimulationConfig) -> StepResult {
        // 日射補正
        let mut solar_correction = 0.0;
        if config.enable_solar {
            let day = self.day_of_year();
            let hour = self.solar_hour(config.longitude);

            let diurnal = diurnal_temperature_correction(day, config.latitude, hour);
            let seasonal = seasonal_temperature_correction(day, config.latitude);
            solar_correction = diurnal + seasonal;

            self.atmosphere.temperature_c += solar_correction * (config.dt_s / 3600.0).min(1.0);
        }

        // 大気-海洋結合ステップ
        coupled_step(
            &mut self.atmosphere,
            &mut self.ocean,
            config.dt_s,
            config.mixed_layer_depth_m,
        );

        // 結合状態を計算
        let coupled = compute_coupled_state(self.atmosphere, self.ocean);

        // 異常検出
        let anomaly = if config.enable_anomaly_detection {
            let ts = self.current_timestamp_ns();
            let detected = detect_anomaly(
                &self.atmosphere,
                &self.baseline,
                config.latitude,
                config.longitude,
                ts,
            );
            if let Some(ref a) = detected {
                self.anomalies.push(*a);
            }
            detected.map(|a| a.kind)
        } else {
            None
        };

        self.elapsed_s += config.dt_s;
        self.step_count += 1;

        StepResult {
            coupled,
            anomaly,
            solar_correction,
        }
    }

    /// N ステップ連続実行。最後の結果を返す。
    pub fn run(&mut self, config: &SimulationConfig, steps: u64) -> Option<StepResult> {
        let mut last = None;
        for _ in 0..steps {
            last = Some(self.step(config));
        }
        last
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn default_state() -> SimulationState {
        let atm = AtmosphericState {
            temperature_c: 20.0,
            pressure_hpa: 1013.25,
            humidity_pct: 60.0,
            wind_velocity_ms: [5.0, 2.0, 0.0],
            density_kg_m3: 1.225,
        };
        let ocean = OceanState {
            temperature_c: 22.0,
            salinity_psu: 35.0,
            current_velocity_ms: [0.3, 0.1, 0.0],
            pressure_bar: 1.0,
            density_kg_m3: 1025.0,
        };
        SimulationState::new(atm, ocean, 180u64 * 24 * 3600 * 1_000_000_000)
    }

    #[test]
    fn initial_state() {
        let state = default_state();
        assert_eq!(state.step_count, 0);
        assert!((state.elapsed_s - 0.0).abs() < 1e-10);
        assert!(state.anomalies.is_empty());
    }

    #[test]
    fn single_step() {
        let mut state = default_state();
        let config = SimulationConfig::default();
        let result = state.step(&config);
        assert_eq!(state.step_count, 1);
        assert!((state.elapsed_s - config.dt_s).abs() < 1e-10);
        assert!(result.coupled.net_heat_flux_wm2.is_finite());
    }

    #[test]
    fn multiple_steps() {
        let mut state = default_state();
        let config = SimulationConfig::default();
        let result = state.run(&config, 10);
        assert!(result.is_some());
        assert_eq!(state.step_count, 10);
        assert!(10.0f64.mul_add(-config.dt_s, state.elapsed_s).abs() < 1e-6);
    }

    #[test]
    fn run_zero_steps() {
        let mut state = default_state();
        let config = SimulationConfig::default();
        let result = state.run(&config, 0);
        assert!(result.is_none());
        assert_eq!(state.step_count, 0);
    }

    #[test]
    fn timestamp_advances() {
        let mut state = default_state();
        let config = SimulationConfig::default();
        let ts0 = state.current_timestamp_ns();
        state.step(&config);
        let ts1 = state.current_timestamp_ns();
        assert!(ts1 > ts0);
    }

    #[test]
    fn day_of_year_valid() {
        let state = default_state();
        let day = state.day_of_year();
        assert!((1..=365).contains(&day), "Day should be 1-365, got {day}");
    }

    #[test]
    fn solar_hour_valid() {
        let state = default_state();
        let hour = state.solar_hour(139.0);
        assert!(
            (0.0..24.0).contains(&hour),
            "Solar hour should be 0-24, got {hour:.2}"
        );
    }

    #[test]
    fn solar_correction_applied() {
        let mut state = default_state();
        let config = SimulationConfig {
            enable_solar: true,
            ..Default::default()
        };
        let result = state.step(&config);
        // 日射補正が計算されていることを確認
        assert!(result.solar_correction.is_finite());
    }

    #[test]
    fn no_solar_when_disabled() {
        let mut state = default_state();
        let config = SimulationConfig {
            enable_solar: false,
            ..Default::default()
        };
        let result = state.step(&config);
        assert!(
            (result.solar_correction - 0.0).abs() < 1e-10,
            "Solar correction should be 0 when disabled"
        );
    }

    #[test]
    fn anomaly_detection_works() {
        // 大幅な加温で異常を発生させる
        let atm = AtmosphericState {
            temperature_c: 50.0, // ベースラインより30度高い
            pressure_hpa: 1013.25,
            humidity_pct: 60.0,
            wind_velocity_ms: [5.0, 2.0, 0.0],
            density_kg_m3: 1.225,
        };
        let ocean = OceanState {
            temperature_c: 22.0,
            salinity_psu: 35.0,
            current_velocity_ms: [0.0; 3],
            pressure_bar: 1.0,
            density_kg_m3: 1025.0,
        };
        let mut state = SimulationState::new(atm, ocean, 0);
        // ベースラインを低温に設定
        state.baseline.temperature_c = 20.0;

        let config = SimulationConfig {
            enable_solar: false,
            enable_anomaly_detection: true,
            ..Default::default()
        };
        let result = state.step(&config);
        assert_eq!(result.anomaly, Some(AnomalyKind::HeatWave));
        assert_eq!(state.anomalies.len(), 1);
    }

    #[test]
    fn no_anomaly_normal_conditions() {
        let mut state = default_state();
        let config = SimulationConfig {
            enable_solar: false,
            enable_anomaly_detection: true,
            ..Default::default()
        };
        let result = state.step(&config);
        assert!(result.anomaly.is_none());
    }

    #[test]
    fn temperature_evolves() {
        let mut state = default_state();
        let config = SimulationConfig {
            enable_solar: false,
            ..Default::default()
        };
        let t0 = state.atmosphere.temperature_c;
        state.run(&config, 24); // 24時間
                                // 海洋が暖かいので大気は温まるはず
        assert!(
            (state.atmosphere.temperature_c - t0).abs() > 1e-10,
            "Temperature should have changed"
        );
    }

    #[test]
    fn default_config() {
        let config = SimulationConfig::default();
        assert!((config.dt_s - 3600.0).abs() < 1e-10);
        assert!(config.enable_solar);
        assert!(config.enable_anomaly_detection);
    }
}
