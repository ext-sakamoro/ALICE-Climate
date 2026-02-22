# Contributing to ALICE-Climate

## Build

```bash
cargo build
```

## Test

```bash
cargo test
```

## Lint

```bash
cargo clippy -- -W clippy::all
cargo fmt -- --check
cargo doc --no-deps 2>&1 | grep warning
```

## Design Constraints

- **Continuous fields**: atmospheric and oceanic states are continuous mathematical functions, not discrete grids.
- **Deterministic**: all field evaluations are pure functions of `(lat, lon, alt, time)` — no hidden state.
- **No external dependencies**: zero crate deps — all physics models are self-contained.
- **Physically grounded**: standard atmosphere model, barometric formula, Coriolis-based wind fields.
