use std::path::Path;
use std::sync::atomic::{AtomicU32, Ordering};

use rsmith::config::Config;

static COUNTER: AtomicU32 = AtomicU32::new(0);

/// Helper: write a temp TOML file with unique name, try to load it, return the error string (if any).
fn load_toml(content: &str) -> Result<Config, String> {
    let id = COUNTER.fetch_add(1, Ordering::Relaxed);
    let dir = std::env::temp_dir().join("rsmith_test_config");
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join(format!("test_{}.toml", id));
    std::fs::write(&path, content).unwrap();
    let result = Config::load(Path::new(&path)).map_err(|e| e.to_string());
    let _ = std::fs::remove_file(&path);
    result
}

/// Minimal valid config for baseline.
const MINIMAL_CONFIG: &str = r#"
[system]
structure = "test.xyz"
format = "xyz"

[data]

[rmc]
"#;

#[test]
fn minimal_config_is_valid() {
    let result = load_toml(MINIMAL_CONFIG);
    assert!(
        result.is_ok(),
        "Minimal config should parse: {:?}",
        result.err()
    );
}

#[test]
fn unknown_system_field() {
    let toml = r#"
[system]
structure = "test.xyz"
format = "xyz"
typo = true

[data]
[rmc]
"#;
    let err = load_toml(toml).unwrap_err();
    assert!(
        err.contains("unknown field"),
        "Expected unknown field error, got: {err}"
    );
    assert!(
        err.contains("typo"),
        "Error should mention the typo field, got: {err}"
    );
}

#[test]
fn unknown_rmc_field() {
    let toml = r#"
[system]
structure = "test.xyz"
format = "xyz"

[data]

[rmc]
max_mooves = 100
"#;
    let err = load_toml(toml).unwrap_err();
    assert!(
        err.contains("unknown field"),
        "Expected unknown field error, got: {err}"
    );
    assert!(
        err.contains("max_mooves"),
        "Error should mention the typo, got: {err}"
    );
}

#[test]
fn unknown_epsr_field() {
    let toml = r#"
[system]
structure = "test.xyz"
format = "xyz"

[data]

[rmc]

[epsr]
feedbak = 0.5
"#;
    let err = load_toml(toml).unwrap_err();
    assert!(
        err.contains("unknown field"),
        "Expected unknown field error, got: {err}"
    );
    assert!(
        err.contains("feedbak"),
        "Error should mention the typo, got: {err}"
    );
}

#[test]
fn unknown_potential_field() {
    let toml = r#"
[system]
structure = "test.xyz"
format = "xyz"

[data]

[rmc]

[potential]
weght = 2.0
"#;
    let err = load_toml(toml).unwrap_err();
    assert!(
        err.contains("unknown field"),
        "Expected unknown field error, got: {err}"
    );
    assert!(
        err.contains("weght"),
        "Error should mention the typo, got: {err}"
    );
}

#[test]
fn unknown_sq_field() {
    let toml = r#"
[system]
structure = "test.xyz"
format = "xyz"

[data]

[rmc]

[sq]
qmaxx = 20.0
"#;
    let err = load_toml(toml).unwrap_err();
    assert!(
        err.contains("unknown field"),
        "Expected unknown field error, got: {err}"
    );
    assert!(
        err.contains("qmaxx"),
        "Error should mention the typo, got: {err}"
    );
}

#[test]
fn unknown_constraints_field() {
    let toml = r#"
[system]
structure = "test.xyz"
format = "xyz"

[data]

[rmc]

[constraints]
min_distanc = {}
"#;
    let err = load_toml(toml).unwrap_err();
    assert!(
        err.contains("unknown field"),
        "Expected unknown field error, got: {err}"
    );
    assert!(
        err.contains("min_distanc"),
        "Error should mention the typo, got: {err}"
    );
}

#[test]
fn unknown_dataset_field() {
    let toml = r#"
[system]
structure = "test.xyz"
format = "xyz"

[data]
[data.xray_sq]
file = "test.sq"
weigth = 1.0

[rmc]
"#;
    let err = load_toml(toml).unwrap_err();
    assert!(
        err.contains("unknown field"),
        "Expected unknown field error, got: {err}"
    );
    assert!(
        err.contains("weigth"),
        "Error should mention the typo, got: {err}"
    );
}

#[test]
fn unknown_top_level_section() {
    let toml = r#"
[system]
structure = "test.xyz"
format = "xyz"

[data]

[rmc]

[rmcc]
max_moves = 100
"#;
    let err = load_toml(toml).unwrap_err();
    assert!(
        err.contains("unknown field"),
        "Expected unknown field error, got: {err}"
    );
    assert!(
        err.contains("rmcc"),
        "Error should mention the typo, got: {err}"
    );
}

#[test]
fn valid_epsr_fields_accepted() {
    let toml = r#"
[system]
structure = "test.xyz"
format = "xyz"

[data]

[rmc]

[epsr]
iterations = 10
feedback = 0.2
smooth_sigma = 0.02
moves_per_iteration = 200000
temperature = 0.025
min_r = 1.0
convergence = 0.0001
ep_restart = "../prev"
"#;
    let result = load_toml(toml);
    assert!(
        result.is_ok(),
        "Valid EPSR config should parse: {:?}",
        result.err()
    );
}

#[test]
fn valid_pedone_fields_accepted() {
    let toml = r#"
[system]
structure = "test.xyz"
format = "xyz"

[data]

[rmc]

[potential]
weight = 2.0
cutoff = 10.0

[[potential.pedone]]
pair = "Si-O"
D0 = 0.34
alpha = 2.0
r0 = 2.1
C0 = 1.0
"#;
    let result = load_toml(toml);
    assert!(
        result.is_ok(),
        "Valid Pedone config should parse: {:?}",
        result.err()
    );
}
