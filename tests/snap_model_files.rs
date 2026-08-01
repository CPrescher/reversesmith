use std::path::{Path, PathBuf};

use rsmith::ml_potential::snap_native::SnapModelFiles;

fn test_data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/snap")
        .join(name)
}

#[test]
fn loads_synthetic_multi_element_model_and_mapping() {
    let elements = vec!["O".to_string(), "Si".to_string(), "O".to_string()];
    let model = SnapModelFiles::load(
        &test_data("linear_two_element.snapcoeff"),
        &test_data("linear_two_element.snapparam"),
        &elements,
    )
    .unwrap();

    assert_eq!(model.type_to_element, vec![1, 0, 1]);
    assert_eq!(model.coefficients.elements[0].name, "Si");
    assert_eq!(model.coefficients.ncoeff, 6);
    assert_eq!(model.parameters.twojmax, 2);
    assert_eq!(model.descriptor_count(), 5);
}

#[test]
fn loads_lammps_example_models_when_available() {
    let Some(directory) = std::env::var_os("RSMITH_LAMMPS_POTENTIALS") else {
        eprintln!("skipping LAMMPS SNAP parser examples: RSMITH_LAMMPS_POTENTIALS is not set");
        return;
    };
    let directory = PathBuf::from(directory);

    let silicon = SnapModelFiles::load(
        &directory.join("Si_Zuo_JPCA2020.snapcoeff"),
        &directory.join("Si_Zuo_JPCA2020.snapparam"),
        &["Si".to_string()],
    )
    .unwrap();
    assert_eq!(silicon.coefficients.elements[0].name, "Si");
    assert_eq!(silicon.coefficients.ncoeff, 56);

    let inp = SnapModelFiles::load(
        &directory.join("InP_JCPA2020.snapcoeff"),
        &directory.join("InP_JCPA2020.snapparam"),
        &["In".to_string(), "P".to_string()],
    )
    .unwrap();
    assert!(inp.parameters.chemflag);
    assert_eq!(inp.type_to_element, vec![0, 1]);
}
