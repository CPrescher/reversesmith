use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    if env::var_os("CARGO_FEATURE_GAP_QUIP").is_none() {
        return;
    }

    println!("cargo:rerun-if-changed=c_src/rsmith_gap_quip_lammps_shim.cpp");
    println!("cargo:rerun-if-changed=c_src/rsmith_gap_quip_local_energy_wrapper.F90");
    println!("cargo:rerun-if-changed=include/rsmith_gap_quip_shim.h");

    println!("cargo:rerun-if-env-changed=QUIP_INCLUDE_DIR");
    println!("cargo:rerun-if-env-changed=QUIP_LIB_DIR");
    println!("cargo:rerun-if-env-changed=QUIP_LIBS");
    println!("cargo:rerun-if-env-changed=PKG_CONFIG_PATH");
    println!("cargo:rerun-if-env-changed=PKG_CONFIG_LIBDIR");
    println!("cargo:rerun-if-env-changed=PKG_CONFIG_SYSROOT_DIR");

    if let (Ok(include_dir), Ok(lib_dirs)) =
        (env::var("QUIP_INCLUDE_DIR"), env::var("QUIP_LIB_DIR"))
    {
        link_from_env(&include_dir, &lib_dirs);
        return;
    }

    for package in ["rsmith-gap-quip", "quip-rsmith"] {
        if link_from_pkg_config(package) {
            return;
        }
    }

    panic!(
        "gap-quip feature requires either QUIP_INCLUDE_DIR/QUIP_LIB_DIR/QUIP_LIBS \
         or a pkg-config package named rsmith-gap-quip on PKG_CONFIG_PATH"
    );
}

fn link_from_env(include_dir: &str, lib_dirs: &str) {
    println!("cargo:include={include_dir}");
    for lib_dir in split_paths(lib_dirs) {
        println!("cargo:rustc-link-search=native={}", lib_dir.display());
    }

    let libs = env::var("QUIP_LIBS").unwrap_or_else(|_| "rsmith_quip_gap_shim quip".to_string());
    for lib in split_words(&libs) {
        println!("cargo:rustc-link-lib={lib}");
    }
}

fn link_from_pkg_config(package: &str) -> bool {
    let output = Command::new("pkg-config")
        .args(["--libs", "--cflags", package])
        .output();

    let Ok(output) = output else {
        return false;
    };
    if !output.status.success() {
        return false;
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    let tokens = split_words(&stdout);
    emit_pkg_config_flags(&tokens);
    true
}

fn emit_pkg_config_flags(tokens: &[String]) {
    let mut iter = tokens.iter();
    while let Some(token) = iter.next() {
        if let Some(path) = token.strip_prefix("-I") {
            if path.is_empty() {
                if let Some(path) = iter.next() {
                    println!("cargo:include={path}");
                }
            } else {
                println!("cargo:include={path}");
            }
        } else if let Some(path) = token.strip_prefix("-L") {
            if path.is_empty() {
                if let Some(path) = iter.next() {
                    println!("cargo:rustc-link-search=native={path}");
                }
            } else {
                println!("cargo:rustc-link-search=native={path}");
            }
        } else if let Some(lib) = token.strip_prefix("-l") {
            if lib.is_empty() {
                if let Some(lib) = iter.next() {
                    println!("cargo:rustc-link-lib={lib}");
                }
            } else {
                println!("cargo:rustc-link-lib={lib}");
            }
        } else if token == "-framework" {
            if let Some(framework) = iter.next() {
                println!("cargo:rustc-link-lib=framework={framework}");
            }
        } else if let Some(path) = token.strip_prefix("-F") {
            if path.is_empty() {
                if let Some(path) = iter.next() {
                    println!("cargo:rustc-link-search=framework={path}");
                }
            } else {
                println!("cargo:rustc-link-search=framework={path}");
            }
        } else if token.starts_with("-Wl,") || token == "-pthread" {
            println!("cargo:rustc-link-arg={token}");
        }
    }
}

fn split_paths(value: &str) -> Vec<PathBuf> {
    value
        .split(|c: char| c == ',' || c.is_whitespace())
        .filter(|s| !s.is_empty())
        .flat_map(env::split_paths)
        .collect()
}

fn split_words(value: &str) -> Vec<String> {
    value
        .split(|c: char| c == ',' || c.is_whitespace())
        .filter(|s| !s.is_empty())
        .map(ToOwned::to_owned)
        .collect()
}
