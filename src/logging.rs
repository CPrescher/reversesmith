use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::{Mutex, OnceLock};
use std::time::SystemTime;

static LOG_FILE: OnceLock<Mutex<BufWriter<File>>> = OnceLock::new();

/// Create/overwrite `reversesmith.log` in the current working directory.
/// Panics if the file cannot be created (fail fast).
pub fn init_log_file() {
    let file = File::create("reversesmith.log").expect("Failed to create reversesmith.log");
    let writer = BufWriter::new(file);
    LOG_FILE
        .set(Mutex::new(writer))
        .expect("init_log_file called more than once");

    // Write header
    let timestamp = format_utc_timestamp();
    let cwd = std::env::current_dir()
        .map(|p| p.display().to_string())
        .unwrap_or_else(|_| "<unknown>".to_string());
    writeln_to_log(&format!(
        "# reversesmith log — {} UTC\n# cwd: {}\n",
        timestamp, cwd
    ));
}

/// Flush the log buffer. Call at the end of main().
pub fn flush_log_file() {
    if let Some(mtx) = LOG_FILE.get() {
        if let Ok(mut w) = mtx.lock() {
            let _ = w.flush();
        }
    }
}

/// Write a string to the log file (no newline). Errors are silently ignored.
pub fn write_to_log(s: &str) {
    if let Some(mtx) = LOG_FILE.get() {
        if let Ok(mut w) = mtx.lock() {
            let _ = w.write_all(s.as_bytes());
        }
    }
}

/// Write a string + newline to the log file. Errors are silently ignored.
pub fn writeln_to_log(s: &str) {
    if let Some(mtx) = LOG_FILE.get() {
        if let Ok(mut w) = mtx.lock() {
            let _ = w.write_all(s.as_bytes());
            let _ = w.write_all(b"\n");
        }
    }
}

/// Format current UTC time as "YYYY-MM-DD HH:MM:SS" using only std.
fn format_utc_timestamp() -> String {
    let dur = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .unwrap_or_default();
    let secs = dur.as_secs();

    // Days since epoch
    let days = secs / 86400;
    let time_of_day = secs % 86400;
    let h = time_of_day / 3600;
    let m = (time_of_day % 3600) / 60;
    let s = time_of_day % 60;

    // Civil date from day count (algorithm from Howard Hinnant)
    let z = days as i64 + 719468;
    let era = z.div_euclid(146097);
    let doe = z.rem_euclid(146097) as u64; // day of era [0, 146096]
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = if month <= 2 { y + 1 } else { y };

    format!(
        "{:04}-{:02}-{:02} {:02}:{:02}:{:02}",
        year, month, d, h, m, s
    )
}

/// Like `println!` but also writes to the log file.
#[macro_export]
macro_rules! log_println {
    () => {
        println!();
        $crate::logging::writeln_to_log("");
    };
    ($($arg:tt)*) => {
        {
            let msg = format!($($arg)*);
            println!("{}", msg);
            $crate::logging::writeln_to_log(&msg);
        }
    };
}

/// Like `eprintln!` but also writes to the log file.
#[macro_export]
macro_rules! log_eprintln {
    ($($arg:tt)*) => {
        {
            let msg = format!($($arg)*);
            eprintln!("{}", msg);
            $crate::logging::writeln_to_log(&msg);
        }
    };
}

/// Like `print!` but also writes to the log file (no newline).
#[macro_export]
macro_rules! log_print {
    ($($arg:tt)*) => {
        {
            let msg = format!($($arg)*);
            print!("{}", msg);
            $crate::logging::write_to_log(&msg);
        }
    };
}
