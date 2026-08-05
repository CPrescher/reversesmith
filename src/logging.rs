use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Mutex, OnceLock};
use std::time::{Duration, Instant, SystemTime};

const LOG_FLUSH_INTERVAL: Duration = Duration::from_secs(1);

static LOG_FILE: OnceLock<Mutex<LogWriter<File>>> = OnceLock::new();
static QUIET: AtomicBool = AtomicBool::new(false);

#[derive(Debug)]
struct LogWriter<W: Write> {
    writer: BufWriter<W>,
    last_flush: Instant,
}

impl<W: Write> LogWriter<W> {
    fn new(writer: W, now: Instant) -> Self {
        Self {
            writer: BufWriter::new(writer),
            last_flush: now,
        }
    }

    fn write(&mut self, s: &str, newline: bool, now: Instant) {
        let _ = self.writer.write_all(s.as_bytes());
        if newline {
            let _ = self.writer.write_all(b"\n");
        }
        self.flush_if_due(now);
    }

    fn flush_if_due(&mut self, now: Instant) {
        if now.saturating_duration_since(self.last_flush) >= LOG_FLUSH_INTERVAL {
            self.flush(now);
        }
    }

    fn flush(&mut self, now: Instant) {
        // Advance the deadline even if flushing fails so a broken destination
        // cannot turn every subsequent log write into another flush attempt.
        self.last_flush = now;
        let _ = self.writer.flush();
    }
}

/// Suppress terminal output (log file only).
pub fn set_quiet(quiet: bool) {
    QUIET.store(quiet, Ordering::Relaxed);
}

/// Check whether quiet mode is active.
pub fn is_quiet() -> bool {
    QUIET.load(Ordering::Relaxed)
}

/// Create/overwrite a log file in the given directory.
/// Panics if the file cannot be created (fail fast).
pub fn init_log_file_in(dir: &Path, name: &str) {
    let path = dir.join(name);
    let file = File::create(&path)
        .unwrap_or_else(|e| panic!("Failed to create {}: {}", path.display(), e));
    let writer = LogWriter::new(file, Instant::now());
    LOG_FILE
        .set(Mutex::new(writer))
        .expect("init_log_file called more than once");

    // Write header
    let timestamp = format_utc_timestamp();
    let cwd = std::env::current_dir()
        .map(|p| p.display().to_string())
        .unwrap_or_else(|_| "<unknown>".to_string());
    writeln_to_log(&format!(
        "# rsmith log — {} UTC\n# cwd: {}\n",
        timestamp, cwd
    ));
    flush_log_file();
}

/// Create/overwrite `rsmith.log` in the current working directory.
/// Panics if the file cannot be created (fail fast).
pub fn init_log_file() {
    init_log_file_in(Path::new("."), "rsmith.log");
}

/// Flush the log buffer immediately. Logging also flushes at most once per
/// second while messages are being written.
pub fn flush_log_file() {
    if let Some(mtx) = LOG_FILE.get() {
        if let Ok(mut w) = mtx.lock() {
            w.flush(Instant::now());
        }
    }
}

/// Write a string to the log file (no newline). Errors are silently ignored.
pub fn write_to_log(s: &str) {
    if let Some(mtx) = LOG_FILE.get() {
        if let Ok(mut w) = mtx.lock() {
            w.write(s, false, Instant::now());
        }
    }
}

/// Write a string + newline to the log file. Errors are silently ignored.
pub fn writeln_to_log(s: &str) {
    if let Some(mtx) = LOG_FILE.get() {
        if let Ok(mut w) = mtx.lock() {
            w.write(s, true, Instant::now());
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
        if !$crate::logging::is_quiet() {
            println!();
        }
        $crate::logging::writeln_to_log("");
    };
    ($($arg:tt)*) => {
        {
            let msg = format!($($arg)*);
            if !$crate::logging::is_quiet() {
                println!("{}", msg);
            }
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
            if !$crate::logging::is_quiet() {
                eprintln!("{}", msg);
            }
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
            if !$crate::logging::is_quiet() {
                print!("{}", msg);
            }
            $crate::logging::write_to_log(&msg);
        }
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;
    use std::io;
    use std::rc::Rc;

    #[derive(Default)]
    struct RecordingState {
        bytes: Vec<u8>,
        flushes: usize,
    }

    struct RecordingWriter(Rc<RefCell<RecordingState>>);

    impl Write for RecordingWriter {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            self.0.borrow_mut().bytes.extend_from_slice(buf);
            Ok(buf.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            self.0.borrow_mut().flushes += 1;
            Ok(())
        }
    }

    #[test]
    fn periodic_flush_includes_latest_message_and_is_rate_limited() {
        let state = Rc::new(RefCell::new(RecordingState::default()));
        let start = Instant::now();
        let mut log = LogWriter::new(RecordingWriter(Rc::clone(&state)), start);

        log.write("first", true, start + LOG_FLUSH_INTERVAL / 2);
        assert_eq!(state.borrow().flushes, 0);
        assert!(state.borrow().bytes.is_empty());

        log.write("second", true, start + LOG_FLUSH_INTERVAL);
        assert_eq!(state.borrow().flushes, 1);
        assert_eq!(state.borrow().bytes, b"first\nsecond\n");

        log.write(
            "third",
            true,
            start + LOG_FLUSH_INTERVAL + Duration::from_millis(1),
        );
        assert_eq!(state.borrow().flushes, 1);
    }

    #[test]
    fn explicit_flush_makes_buffered_messages_visible() {
        let state = Rc::new(RefCell::new(RecordingState::default()));
        let start = Instant::now();
        let mut log = LogWriter::new(RecordingWriter(Rc::clone(&state)), start);

        log.write("ready", true, start);
        log.flush(start);

        assert_eq!(state.borrow().flushes, 1);
        assert_eq!(state.borrow().bytes, b"ready\n");
    }
}
