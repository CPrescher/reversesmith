use std::fs;
use std::thread;
use std::time::{Duration, SystemTime};

#[test]
fn quiet_mode_still_periodically_flushes_to_disk() {
    let unique = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("rsmith-logging-{}-{unique}", std::process::id()));
    fs::create_dir(&dir).unwrap();

    rsmith::logging::set_quiet(true);
    rsmith::logging::init_log_file_in(&dir, "quiet.log");
    rsmith::log_println!("quiet message before interval");
    thread::sleep(Duration::from_millis(1020));
    rsmith::log_println!("quiet progress message");

    // Read before an explicit flush or process shutdown: the second write
    // must have made both buffered messages visible.
    let contents = fs::read_to_string(dir.join("quiet.log")).unwrap();

    assert!(contents.contains("quiet message before interval"));
    assert!(contents.contains("quiet progress message"));

    // The process-wide logger intentionally keeps the file open. Cleanup is
    // best-effort for platforms that permit unlinking it.
    let _ = fs::remove_file(dir.join("quiet.log"));
    let _ = fs::remove_dir(dir);
}
