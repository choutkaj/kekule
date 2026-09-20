//! Refresh the optional local report view without changing the scientific exit status.
use super::*;

pub(super) fn publish(opts: &Options) {
    if opts.generate {
        return;
    }
    if let Err(error) = refresh(opts) {
        eprintln!(
            "warning: dashboard was not refreshed: {error}\nReport remains at {}",
            opts.output.display()
        );
    }
}

fn refresh(opts: &Options) -> Result<(), Box<dyn Error>> {
    fs::create_dir_all(&opts.runs_dir)?;
    // Explicit --output still owns its case records. Archive only its small summary.
    if opts
        .output
        .extension()
        .is_none_or(|extension| extension != "json")
        || opts.output.canonicalize()?.parent() != Some(opts.runs_dir.canonicalize()?.as_path())
    {
        let archive = opts.runs_dir.join(format!(
            "run-{}-{}.json",
            opts.started_at_unix_ms,
            process::id()
        ));
        ReportFile::new(&archive).write(|file| {
            std::io::copy(&mut fs::File::open(&opts.output)?, file)?;
            Ok(())
        })?;
    }
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let configured = env::var_os("KEKULE_DASHBOARD_PYTHON")
        .map(PathBuf::from)
        .or_else(|| {
            fs::read_to_string(root.join(".dashboard-python"))
                .ok()
                .map(|value| PathBuf::from(value.trim()))
                .filter(|path| !path.as_os_str().is_empty())
        });
    let mut candidates = Vec::new();
    if let Some(path) = configured {
        candidates.push((path, false));
    } else {
        if let Some(path) = opts
            .python
            .clone()
            .or_else(|| env::var_os("KEKULE_WRITER_PYTHON").map(PathBuf::from))
        {
            candidates.push((path, false));
        }
        candidates.extend([("python3".into(), false), ("python".into(), false)]);
        if cfg!(windows) {
            candidates.push(("py".into(), true));
        }
    }
    for (python, launcher) in candidates {
        let mut command = process::Command::new(&python);
        if launcher {
            command.arg("-3");
        }
        let available = command
            .args(["-c", "import sys; sys.exit(sys.version_info < (3, 11))"])
            .output()
            .is_ok_and(|output| output.status.success());
        if !available {
            continue;
        }
        let mut command = process::Command::new(python);
        if launcher {
            command.arg("-3");
        }
        let output = command
            .arg(root.join("dashboard.py"))
            .arg("--runs-dir")
            .arg(&opts.runs_dir)
            .output()?;
        if !output.status.success() {
            return Err(boxed_error(String::from_utf8_lossy(&output.stderr).trim()));
        }
        if !output.stderr.is_empty() {
            eprint!("{}", String::from_utf8_lossy(&output.stderr));
        }
        println!("Dashboard: {}", opts.runs_dir.join("index.html").display());
        return Ok(());
    }
    Err(boxed_error("Python 3.11+ is needed; set KEKULE_DASHBOARD_PYTHON or put its executable path in benchmarks/.dashboard-python"))
}
