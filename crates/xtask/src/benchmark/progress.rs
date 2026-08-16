use std::io::{self, IsTerminal, Write};
use std::sync::{Arc, Mutex};

const PROGRESS_BAR_WIDTH: usize = 24;

pub(crate) fn progress_bar(completed: usize, total: usize) -> String {
    let clamped = completed.min(total);
    let filled = (clamped * PROGRESS_BAR_WIDTH)
        .checked_div(total)
        .unwrap_or(PROGRESS_BAR_WIDTH);
    let percent = (clamped * 100).checked_div(total).unwrap_or(100);
    format!(
        "[{}{}] {}/{} {:>3}%",
        "#".repeat(filled),
        "-".repeat(PROGRESS_BAR_WIDTH - filled),
        clamped,
        total,
        percent
    )
}

pub(crate) fn benchmark_worker_count(jobs: usize, fixture_count: usize) -> usize {
    jobs.max(1).min(fixture_count.max(1))
}

pub(crate) struct BenchmarkProgress {
    total_targets: usize,
    completed_targets: usize,
    matched_targets: usize,
    failed_targets: usize,
}

impl BenchmarkProgress {
    pub(crate) fn start(total_targets: usize, jobs: usize) -> Self {
        println!("benchmark plan: {total_targets} target(s), {jobs} worker(s)");
        println!("overall {}", progress_bar(0, total_targets));
        Self {
            total_targets,
            completed_targets: 0,
            matched_targets: 0,
            failed_targets: 0,
        }
    }

    pub(crate) fn target_start(&self, index: usize, benchmark_id: &str, corpus: &str) {
        println!();
        println!("[{index}/{}] {benchmark_id} [{corpus}]", self.total_targets);
    }

    pub(crate) fn manifest(&self, reference_tool: &str, reference_version: &str) {
        if reference_version
            .to_ascii_lowercase()
            .starts_with(&reference_tool.to_ascii_lowercase())
        {
            println!("  manifest: {reference_version}");
        } else {
            println!("  manifest: {reference_tool} {reference_version}");
        }
    }

    pub(crate) fn target_match(&mut self, match_count: usize, fixture_count: usize) {
        println!("  result: match ({match_count}/{fixture_count} fixture(s))");
        self.completed_targets += 1;
        self.matched_targets += 1;
        self.print_overall();
    }

    pub(crate) fn target_differences(
        &mut self,
        difference_count: usize,
        match_count: usize,
        fixture_count: usize,
    ) {
        println!(
            "  result: differences ({difference_count} different, {match_count}/{fixture_count} matched)"
        );
        self.completed_targets += 1;
        self.failed_targets += 1;
        self.print_overall();
    }

    pub(crate) fn target_error(&mut self) {
        println!("  result: error");
        self.completed_targets += 1;
        self.failed_targets += 1;
        self.print_overall();
    }

    fn print_overall(&self) {
        println!(
            "overall {} matched={} differences-or-errors={}",
            progress_bar(self.completed_targets, self.total_targets),
            self.matched_targets,
            self.failed_targets
        );
    }
}

#[derive(Clone)]
pub(crate) struct FixtureProgress {
    inner: Arc<Mutex<FixtureProgressState>>,
}

struct FixtureProgressState {
    total: usize,
    completed: usize,
    worker_count: usize,
    interactive: bool,
    last_filled: usize,
}

impl FixtureProgress {
    pub(crate) fn start(total: usize, worker_count: usize) -> Self {
        let progress = Self {
            inner: Arc::new(Mutex::new(FixtureProgressState {
                total,
                completed: 0,
                worker_count,
                interactive: io::stdout().is_terminal(),
                last_filled: 0,
            })),
        };
        progress.print_current(true);
        progress
    }

    pub(crate) fn fixture_finished(&self) {
        let mut state = self
            .inner
            .lock()
            .expect("benchmark fixture progress lock should not be poisoned");
        state.completed = (state.completed + 1).min(state.total);
        let filled = progress_filled_segments(state.completed, state.total);
        if filled != state.last_filled || state.completed == 1 || state.completed == state.total {
            state.last_filled = filled;
            print_fixture_progress(&state, false);
        }
    }

    pub(crate) fn finish(&self) {
        let state = self
            .inner
            .lock()
            .expect("benchmark fixture progress lock should not be poisoned");
        if state.interactive {
            println!();
        }
    }

    fn print_current(&self, start: bool) {
        let state = self
            .inner
            .lock()
            .expect("benchmark fixture progress lock should not be poisoned");
        print_fixture_progress(&state, start);
    }
}

fn progress_filled_segments(completed: usize, total: usize) -> usize {
    (completed.min(total) * PROGRESS_BAR_WIDTH)
        .checked_div(total)
        .unwrap_or(PROGRESS_BAR_WIDTH)
}

fn print_fixture_progress(state: &FixtureProgressState, start: bool) {
    let line = format!(
        "  fixtures {} ({} worker(s))",
        progress_bar(state.completed, state.total),
        state.worker_count
    );
    if state.interactive {
        print!("\r{line}");
        io::stdout()
            .flush()
            .expect("benchmark progress should flush");
    } else if start || state.completed > 0 {
        println!("{line}");
    }
}
