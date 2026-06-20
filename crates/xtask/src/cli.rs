use crate::*;

pub(crate) fn run() -> Result<(), Box<dyn Error>> {
    let mut args = env::args().skip(1);
    match args.next().as_deref() {
        Some("dashboard") => crate::dashboard::dashboard(args.collect()),
        Some("validate") => validate(args.collect()),
        Some("corpus") => crate::corpus::corpus(args.collect()),
        Some("features") => list_features(),
        Some("skills") => skills(args.collect()),
        _ => {
            print_help();
            Ok(())
        }
    }
}

pub(crate) fn print_help() {
    eprintln!(
        "usage:\n  cargo xtask dashboard [--check]\n  cargo xtask validate --feature FEATURE_ID|all [--corpus CORPUS_ID|all] [--update]\n  cargo xtask corpus check --corpus CORPUS_ID|all [--require-data]\n  cargo xtask skills --check\n  cargo xtask features"
    );
}

pub(crate) fn value_after_flag<'a>(args: &'a [String], flag: &str) -> Option<&'a str> {
    args.windows(2)
        .find(|window| window[0] == flag)
        .map(|window| window[1].as_str())
}

pub(crate) fn validate_args(args: &[String]) -> Result<(), Box<dyn Error>> {
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "--feature" | "--corpus" => {
                if index + 1 >= args.len() {
                    return Err(boxed_error(format!("missing value after {}", args[index])));
                }
                index += 2;
            }
            "--update" => index += 1,
            arg => return Err(boxed_error(format!("unknown validate argument: {arg}"))),
        }
    }
    Ok(())
}
