#![no_main]
use kekule::select::{self, SelectParseOptions};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    let options = SelectParseOptions::default()
        .with_max_input_bytes(4096)
        .expect("positive fuzz bounds")
        .with_max_depth(64)
        .expect("positive fuzz bounds")
        .with_host_keywords(["displayed", "selected"])
        .expect("valid host keywords");
    match select::parse(input, options.clone()) {
        Ok(expression) => {
            // Canonical text parses back to an equal expression and is a
            // fixed point of printing.
            let canonical = expression.to_string();
            let reparsed = select::parse(&canonical, options)
                .unwrap_or_else(|error| panic!("canonical {canonical:?} fails: {error}"));
            assert!(reparsed == expression, "canonical {canonical:?} changes meaning");
            assert_eq!(reparsed.to_string(), canonical);
        }
        Err(error) => {
            let span = error.span().expect("selection errors carry a span");
            assert!(span.start <= span.end && span.end <= input.len());
            assert!(input.is_char_boundary(span.start) && input.is_char_boundary(span.end));
        }
    }
});
