#![no_main]
use kekule::{query, smiles, substructure::*};
use libfuzzer_sys::fuzz_target;
use std::ops::ControlFlow;

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    let Some((pattern, source)) = input.split_once('\n') else {
        return;
    };
    if source.len() > 512 {
        return;
    }
    let Ok(query) = query::parse_smarts_with_options(
        pattern,
        query::SmartsParseOptions {
            max_input_bytes: 1024,
            max_atoms: 32,
            max_recursive_depth: 4,
            ..Default::default()
        },
    ) else {
        return;
    };
    let Ok(molecules) = smiles::to_molecules(source) else {
        return;
    };
    for mut target in molecules {
        if target.perceive().is_err() {
            continue;
        }
        let options = SubstructureMatchOptions {
            max_matches: 64,
            max_search_states: 4096,
            max_candidate_pairs: 8192,
            uniquify: false,
            ..Default::default()
        };
        let collected = find_substructure_matches_complete(&target, &query, options);
        let mut streamed = Vec::new();
        let completion = visit_substructure_matches(&target, &query, options, |m| {
            streamed.push(m.clone());
            ControlFlow::Continue(())
        });
        match (collected, completion) {
            (Ok(matches), Ok(MatchCompletion::Complete)) => assert_eq!(matches, streamed),
            (Err(a), Err(b)) => assert_eq!(a, b),
            other => panic!("inconsistent collection and streaming: {other:?}"),
        }
    }
});
