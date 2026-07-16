use crate::*;

pub(crate) fn dashboard(args: Vec<String>) -> Result<(), Box<dyn Error>> {
    let check = args.iter().any(|arg| arg == "--check");
    let features = read_features()?;
    let statuses = read_validation_statuses(&features)?;
    let corpus_info = read_dashboard_corpus_info()?;
    let rendered = render_dashboard(&features, &statuses, &corpus_info);
    let path = Path::new(DASHBOARD_PATH);

    if check {
        let existing = fs::read_to_string(path)?;
        if normalize_text_line_endings(&existing) != normalize_text_line_endings(&rendered) {
            return Err(boxed_error(
                "features/DASHBOARD.html is out of date; run `cargo xtask dashboard`",
            ));
        }
    } else {
        write_atomic_text(path, &rendered)?;
    }
    Ok(())
}

pub(crate) fn render_dashboard(
    features: &[Feature],
    statuses: &BTreeMap<String, ValidationStatus>,
    corpus_info: &BTreeMap<String, CorpusDashboardInfo>,
) -> String {
    let mut out = String::new();
    out.push_str("<!doctype html>\n");
    out.push_str("<html lang=\"en\">\n");
    out.push_str("<head>\n");
    out.push_str("<meta charset=\"utf-8\">\n");
    out.push_str("<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n");
    out.push_str("<title>Feature Dashboard</title>\n");
    out.push_str("<style>\n");
    out.push_str(":root { color-scheme: light dark; --border: #d0d7de; --head: #f6f8fa; --ok: #1a7f37; --bad: #cf222e; --unknown: #9a6700; --muted: #656d76; --text: #24292f; --bg: #ffffff; }\n");
    out.push_str("@media (prefers-color-scheme: dark) { :root { --border: #30363d; --head: #161b22; --ok: #3fb950; --bad: #ff7b72; --unknown: #d29922; --muted: #8b949e; --text: #c9d1d9; --bg: #0d1117; } }\n");
    out.push_str("body { margin: 24px; background: var(--bg); color: var(--text); font: 14px/1.4 system-ui, -apple-system, Segoe UI, sans-serif; }\n");
    out.push_str("h1 { margin: 0 0 4px; font-size: 24px; }\n");
    out.push_str("p { margin: 0 0 18px; color: var(--muted); }\n");
    out.push_str(".dashboard-wrap { overflow-x: auto; }\n");
    out.push_str("table { border-collapse: collapse; width: 100%; min-width: 980px; }\n");
    out.push_str(
        "th, td { border: 1px solid var(--border); padding: 6px 8px; vertical-align: middle; }\n",
    );
    out.push_str("thead th { position: sticky; top: 0; z-index: 1; height: 168px; background: var(--head); white-space: nowrap; }\n");
    out.push_str("tbody tr:nth-child(even) { background: color-mix(in srgb, var(--head) 45%, transparent); }\n");
    out.push_str("th.text, td.text { text-align: left; }\n");
    out.push_str("th.compact, td.compact, th.rotated, td.marker { text-align: center; }\n");
    out.push_str("th.area, td.area { text-align: left; }\n");
    out.push_str(
        "th.rotated { width: 52px; min-width: 52px; padding: 0; vertical-align: bottom; overflow: hidden; }\n",
    );
    out.push_str("th.rotated button { position: relative; height: 168px; width: 52px; padding: 0; display: block; overflow: hidden; }\n");
    out.push_str("th.rotated .rotated-label { position: absolute; left: calc(50% + 23px); bottom: 12px; width: 144px; height: 46px; display: flex; flex-direction: column; justify-content: center; align-items: flex-start; transform: rotate(-90deg); transform-origin: left bottom; text-align: left; line-height: 1.15; }\n");
    out.push_str("th.rotated .rotated-name, th.rotated .rotated-count { white-space: nowrap; }\n");
    out.push_str(
        "th.rotated .rotated-count { font-size: 12px; font-weight: 650; color: var(--muted); }\n",
    );
    out.push_str(
        "button.sort { all: unset; cursor: pointer; color: inherit; font-weight: 650; }\n",
    );
    out.push_str(
        "button.sort:focus-visible { outline: 2px solid Highlight; outline-offset: 2px; }\n",
    );
    out.push_str("th[aria-sort=\"ascending\"] button.sort::after { content: \" \\25B2\"; font-size: 10px; color: var(--muted); }\n");
    out.push_str("th[aria-sort=\"descending\"] button.sort::after { content: \" \\25BC\"; font-size: 10px; color: var(--muted); }\n");
    out.push_str(".ok { color: var(--ok); font-weight: 700; }\n");
    out.push_str(".bad { color: var(--bad); font-weight: 700; }\n");
    out.push_str(".unknown { color: var(--unknown); font-weight: 700; }\n");
    out.push_str(".count { display: inline-block; min-width: 1.2em; margin-left: 2px; font-size: 12px; font-weight: 650; color: var(--bad); }\n");
    out.push_str(".na { color: var(--muted); }\n");
    out.push_str(".legend { margin-top: -8px; font-size: 12px; }\n");
    out.push_str(".legend span { margin-right: 3px; }\n");
    out.push_str("code { font-family: ui-monospace, SFMono-Regular, Consolas, monospace; font-size: 13px; }\n");
    out.push_str("</style>\n");
    out.push_str("</head>\n");
    out.push_str("<body>\n");
    out.push_str("<h1>Feature Dashboard</h1>\n");
    out.push_str("<p>Generated from feature metadata and recorded per-corpus parity status. Run cargo xtask validate to compare against the current checkout. Do not hand-edit this file.</p>\n");
    out.push_str("<p class=\"legend\"><span class=\"ok\">&#10003;</span>passed <span class=\"bad\">&#10007;</span>failed <span class=\"unknown\">?</span>unknown <span class=\"na\">-</span>not required</p>\n");
    out.push_str("<div class=\"dashboard-wrap\">\n");
    out.push_str("<table id=\"feature-dashboard\">\n");
    out.push_str("<thead>\n<tr>");
    out.push_str("<th class=\"text\" data-sort-type=\"text\"><button class=\"sort\" type=\"button\">Feature</button></th>");
    out.push_str("<th class=\"text\" data-sort-type=\"text\"><button class=\"sort\" type=\"button\">Title</button></th>");
    out.push_str(&area_header());
    out.push_str(&rotated_header("Version", "Version", "number"));
    out.push_str(&rotated_header("Implemented", "Implemented", "number"));
    for corpus in VALIDATION_CORPORA {
        let (visible, title) = corpus_header(corpus.id, corpus.label, corpus_info);
        out.push_str(&rotated_header(&visible, &title, "number"));
    }
    out.push_str("</tr>\n</thead>\n<tbody>\n");
    for feature in features {
        let status = statuses.get(&feature.id);
        out.push_str(&format!(
            "<tr><td class=\"text\" data-sort-value=\"{0}\"><code>{0}</code></td>",
            escape_html(&feature.id)
        ));
        out.push_str(&format!(
            "<td class=\"text\" data-sort-value=\"{}\">{}</td>",
            escape_html(&feature.title),
            escape_html(&feature.title)
        ));
        out.push_str(&format!(
            "<td class=\"compact area\" data-sort-value=\"{}\">{}</td>",
            escape_html(&feature.area),
            escape_html(&feature.area)
        ));
        out.push_str(&format!(
            "<td class=\"compact\" data-sort-value=\"{}\">{}</td>",
            feature.version, feature.version
        ));
        out.push_str(&format!(
            "<td class=\"marker\" data-sort-value=\"{}\">{}</td>",
            bool_sort_value(feature.implemented),
            dashboard_bool_marker(feature.implemented)
        ));
        for corpus in VALIDATION_CORPORA {
            let applicable = corpus_info
                .get(corpus.id)
                .is_some_and(|info| info.feature_ids.contains(&feature.id));
            out.push_str(&dashboard_corpus_cell(
                feature, status, corpus.id, applicable,
            ));
        }
        out.push_str("</tr>\n");
    }
    out.push_str("</tbody>\n</table>\n</div>\n");
    out.push_str("<script>\n");
    out.push_str("(() => {\n");
    out.push_str("  const table = document.getElementById('feature-dashboard');\n");
    out.push_str("  const tbody = table.tBodies[0];\n");
    out.push_str("  const headers = Array.from(table.tHead.rows[0].cells);\n");
    out.push_str("  const value = (row, index, type) => {\n");
    out.push_str("    const raw = row.cells[index].dataset.sortValue || row.cells[index].textContent.trim();\n");
    out.push_str("    return type === 'number' ? Number(raw) : raw.toLocaleLowerCase();\n");
    out.push_str("  };\n");
    out.push_str("  headers.forEach((header, index) => {\n");
    out.push_str("    const button = header.querySelector('button.sort');\n");
    out.push_str("    if (!button) return;\n");
    out.push_str("    button.addEventListener('click', () => {\n");
    out.push_str("      const ascending = header.getAttribute('aria-sort') !== 'ascending';\n");
    out.push_str("      headers.forEach(other => other.removeAttribute('aria-sort'));\n");
    out.push_str(
        "      header.setAttribute('aria-sort', ascending ? 'ascending' : 'descending');\n",
    );
    out.push_str("      const type = header.dataset.sortType || 'text';\n");
    out.push_str("      const rows = Array.from(tbody.rows);\n");
    out.push_str("      rows.sort((left, right) => {\n");
    out.push_str("        const a = value(left, index, type);\n");
    out.push_str("        const b = value(right, index, type);\n");
    out.push_str("        if (a < b) return ascending ? -1 : 1;\n");
    out.push_str("        if (a > b) return ascending ? 1 : -1;\n");
    out.push_str("        return value(left, 0, 'text').localeCompare(value(right, 0, 'text'));\n");
    out.push_str("      });\n");
    out.push_str("      rows.forEach(row => tbody.appendChild(row));\n");
    out.push_str("    });\n");
    out.push_str("  });\n");
    out.push_str("})();\n");
    out.push_str("</script>\n");
    out.push_str("</body>\n</html>\n");
    out
}

pub(crate) fn corpus_header(
    corpus: &str,
    label: &str,
    corpus_info: &BTreeMap<String, CorpusDashboardInfo>,
) -> (String, String) {
    corpus_info
        .get(corpus)
        .map(|info| {
            (
                format!("{} (n={})", info.id, info.expected_count),
                format!("{}: {} expected case(s)", info.title, info.expected_count),
            )
        })
        .unwrap_or_else(|| (label.to_owned(), label.to_owned()))
}

pub(crate) fn rotated_header(label: &str, title: &str, sort_type: &str) -> String {
    format!(
        "<th class=\"compact rotated\" data-sort-type=\"{}\" title=\"{}\"><button class=\"sort\" type=\"button\" aria-label=\"Sort by {}\"><span class=\"rotated-label\">{}</span></button></th>",
        escape_html(sort_type),
        escape_html(title),
        escape_html(title),
        rotated_label_html(label)
    )
}

pub(crate) fn area_header() -> String {
    "<th class=\"compact area\" data-sort-type=\"text\" title=\"Area\"><button class=\"sort\" type=\"button\" aria-label=\"Sort by Area\">Area</button></th>".to_owned()
}

pub(crate) fn rotated_label_html(label: &str) -> String {
    if let Some((name, count)) = label.split_once(" (n=") {
        return format!(
            "<span class=\"rotated-name\">{}</span><br><span class=\"rotated-count\">(n={}</span>",
            escape_html(name),
            escape_html(count)
        );
    }
    format!("<span class=\"rotated-name\">{}</span>", escape_html(label))
}

pub(crate) fn dashboard_bool_marker(value: bool) -> &'static str {
    if value {
        "<span class=\"ok\" aria-label=\"passed\">&#10003;</span>"
    } else {
        "<span class=\"bad\" aria-label=\"failed\">&#10007;</span>"
    }
}

pub(crate) fn dashboard_corpus_cell(
    feature: &Feature,
    status: Option<&ValidationStatus>,
    corpus: &str,
    applicable: bool,
) -> String {
    let required = feature
        .validation_required
        .iter()
        .any(|required| required == corpus);
    let corpus_status = status.and_then(|status| status.corpora.get(corpus));
    if !required && !applicable && corpus_status.is_none() {
        return "<td class=\"marker\" data-sort-value=\"-1\"><span class=\"na\" aria-label=\"not required\" title=\"not required\">-</span></td>".to_owned();
    }
    if (required && recorded_corpus_passed(feature, status, corpus))
        || (!required && recorded_corpus_status_passed(corpus_status))
    {
        return "<td class=\"marker\" data-sort-value=\"1\"><span class=\"ok\" aria-label=\"passed\" title=\"recorded evidence passed\">&#10003;</span></td>".to_owned();
    }
    let marker = if let Some(corpus_status) = corpus_status {
        if !corpus_status.passed && corpus_status.failed_count > 0 {
            let title = corpus_status
                .first_failure
                .as_deref()
                .map(|failure| {
                    format!(
                        "{} non-passing case(s); first failure: {}",
                        corpus_status.failed_count, failure
                    )
                })
                .unwrap_or_else(|| format!("{} non-passing case(s)", corpus_status.failed_count));
            format!(
                "<span class=\"bad\" aria-label=\"failed: {} non-passing case(s)\" title=\"{}\">&#10007;<span class=\"count\">{}</span></span>",
                corpus_status.failed_count,
                escape_html(&title),
                corpus_status.failed_count
            )
        } else {
            let title = if corpus_status.passed {
                "recorded evidence is stale or incomplete"
            } else {
                "validation did not record fixture-level failures"
            };
            dashboard_unknown_marker(title)
        }
    } else {
        dashboard_unknown_marker("no recorded validation status")
    };
    format!("<td class=\"marker\" data-sort-value=\"0\">{marker}</td>")
}

pub(crate) fn dashboard_unknown_marker(title: &str) -> String {
    format!(
        "<span class=\"unknown\" aria-label=\"unknown\" title=\"{}\">?</span>",
        escape_html(title)
    )
}

pub(crate) fn bool_sort_value(value: bool) -> &'static str {
    if value {
        "1"
    } else {
        "0"
    }
}

pub(crate) fn escape_html(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

pub(crate) fn normalize_text_line_endings(text: &str) -> String {
    text.replace("\r\n", "\n").replace('\r', "\n")
}
