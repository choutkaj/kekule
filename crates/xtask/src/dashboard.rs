use crate::*;

pub(crate) fn dashboard(args: Vec<String>) -> Result<(), Box<dyn Error>> {
    let check = args.iter().any(|arg| arg == "--check");
    let features = read_features()?;
    let statuses = read_validation_statuses(&features)?;
    ensure_validation_flags_synced(&features, &statuses)?;
    let rendered = render_dashboard(&features, &statuses);
    let path = Path::new(DASHBOARD_PATH);

    if check {
        let existing = fs::read_to_string(path)?;
        if existing != rendered {
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
) -> String {
    let mut out = String::new();
    out.push_str("<!doctype html>\n");
    out.push_str("<html lang=\"en\">\n");
    out.push_str("<head>\n");
    out.push_str("<meta charset=\"utf-8\">\n");
    out.push_str("<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n");
    out.push_str("<title>Feature Dashboard</title>\n");
    out.push_str("<style>\n");
    out.push_str(":root { color-scheme: light dark; --border: #d0d7de; --head: #f6f8fa; --ok: #1a7f37; --bad: #cf222e; --muted: #656d76; --text: #24292f; --bg: #ffffff; }\n");
    out.push_str("@media (prefers-color-scheme: dark) { :root { --border: #30363d; --head: #161b22; --ok: #3fb950; --bad: #ff7b72; --muted: #8b949e; --text: #c9d1d9; --bg: #0d1117; } }\n");
    out.push_str("body { margin: 24px; background: var(--bg); color: var(--text); font: 14px/1.4 system-ui, -apple-system, Segoe UI, sans-serif; }\n");
    out.push_str("h1 { margin: 0 0 4px; font-size: 24px; }\n");
    out.push_str("p { margin: 0 0 18px; color: var(--muted); }\n");
    out.push_str(".dashboard-wrap { overflow-x: auto; }\n");
    out.push_str("table { border-collapse: collapse; width: 100%; min-width: 980px; }\n");
    out.push_str(
        "th, td { border: 1px solid var(--border); padding: 6px 8px; vertical-align: middle; }\n",
    );
    out.push_str("thead th { position: sticky; top: 0; z-index: 1; height: 112px; background: var(--head); white-space: nowrap; }\n");
    out.push_str("tbody tr:nth-child(even) { background: color-mix(in srgb, var(--head) 45%, transparent); }\n");
    out.push_str("th.text, td.text { text-align: left; }\n");
    out.push_str("th.compact, td.compact, th.rotated, td.marker { text-align: center; }\n");
    out.push_str(
        "th.rotated { width: 42px; min-width: 42px; padding: 0; vertical-align: bottom; }\n",
    );
    out.push_str("th.rotated button { height: 108px; width: 42px; padding: 0; display: flex; align-items: flex-end; justify-content: center; }\n");
    out.push_str("th.rotated span { display: inline-block; transform: rotate(-60deg); transform-origin: bottom left; width: 96px; text-align: left; }\n");
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
    out.push_str(".na { color: var(--muted); }\n");
    out.push_str("code { font-family: ui-monospace, SFMono-Regular, Consolas, monospace; font-size: 13px; }\n");
    out.push_str("</style>\n");
    out.push_str("</head>\n");
    out.push_str("<body>\n");
    out.push_str("<h1>Feature Dashboard</h1>\n");
    out.push_str("<p>Generated from feature metadata and validation status. Do not hand-edit this file.</p>\n");
    out.push_str("<div class=\"dashboard-wrap\">\n");
    out.push_str("<table id=\"feature-dashboard\">\n");
    out.push_str("<thead>\n<tr>");
    out.push_str("<th class=\"text\" data-sort-type=\"text\"><button class=\"sort\" type=\"button\">Feature</button></th>");
    out.push_str("<th class=\"text\" data-sort-type=\"text\"><button class=\"sort\" type=\"button\">Title</button></th>");
    out.push_str("<th class=\"compact rotated\" data-sort-type=\"text\"><button class=\"sort\" type=\"button\"><span>Area</span></button></th>");
    out.push_str("<th class=\"compact rotated\" data-sort-type=\"number\"><button class=\"sort\" type=\"button\"><span>Version</span></button></th>");
    out.push_str("<th class=\"compact rotated\" data-sort-type=\"number\"><button class=\"sort\" type=\"button\"><span>Implemented</span></button></th>");
    for (_, label) in VALIDATION_CORPORA {
        out.push_str(&format!(
            "<th class=\"rotated\" data-sort-type=\"number\"><button class=\"sort\" type=\"button\"><span>{}</span></button></th>",
            escape_html(label)
        ));
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
            "<td class=\"compact\" data-sort-value=\"{}\">{}</td>",
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
            dashboard_marker(Some(feature.implemented))
        ));
        for (corpus, _) in VALIDATION_CORPORA {
            let marker = if feature
                .validation_required
                .iter()
                .any(|required| required == corpus)
            {
                Some(corpus_passed(feature, status, corpus))
            } else {
                None
            };
            out.push_str(&format!(
                "<td class=\"marker\" data-sort-value=\"{}\">{}</td>",
                optional_bool_sort_value(marker),
                dashboard_marker(marker)
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

pub(crate) fn dashboard_marker(value: Option<bool>) -> &'static str {
    match value {
        Some(true) => "<span class=\"ok\" aria-label=\"passed\">&#10003;</span>",
        Some(false) => "<span class=\"bad\" aria-label=\"failed\">&#10007;</span>",
        None => "<span class=\"na\" aria-label=\"not required\">-</span>",
    }
}

pub(crate) fn bool_sort_value(value: bool) -> &'static str {
    if value {
        "1"
    } else {
        "0"
    }
}

pub(crate) fn optional_bool_sort_value(value: Option<bool>) -> &'static str {
    match value {
        Some(true) => "1",
        Some(false) => "0",
        None => "-1",
    }
}

pub(crate) fn escape_html(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}
