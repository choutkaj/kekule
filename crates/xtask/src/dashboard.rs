use crate::*;

#[derive(Debug, Clone, Copy)]
pub(crate) struct DashboardSection {
    pub(crate) title: &'static str,
    pub(crate) table_id: &'static str,
    pub(crate) domain: FeatureDomain,
    pub(crate) include_corpora: bool,
}

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
    out.push_str(":root { color-scheme: light dark; --border: #d0d7de; --head: #f6f8fa; --ok: #1a7f37; --bad: #cf222e; --unknown: #9a6700; --muted: #656d76; --text: #24292f; --bg: #ffffff; --planned: #656d76; --experimental: #9a6700; --supported: #1a7f37; --deprecated: #8250df; }\n");
    out.push_str("@media (prefers-color-scheme: dark) { :root { --border: #30363d; --head: #161b22; --ok: #3fb950; --bad: #ff7b72; --unknown: #d29922; --muted: #8b949e; --text: #c9d1d9; --bg: #0d1117; --planned: #8b949e; --experimental: #d29922; --supported: #3fb950; --deprecated: #d2a8ff; } }\n");
    out.push_str("body { margin: 24px; background: var(--bg); color: var(--text); font: 14px/1.4 system-ui, -apple-system, Segoe UI, sans-serif; }\n");
    out.push_str("h1 { margin: 0 0 4px; font-size: 24px; }\n");
    out.push_str("h2 { margin: 30px 0 4px; font-size: 20px; }\n");
    out.push_str("p { margin: 0 0 18px; color: var(--muted); }\n");
    out.push_str(".section-description { margin-bottom: 6px; }\n");
    out.push_str(".reference { margin-bottom: 4px; color: var(--text); }\n");
    out.push_str(".supplemental { margin-bottom: 10px; font-size: 12px; }\n");
    out.push_str(".dashboard-wrap { overflow-x: auto; }\n");
    out.push_str(".graph-wrap { overflow-x: auto; padding: 8px; border: 1px solid var(--border); border-radius: 8px; background: color-mix(in srgb, var(--head) 55%, transparent); }\n");
    out.push_str(".feature-graph { display: block; min-width: 100%; height: auto; }\n");
    out.push_str(".graph-edge { fill: none; stroke: var(--muted); stroke-width: 1.2; stroke-opacity: .42; }\n");
    out.push_str(".graph-node rect { stroke-width: 1.5; rx: 6; }\n");
    out.push_str(".graph-node text { fill: var(--text); font: 11px ui-monospace, SFMono-Regular, Consolas, monospace; pointer-events: none; }\n");
    out.push_str(".graph-layer-label { fill: var(--muted); font: 11px system-ui, -apple-system, Segoe UI, sans-serif; }\n");
    out.push_str(".graph-node.status-planned rect { fill: color-mix(in srgb, var(--planned) 12%, var(--bg)); stroke: var(--planned); stroke-dasharray: 4 3; }\n");
    out.push_str(".graph-node.status-experimental rect { fill: color-mix(in srgb, var(--experimental) 13%, var(--bg)); stroke: var(--experimental); }\n");
    out.push_str(".graph-node.status-supported rect { fill: color-mix(in srgb, var(--supported) 12%, var(--bg)); stroke: var(--supported); }\n");
    out.push_str(".graph-node.status-deprecated rect { fill: color-mix(in srgb, var(--deprecated) 13%, var(--bg)); stroke: var(--deprecated); }\n");
    out.push_str("table { border-collapse: collapse; width: 100%; min-width: 980px; }\n");
    out.push_str("table.infrastructure-table { min-width: 720px; }\n");
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
    out.push_str(".feature-status { display: inline-block; padding: 1px 6px; border: 1px solid currentColor; border-radius: 999px; font-size: 12px; font-weight: 650; white-space: nowrap; }\n");
    out.push_str(".status-planned { color: var(--planned); }\n");
    out.push_str(".status-experimental { color: var(--experimental); }\n");
    out.push_str(".status-supported { color: var(--supported); }\n");
    out.push_str(".status-deprecated { color: var(--deprecated); }\n");
    out.push_str("code { font-family: ui-monospace, SFMono-Regular, Consolas, monospace; font-size: 13px; }\n");
    out.push_str("</style>\n");
    out.push_str("</head>\n");
    out.push_str("<body>\n");
    out.push_str("<h1>Feature Dashboard</h1>\n");
    out.push_str("<p>Generated from feature metadata and recorded per-corpus parity status. Run cargo xtask validate to compare against the current checkout. Do not hand-edit this file.</p>\n");
    out.push_str("<p class=\"legend\"><span class=\"ok\">&#10003;</span>passed <span class=\"bad\">&#10007;</span>failed <span class=\"unknown\">?</span>unknown <span class=\"na\">-</span>not required</p>\n");
    out.push_str("<p>Features shared by both molecular domains are intentionally shown in both chemistry tables. The mixed smoke corpus also appears in both tables; only applicable feature rows carry parity status.</p>\n");
    render_feature_graph(&mut out, features);
    render_dashboard_section(
        &mut out,
        DashboardSection {
            title: "Small molecules",
            table_id: "small-molecules-dashboard",
            domain: FeatureDomain::SmallMolecule,
            include_corpora: true,
        },
        features,
        statuses,
        corpus_info,
    );
    render_dashboard_section(
        &mut out,
        DashboardSection {
            title: "Macromolecules",
            table_id: "macromolecules-dashboard",
            domain: FeatureDomain::Macromolecule,
            include_corpora: true,
        },
        features,
        statuses,
        corpus_info,
    );
    render_dashboard_section(
        &mut out,
        DashboardSection {
            title: "Infrastructure and harness",
            table_id: "infrastructure-dashboard",
            domain: FeatureDomain::Infrastructure,
            include_corpora: false,
        },
        features,
        statuses,
        corpus_info,
    );
    out.push_str("<script>\n");
    out.push_str("(() => {\n");
    out.push_str("  document.querySelectorAll('table.feature-dashboard').forEach(table => {\n");
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
    out.push_str("  });\n");
    out.push_str("})();\n");
    out.push_str("</script>\n");
    out.push_str("</body>\n</html>\n");
    out
}

pub(crate) fn render_feature_graph(out: &mut String, features: &[Feature]) {
    const NODE_WIDTH: usize = 210;
    const NODE_HEIGHT: usize = 34;
    const COLUMN_GAP: usize = 46;
    const ROW_GAP: usize = 18;
    const MARGIN_X: usize = 24;
    const MARGIN_TOP: usize = 40;
    const MARGIN_BOTTOM: usize = 22;

    let layers = feature_dependency_layers(features)
        .expect("dashboard features should form a validated dependency graph");
    let max_rows = layers.iter().map(Vec::len).max().unwrap_or(0).max(1);
    let width =
        MARGIN_X * 2 + layers.len() * NODE_WIDTH + layers.len().saturating_sub(1) * COLUMN_GAP;
    let height =
        MARGIN_TOP + max_rows * NODE_HEIGHT + max_rows.saturating_sub(1) * ROW_GAP + MARGIN_BOTTOM;
    let mut positions = BTreeMap::new();
    for (layer_index, layer) in layers.iter().enumerate() {
        let layer_height = layer.len() * NODE_HEIGHT + layer.len().saturating_sub(1) * ROW_GAP;
        let available_height = max_rows * NODE_HEIGHT + max_rows.saturating_sub(1) * ROW_GAP;
        let layer_offset = (available_height.saturating_sub(layer_height)) / 2;
        for (row_index, feature) in layer.iter().enumerate() {
            positions.insert(
                feature.id.as_str(),
                (
                    MARGIN_X + layer_index * (NODE_WIDTH + COLUMN_GAP),
                    MARGIN_TOP + layer_offset + row_index * (NODE_HEIGHT + ROW_GAP),
                ),
            );
        }
    }

    out.push_str("<section id=\"feature-dependency-graph\">\n");
    out.push_str("<h2>Feature dependency graph</h2>\n");
    out.push_str("<p class=\"section-description\">Generated from <code>depends_on</code>. Arrows point from a prerequisite to the feature that depends on it; columns are deterministic dependency layers.</p>\n");
    out.push_str("<p class=\"legend\">");
    for status in [
        FeatureStatus::Planned,
        FeatureStatus::Experimental,
        FeatureStatus::Supported,
        FeatureStatus::Deprecated,
    ] {
        out.push_str(&dashboard_status_marker(status));
        out.push(' ');
    }
    out.push_str("</p>\n");
    out.push_str("<div class=\"graph-wrap\">\n");
    out.push_str(&format!(
        "<svg class=\"feature-graph\" viewBox=\"0 0 {width} {height}\" width=\"{width}\" height=\"{height}\" role=\"img\" aria-labelledby=\"feature-graph-title feature-graph-description\">\n"
    ));
    out.push_str("<title id=\"feature-graph-title\">Feature dependency graph</title>\n");
    out.push_str("<desc id=\"feature-graph-description\">Directed acyclic graph of repository features. Arrows lead from dependencies to dependents.</desc>\n");
    out.push_str("<defs><marker id=\"feature-graph-arrow\" markerWidth=\"7\" markerHeight=\"7\" refX=\"6\" refY=\"3.5\" orient=\"auto\" markerUnits=\"strokeWidth\"><path d=\"M0,0 L7,3.5 L0,7 Z\" fill=\"#8c959f\"/></marker></defs>\n");
    for feature in features {
        let Some((end_x, end_y)) = positions.get(feature.id.as_str()).copied() else {
            continue;
        };
        for dependency in &feature.depends_on {
            let Some((start_x, start_y)) = positions.get(dependency.as_str()).copied() else {
                continue;
            };
            let start_x = start_x + NODE_WIDTH;
            let start_y = start_y + NODE_HEIGHT / 2;
            let end_y = end_y + NODE_HEIGHT / 2;
            let bend_x = (start_x + end_x) / 2;
            out.push_str(&format!(
                "<path class=\"graph-edge\" marker-end=\"url(#feature-graph-arrow)\" d=\"M {start_x} {start_y} C {bend_x} {start_y}, {bend_x} {end_y}, {end_x} {end_y}\"/>\n"
            ));
        }
    }
    for (layer_index, layer) in layers.iter().enumerate() {
        let x = MARGIN_X + layer_index * (NODE_WIDTH + COLUMN_GAP);
        out.push_str(&format!(
            "<text class=\"graph-layer-label\" x=\"{x}\" y=\"18\">layer {layer_index}</text>\n"
        ));
        for feature in layer {
            let (x, y) = positions[feature.id.as_str()];
            let id = escape_html(&feature.id);
            let title = escape_html(&format!(
                "{} — {}; depends on: {}",
                feature.title,
                feature.status.as_str(),
                if feature.depends_on.is_empty() {
                    "none".to_owned()
                } else {
                    feature.depends_on.join(", ")
                }
            ));
            out.push_str(&format!(
                "<a href=\"./{id}/feature.md\"><g class=\"graph-node status-{}\" transform=\"translate({x} {y})\"><title>{title}</title><rect width=\"{NODE_WIDTH}\" height=\"{NODE_HEIGHT}\" rx=\"6\"/><text x=\"10\" y=\"21\">{id}</text></g></a>\n",
                feature.status.as_str()
            ));
        }
    }
    out.push_str("</svg>\n</div>\n</section>\n");
}

pub(crate) fn render_dashboard_section(
    out: &mut String,
    section: DashboardSection,
    features: &[Feature],
    statuses: &BTreeMap<String, ValidationStatus>,
    corpus_info: &BTreeMap<String, CorpusDashboardInfo>,
) {
    let corpora = if section.include_corpora {
        dashboard_corpora(section.domain, corpus_info)
    } else {
        Vec::new()
    };
    out.push_str(&format!(
        "<section>\n<h2>{}</h2>\n",
        escape_html(section.title)
    ));
    if section.include_corpora {
        let description = match section.domain {
            FeatureDomain::SmallMolecule => {
                "Small-molecule features with small-molecule validation corpora."
            }
            FeatureDomain::Macromolecule => {
                "Macromolecular features with PDB-derived validation corpora."
            }
            FeatureDomain::Infrastructure => unreachable!("infrastructure has no corpora"),
        };
        out.push_str(&format!(
            "<p class=\"section-description\">{}</p>\n",
            escape_html(description)
        ));
        out.push_str(&dashboard_reference_summary(
            section.domain,
            features,
            &corpora,
            corpus_info,
        ));
    } else {
        out.push_str("<p class=\"section-description\">Repository feature-registry and validation-harness capabilities. These rows do not use an external chemistry reference codebase.</p>\n");
    }
    out.push_str("<div class=\"dashboard-wrap\">\n");
    let table_class = if section.include_corpora {
        "feature-dashboard"
    } else {
        "feature-dashboard infrastructure-table"
    };
    out.push_str(&format!(
        "<table id=\"{}\" class=\"{}\">\n",
        escape_html(section.table_id),
        table_class
    ));
    out.push_str("<thead>\n<tr>");
    out.push_str("<th class=\"text\" data-sort-type=\"text\"><button class=\"sort\" type=\"button\">Feature</button></th>");
    out.push_str("<th class=\"text\" data-sort-type=\"text\"><button class=\"sort\" type=\"button\">Title</button></th>");
    out.push_str(&area_header());
    out.push_str(&rotated_header("Version", "Version", "number"));
    out.push_str("<th class=\"compact\" data-sort-type=\"number\" title=\"Release status\"><button class=\"sort\" type=\"button\" aria-label=\"Sort by Status\">Status</button></th>");
    for corpus in &corpora {
        let (visible, title) = corpus_header(corpus.id, corpus.label, corpus_info);
        out.push_str(&rotated_header(&visible, &title, "number"));
    }
    out.push_str("</tr>\n</thead>\n<tbody>\n");
    for feature in features
        .iter()
        .filter(|feature| feature.domains.contains(&section.domain))
    {
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
            "<td class=\"compact\" data-sort-value=\"{}\">{}</td>",
            feature.status.sort_value(),
            dashboard_status_marker(feature.status)
        ));
        for corpus in &corpora {
            let info = corpus_info
                .get(corpus.id)
                .expect("dashboard corpus metadata should be available");
            let reference = info.features.get(&feature.id);
            let corpus_status = status.and_then(|status| status.corpora.get(corpus.id));
            let reference_tool = reference
                .map(|reference| reference.reference_tool.as_str())
                .or_else(|| corpus_status.map(|status| status.reference_tool.as_str()));
            let domain_applicable = match info.kind {
                CorpusKind::SmallMolecule => section.domain == FeatureDomain::SmallMolecule,
                CorpusKind::Macromolecule => section.domain == FeatureDomain::Macromolecule,
                CorpusKind::Mixed => reference_tool.map_or_else(
                    || feature.domains.len() == 1 && feature.domains.contains(&section.domain),
                    |tool| dashboard_reference_applies_to_domain(section.domain, info.kind, tool),
                ),
            };
            out.push_str(&dashboard_corpus_cell(
                feature,
                status,
                corpus.id,
                reference,
                domain_applicable,
            ));
        }
        out.push_str("</tr>\n");
    }
    out.push_str("</tbody>\n</table>\n</div>\n</section>\n");
}

pub(crate) fn dashboard_corpora(
    domain: FeatureDomain,
    corpus_info: &BTreeMap<String, CorpusDashboardInfo>,
) -> Vec<&'static ValidationCorpus> {
    VALIDATION_CORPORA
        .iter()
        .filter(|corpus| {
            corpus_info.get(corpus.id).is_some_and(|info| match domain {
                FeatureDomain::SmallMolecule => {
                    matches!(info.kind, CorpusKind::SmallMolecule | CorpusKind::Mixed)
                }
                FeatureDomain::Macromolecule => {
                    matches!(info.kind, CorpusKind::Macromolecule | CorpusKind::Mixed)
                }
                FeatureDomain::Infrastructure => false,
            })
        })
        .collect()
}

pub(crate) fn dashboard_reference_summary(
    domain: FeatureDomain,
    features: &[Feature],
    corpora: &[&ValidationCorpus],
    corpus_info: &BTreeMap<String, CorpusDashboardInfo>,
) -> String {
    let feature_ids = features
        .iter()
        .filter(|feature| feature.domains.contains(&domain))
        .map(|feature| feature.id.as_str())
        .collect::<BTreeSet<_>>();
    let mut codebases = BTreeSet::new();
    let mut dssp_executables = BTreeSet::new();
    let mut supplemental = BTreeSet::new();
    for corpus in corpora {
        let Some(info) = corpus_info.get(corpus.id) else {
            continue;
        };
        for (feature_id, reference) in &info.features {
            if !feature_ids.contains(feature_id.as_str()) {
                continue;
            }
            if !dashboard_reference_applies_to_domain(domain, info.kind, &reference.reference_tool)
            {
                continue;
            }
            match reference.reference_tool.as_str() {
                "rdkit" => {
                    codebases.insert(dashboard_reference_label(
                        &reference.reference_tool,
                        &reference.reference_version,
                    ));
                }
                "biopython" => {
                    if let Some((biopython, mkdssp)) =
                        reference.reference_version.split_once(" / mkdssp version ")
                    {
                        codebases.insert(version_label("Biopython", biopython));
                        dssp_executables.insert(version_label("mkdssp", mkdssp));
                    } else {
                        codebases.insert(dashboard_reference_label(
                            &reference.reference_tool,
                            &reference.reference_version,
                        ));
                    }
                }
                tool if tool.ends_with("-manual-semantic") => {
                    supplemental.insert(reference.reference_version.clone());
                }
                _ => {}
            }
        }
    }
    let mut out = String::new();
    out.push_str("<p class=\"reference\"><strong>Reference codebase:</strong> ");
    if codebases.is_empty() {
        out.push_str("none recorded");
    } else {
        out.push_str(&escape_html(
            &codebases.into_iter().collect::<Vec<_>>().join("; "),
        ));
    }
    if !dssp_executables.is_empty() {
        out.push_str("; <strong>DSSP executable:</strong> ");
        out.push_str(&escape_html(
            &dssp_executables.into_iter().collect::<Vec<_>>().join("; "),
        ));
    }
    out.push_str("</p>\n");
    if !supplemental.is_empty() {
        out.push_str(&format!(
            "<p class=\"supplemental\"><strong>Supplemental semantic references:</strong> {}</p>\n",
            escape_html(&supplemental.into_iter().collect::<Vec<_>>().join("; "))
        ));
    }
    out
}

pub(crate) fn dashboard_reference_applies_to_domain(
    domain: FeatureDomain,
    corpus_kind: CorpusKind,
    reference_tool: &str,
) -> bool {
    match corpus_kind {
        CorpusKind::SmallMolecule => domain == FeatureDomain::SmallMolecule,
        CorpusKind::Macromolecule => domain == FeatureDomain::Macromolecule,
        CorpusKind::Mixed => match reference_tool {
            "biopython" => domain == FeatureDomain::Macromolecule,
            "rdkit" | "planned-rdkit" => domain == FeatureDomain::SmallMolecule,
            tool if tool.ends_with("-manual-semantic") => domain == FeatureDomain::SmallMolecule,
            _ => false,
        },
    }
}

pub(crate) fn version_label(name: &str, version: &str) -> String {
    let without_name = version.strip_prefix(name).map(str::trim).unwrap_or(version);
    let without_v = without_name.strip_prefix('v').unwrap_or(without_name);
    format!("{name} v{without_v}")
}

pub(crate) fn dashboard_reference_label(tool: &str, version: &str) -> String {
    match tool {
        "rdkit" => version_label("RDKit", version),
        "biopython" => {
            if let Some((biopython, mkdssp)) = version.split_once(" / mkdssp version ") {
                format!(
                    "{} / {}",
                    version_label("Biopython", biopython),
                    version_label("mkdssp", mkdssp)
                )
            } else {
                version_label("Biopython", version)
            }
        }
        tool if tool.ends_with("-manual-semantic") => version.to_owned(),
        _ if version
            .to_ascii_lowercase()
            .starts_with(&tool.to_ascii_lowercase()) =>
        {
            version.to_owned()
        }
        _ => format!("{tool} {version}"),
    }
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

pub(crate) fn dashboard_status_marker(status: FeatureStatus) -> String {
    let status = status.as_str();
    format!("<span class=\"feature-status status-{status}\">{status}</span>")
}

pub(crate) fn dashboard_corpus_cell(
    feature: &Feature,
    status: Option<&ValidationStatus>,
    corpus: &str,
    manifest_reference: Option<&CorpusFeatureDashboardInfo>,
    domain_applicable: bool,
) -> String {
    if !domain_applicable {
        return "<td class=\"marker\" data-sort-value=\"-1\"><span class=\"na\" aria-label=\"not required\" title=\"not required\">-</span></td>".to_owned();
    }
    let required = feature
        .validation_required
        .iter()
        .any(|required| required == corpus);
    let corpus_status = status.and_then(|status| status.corpora.get(corpus));
    if !required && manifest_reference.is_none() && corpus_status.is_none() {
        return "<td class=\"marker\" data-sort-value=\"-1\"><span class=\"na\" aria-label=\"not required\" title=\"not required\">-</span></td>".to_owned();
    }
    let reference = corpus_status
        .map(|status| {
            (
                status.reference_tool.as_str(),
                status.reference_version.as_str(),
            )
        })
        .or_else(|| {
            manifest_reference.map(|reference| {
                (
                    reference.reference_tool.as_str(),
                    reference.reference_version.as_str(),
                )
            })
        });
    if (required && recorded_corpus_passed(feature, status, corpus))
        || (!required && recorded_corpus_status_passed(corpus_status))
    {
        let title = dashboard_cell_title("recorded evidence passed", reference);
        return format!(
            "<td class=\"marker\" data-sort-value=\"1\"><span class=\"ok\" aria-label=\"passed\" title=\"{}\">&#10003;</span></td>",
            escape_html(&title)
        );
    }
    let marker = if let Some(corpus_status) = corpus_status {
        if !corpus_status.passed && corpus_status.failed_count > 0 {
            let failure_title = corpus_status
                .first_failure
                .as_deref()
                .map(|failure| {
                    format!(
                        "{} non-passing case(s); first failure: {}",
                        corpus_status.failed_count, failure
                    )
                })
                .unwrap_or_else(|| format!("{} non-passing case(s)", corpus_status.failed_count));
            let title = dashboard_cell_title(&failure_title, reference);
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
            dashboard_unknown_marker(&dashboard_cell_title(title, reference))
        }
    } else {
        dashboard_unknown_marker(&dashboard_cell_title(
            "no recorded validation status",
            reference,
        ))
    };
    format!("<td class=\"marker\" data-sort-value=\"0\">{marker}</td>")
}

pub(crate) fn dashboard_cell_title(base: &str, reference: Option<(&str, &str)>) -> String {
    reference.map_or_else(
        || base.to_owned(),
        |(tool, version)| {
            format!(
                "{base}; reference: {}",
                dashboard_reference_label(tool, version)
            )
        },
    )
}

pub(crate) fn dashboard_unknown_marker(title: &str) -> String {
    format!(
        "<span class=\"unknown\" aria-label=\"unknown\" title=\"{}\">?</span>",
        escape_html(title)
    )
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
