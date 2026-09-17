"use strict";
let data = JSON.parse(document.getElementById("benchmark-data").textContent);
const $ = id => document.getElementById(id);
const number = value => value.toLocaleString("en-US");
const node = (tag, cls, text) => {
  const element = document.createElement(tag);
  if (cls) element.className = cls;
  if (text !== undefined) element.textContent = text;
  return element;
};
const outcomes = ["agree", "disagree", "error"];
let active = 0;
const current = () => data.runs[active];
const index = () => new Map(current().results.map(row => [`${row.dataset}/${row.feature}`, row]));
function segments(container, counts, classes, total) {
  container.replaceChildren();
  counts.forEach((value, i) => {
    if (!value || !total) return;
    const bar = node("span", `segment ${classes[i]}`);
    bar.style.width = `${100 * value / total}%`;
    container.append(bar);
  });
}
function renderCorpora() {
  $("corpora").replaceChildren(...data.catalog.datasets.map(dataset => {
    const row = node("tr");
    const name = node("th", "", dataset.id); name.scope = "row";
    row.append(name, node("td", "numeric", number(dataset.source_ids)),
      node("td", "", dataset.formats.join(", ") || "—"),
      node("td", "", dataset.description));
    return row;
  }));
}
function renderMatrix() {
  const rows = current() ? index() : new Map();
  const header = node("tr"); header.append(node("th", "", "Feature"));
  for (const dataset of data.catalog.datasets) {
    const th = node("th", "", dataset.id); th.scope = "col";
    header.append(th);
  }
  $("matrix").tHead.replaceChildren(header);
  const body = $("matrix").tBodies[0]; body.replaceChildren();
  for (const feature of data.catalog.features.filter(name => name.includes($("search").value.toLowerCase().trim()))) {
    const tr = node("tr"); const label = node("th", "feature-label", feature);
    label.scope = "row"; tr.append(label);
    for (const dataset of data.catalog.datasets) {
      const key = `${dataset.id}/${feature}`; const row = rows.get(key);
      const td = node("td"); const cell = node("div", `cell${row ? "" : " not-run"}`);
      cell.setAttribute("role", "img");
      cell.setAttribute("aria-label", `${dataset.id}, ${feature}: ${row ? `${row.agrees} agreeing of ${row.cases} applicable cases; ${row.coverage} selection` : "not reported in this run"}`);
      if (row) {
        const top = node("div", "cell-top"); top.append(node("span", "", row.cases ? `${number(row.agrees)} / ${number(row.cases)}` : "No supplied format"));
        top.append(node("span", `scope ${row.coverage === "sampled" ? "sample-label" : row.coverage === "stale" ? "stale-label" : ""}`, row.coverage));
        cell.append(top);
        if (row.cases) { const track = node("div", "cell-track"); segments(track, [row.agrees, row.disagrees, row.errors], outcomes, row.cases); cell.append(track); }
        else cell.append(node("span", "cell-note", `${row.not_applicable} not applicable`));
        cell.title = `${row.source_ids} selected IDs; ${row.agrees} agrees; ${row.disagrees} disagreements; ${row.errors} errors; ${row.not_applicable} not applicable`;
      } else cell.append(node("span", "", "— Unrun"));
      td.append(cell); tr.append(td);
    }
    body.append(tr);
  }
  $("no-matches").hidden = body.rows.length !== 0;
}
function renderProvenance() {
  if (!current()) { $("provenance").textContent = "No benchmark runs recorded."; return; }
  const run = current(); const identity = run.implementation;
  const list = node("dl");
  const entries = [["Source report", run.name], ["Report SHA-256", run.sha256], ["Run completion", run.complete ? "Finished" : "Incomplete"],
    ["Comparison status", run.passed ? "All measured cases agree" : "Not passing"], ["Git revision", identity.revision || "Unavailable"],
    ["Worktree", identity.dirty === null ? "Unknown" : identity.dirty ? "Modified" : "Clean"],
    ...["working_tree_status_sha256", "executable_sha256", "reference_code_sha256", "contract_sha256"].map(key => [key.replaceAll("_", " "), identity[key] || "Unavailable"])];
  for (const [label, value] of entries) list.append(node("dt", "", label), node("dd", "", value));
  $("provenance").replaceChildren(list, node("p", "caption", "No run timestamp is inferred from file modification time. The report hash identifies the exact source; local paths and case payloads are omitted from this export."));
}
function renderRuns() {
  $("run").replaceChildren(...data.runs.map((run, i) => {
    const label = run.started_at_unix_ms === null ? run.name :
      new Date(run.started_at_unix_ms).toISOString().replace("T", " ").replace("Z", " UTC");
    const option = node("option", "", label); option.value = i; return option;
  }));
  if (!data.runs.length) $("run").append(node("option", "", "No runs yet"));
  $("run").value = "0";
  $("run").disabled = $("download").disabled = !data.runs.length;
}
function render() { renderMatrix(); renderProvenance(); }
$("run").addEventListener("change", () => { active = Number($("run").value); render(); });
$("search").addEventListener("input", renderMatrix);
$("download").addEventListener("click", () => {
  if (!current()) return;
  const url = URL.createObjectURL(new Blob([JSON.stringify({export_schema: 1, run: current()}, null, 2)], {type: "application/json"}));
  const link = node("a"); link.href = url; link.download = "kekule-benchmark-summary.json";
  document.body.append(link); link.click(); link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
});
renderRuns();
renderCorpora();
render();
if (data.live) {
  window.updateKekuleBenchmarks = next => {
    if (next.signature === data.signature) return;
    data = next; active = 0; $("search").value = "";
    renderRuns(); renderCorpora(); render();
  };
  function checkForRuns() {
    const script = document.createElement("script");
    script.src = `dashboard-data.js?${Date.now()}`;
    script.onload = script.onerror = () => { script.remove(); setTimeout(checkForRuns, 5000); };
    document.body.append(script);
  }
  checkForRuns();
}
