// Browser-independent regression for the local report feed. No external packages.
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const assert = require('node:assert/strict');

class Element {
  constructor(tag) { this.tag = tag; this.children = []; this.attributes = {}; this.events = {}; this.style = {}; this.value = ''; this._text = ''; }
  set textContent(value) { this._text = String(value); this.children = []; }
  get textContent() { return this._text + this.children.map(child => child.textContent).join(''); }
  append(...children) { children.forEach(child => { child.parent = this; this.children.push(child); }); }
  replaceChildren(...children) { this._text = ''; this.children = []; this.append(...children); }
  setAttribute(name, value) { this.attributes[name] = value; }
  addEventListener(name, action) { this.events[name] = action; }
  get rows() { return this.children.filter(child => child.tag === 'tr'); }
  click() { this.events.click?.(); }
  remove() { this.parent.children = this.parent.children.filter(child => child !== this); }
}
const html = fs.readFileSync(path.join(__dirname, 'dashboard/index.html'), 'utf8');
const elements = Object.fromEntries([...html.matchAll(/id="([^"]+)"/g)].map(match => [match[1], new Element('div')]));
elements.matrix.tHead = new Element('thead'); elements.matrix.tBodies = [new Element('tbody')];
const catalog = {datasets: [{id: 'example', source_ids: 4, formats: ['SMILES'], description: 'Unit regression'}], features: ['io.smiles.parse']};
const snapshot = (signature, runs) => ({catalog, runs, signature, live: true});
const run = (name, started, agrees) => ({
  name, started_at_unix_ms: started, sha256: name, complete: true, passed: agrees === 4,
  implementation: {revision: null, dirty: false},
  results: [{dataset: 'example', feature: 'io.smiles.parse', source_ids: 4, cases: 4,
    agrees, exact_agrees: 1, disagrees: 4 - agrees, errors: 0, not_applicable: 0, coverage: 'full'}],
});
elements['benchmark-data'].textContent = JSON.stringify(snapshot('empty', []));
const document = {body: new Element('body'), createElement: tag => new Element(tag),
  getElementById: id => { assert.ok(elements[id], `missing element ${id}`); return elements[id]; }};
let download;
const timers = [];
const context = {document, window: {}, Blob, Date,
  URL: {createObjectURL(blob) { download = blob; return 'blob:test'; }, revokeObjectURL() {}},
  setTimeout(action, delay) { timers.push({action, delay}); }};
vm.runInNewContext(fs.readFileSync(path.join(__dirname, 'dashboard/app.js'), 'utf8'), context);
assert.equal(elements.run.disabled, true);
assert.equal(elements.download.disabled, true);
assert.match(elements.provenance.textContent, /No benchmark runs/);
assert.match(document.body.children[0].src, /^dashboard-data\.js\?\d+$/);
const update = context.window.updateKekuleBenchmarks;
const older = run('older.json', 1000, 2);
const newer = run('newer.json', 2000, 3);
update(snapshot('first', [newer, older]));
assert.equal(elements.run.disabled, false);
assert.equal(elements.run.children.length, 2);
assert.equal(elements.run.value, '0');
assert.match(elements.run.children[0].textContent, /1970-01-01 00:00:02.000 UTC/);
const cell = () => elements.matrix.tBodies[0].rows[0].children[1].children[0];
assert.match(cell().textContent, /3 \/ 4/);
assert.deepEqual(cell().children[1].children.map(bar => [bar.className, bar.style.width]),
  [['segment agree', '75%'], ['segment disagree', '25%']]);
elements.run.value = '1'; elements.run.events.change();
elements.search.value = 'io.smiles';
update(snapshot('first', [newer, older]));
assert.equal(elements.run.value, '1');
assert.equal(elements.search.value, 'io.smiles');
assert.match(cell().textContent, /2 \/ 4/);
document.body.children[0].onload();
assert.equal(document.body.children.length, 0);
assert.equal(timers[0].delay, 5000);
timers.shift().action();
assert.equal(document.body.children.length, 1);
document.body.children[0].onerror();
assert.equal(document.body.children.length, 0);
assert.equal(timers[0].delay, 5000);
update(snapshot('second', [run('latest.json', 3000, 4), newer, older]));
assert.equal(elements.run.value, '0');
assert.equal(elements.search.value, '');
assert.equal(elements.run.children.length, 3);
assert.match(cell().textContent, /4 \/ 4/);
elements.download.click();
download.text().then(text => {
  const exported = JSON.parse(text);
  assert.equal(exported.run.name, 'latest.json');
  assert.equal(exported.run.results[0].exact_agrees, 1);
  console.log('Dashboard live-update checks passed.');
});
