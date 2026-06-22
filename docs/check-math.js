#!/usr/bin/env node
//
// check-math.js — render every math block in the built docs through KaTeX and
// report any that fail to parse (or trip a strict-mode warning).
//
// Documenter renders math with KaTeX, which is stricter than the MathJax used
// by Obsidian/HackMD. This catches LaTeX that those tolerate but KaTeX rejects
// — unbraced multi-letter subscripts (`B_\max`), wrong alignat column counts,
// stray non-breaking spaces, unsupported macros, etc. — which otherwise show
// up as raw, unrendered LaTeX in the browser.
//
// Usage:
//   julia --project=docs docs/make.jl     # 1. build the docs
//   node docs/check-math.js               # 2. validate the rendered math
//
// Exit codes: 0 = all good, 1 = one or more blocks failed, 2 = setup problem.
//
// Requires KaTeX, found automatically from either:
//   - npm:  npm install katex
//   - apt:  sudo apt install katex   (provides /usr/share/nodejs/katex)

const fs = require('fs');
const path = require('path');

let katex;
for (const m of ['katex', '/usr/share/nodejs/katex']) {
  try { katex = require(m); break; } catch (_) { /* try next */ }
}
if (!katex) {
  console.error('error: KaTeX not found. Install with `npm install katex` or `sudo apt install katex`.');
  process.exit(2);
}

const docsDir  = __dirname;
const buildDir = path.join(docsDir, 'build');
const makeJl   = path.join(docsDir, 'make.jl');

if (!fs.existsSync(buildDir)) {
  console.error(`error: ${buildDir} not found — build the docs first:\n  julia --project=docs docs/make.jl`);
  process.exit(2);
}

// Macros: parse the `raw"\x" => raw"\y"` pairs straight out of make.jl, so this
// always validates against the exact macro set the site ships with.
const macros = {};
for (const m of fs.readFileSync(makeJl, 'utf8').matchAll(/raw"(\\[^"]*)"\s*=>\s*raw"([^"]*)"/g)) {
  macros[m[1]] = m[2];
}

function unescapeHtml(s) {
  return s.replace(/&lt;/g, '<').replace(/&gt;/g, '>').replace(/&quot;/g, '"')
          .replace(/&#39;/g, "'").replace(/&#(\d+);/g, (_, n) => String.fromCharCode(+n))
          .replace(/&amp;/g, '&'); // ampersand last
}

const files = fs.readdirSync(buildDir).filter(f => f.endsWith('.html'));
let total = 0;
const failures = [];
const warnings = [];

for (const f of files) {
  const html = fs.readFileSync(path.join(buildDir, f), 'utf8');
  const blocks = [];
  // Documenter emits display math as <p class="math-container">\[ … \]</p>
  for (const m of html.matchAll(/<p class="math-container">([\s\S]*?)<\/p>/g)) {
    blocks.push({ display: true, raw: m[1] });
  }
  // …and inline math as <span>$ … $</span>
  for (const m of html.matchAll(/<span>\$([\s\S]*?)\$<\/span>/g)) {
    blocks.push({ display: false, raw: m[1] });
  }

  for (const b of blocks) {
    total++;
    let tex = unescapeHtml(b.raw);
    if (b.display) tex = tex.replace(/^\s*\\\[/, '').replace(/\\\]\s*$/, '');
    const blockWarnings = [];
    try {
      katex.renderToString(tex, {
        displayMode: b.display,
        throwOnError: true,
        macros: Object.assign({}, macros), // fresh copy: renderToString mutates it
        strict: (code, msg) => { blockWarnings.push(`${code}: ${msg}`); return 'ignore'; },
      });
    } catch (e) {
      failures.push({ file: f, display: b.display, tex, err: String(e.message || e) });
    }
    for (const w of blockWarnings) warnings.push({ file: f, tex, w });
  }
}

const trunc = (s, n = 180) => s.replace(/\s+/g, ' ').trim().slice(0, n);

if (total === 0) {
  console.error('error: found no math blocks — the build may be empty or Documenter\'s HTML format may have changed.');
  process.exit(2);
}

console.log(`Checked ${total} math blocks across ${files.length} pages.`);
console.log(`  failures: ${failures.length}    warnings: ${warnings.length}\n`);

for (const x of failures) {
  console.log(`FAIL [${x.file}] (${x.display ? 'display' : 'inline'})`);
  console.log(`  TeX: ${trunc(x.tex)}`);
  console.log(`  ERR: ${trunc(x.err)}\n`);
}
for (const x of warnings) {
  console.log(`warn [${x.file}]: ${trunc(x.w)}`);
  console.log(`  TeX: ${trunc(x.tex)}`);
}

process.exit(failures.length ? 1 : 0);
