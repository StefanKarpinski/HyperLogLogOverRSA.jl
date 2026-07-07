# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

`HyperLogLogOverRSA.jl` is **two things in one package**:

1. A **reference/demo implementation** (`src/*.jl`) of the HyperLogLog-Over-RSA
   protocol — counting unique clients of a service such that no one, even with
   full server-log access, can identify or track an individual client. It is a
   correctness/clarity demo, *not* a hardened production library.
2. A **[Documenter.jl](https://documenter.juliadocs.org/) site** (`docs/`) that
   is both the API docs and the full prose writeup of the protocol's
   motivation, design, and security arguments. The writeup is the primary
   deliverable; most ongoing work is on `docs/src/*.md`.

The package and the writeup are kept in sync: `docs/src/reference.md` pulls
docstrings from the exported API.

## Commands

```sh
# Run the test suite
julia --project=. -e 'using Pkg; Pkg.test()'

# Normalize smart quotes in docs prose (run after every docs edit)
python3 docs/normalize-quotes.py

# Build the docs (output → docs/build/, gitignored)
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'  # first time only
julia --project=docs docs/make.jl

# Validate the math AFTER building (renders every math block through KaTeX)
node docs/check-math.js
```

- **After every docs edit:** run `python3 docs/normalize-quotes.py` before
  committing. It converts straight ASCII apostrophes/quotes in prose to
  typographic ("smart") quotes, and — critically — reverts any curly apostrophes
  that crept into math contexts (inline `$…$` or fenced ` ```math `) back to
  ASCII, where KaTeX requires them. Skipping this step is how curly-quote
  corruption in math blocks happens.
- **`julia` hangs on a juliaup config lock?** Call the binary directly:
  `~/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia`.
- **Running a single testset:** `test/runtests.jl` is one file of `@testset`
  blocks. Lines `# false &&` before a testset are manual toggles — uncomment to
  `false && @testset ...` and skip an expensive block while iterating. There is
  no built-in per-test selector; `Pkg.test()` runs the whole file.
- **`check-math.js` needs the build first** and finds KaTeX from either
  `npm install katex` or the apt `katex` package (`/usr/share/nodejs/katex`).
  It exits non-zero on any failure, so it can gate CI. Run it after every docs
  edit — it catches KaTeX-vs-MathJax breakage before it ships as raw LaTeX.

## Commits

- Make **coherent, logical commits** — each commit one self-contained change,
  with a useful and apt message that explains the *what* and *why*.
- **Standard git formatting:** a summary line wrapped at **80 chars**, then a
  blank line, then a body wrapped at **80 chars**.
- **Never push** unless explicitly asked to.
- **Co-authorship:** when *you* authored the change, end the commit message with
  a `Co-Authored-By:` trailer for Claude. When the change was authored by Stefan
  and you are only committing it on his behalf, leave the co-author line off.

## Implementation architecture (`src/`)

The whole scheme rides on an RSA modulus with deliberately engineered structure:

```
N = P·Q = (2·B·p + 1)·(2^m·q + 1)   ⟹   ℤ_N* ≅ C_2 × C_B × C_(2^m) × C_(p·q)
```

- `B` (odd) = number of HyperLogLog **buckets**; `m` = max **geometric** sample.
- The factorization `(P, Q, p, q)` is the **secret** held only by the ring owner
  (the offline decoder). `Ring.show` deliberately omits the primes — **never add
  code that prints or logs them.**

The pipeline (each stage is one struct + constructor):

1. **`Ring(B, m, L)`** (`Ring.jl`) — generates the secret ring at bit-length `L`.
   `N`, `P`, `Q`, `λ` are computed lazily from `(p, q, B, m)` via `getproperty`.
   `ring_type(L)` picks `Int64`/`Int128`/`BigInt` by size.
2. **`RingCert(ring)`** (`RingCert.jl`) — a *publishable* certificate: `N`, a
   semigenerator `g`, and square roots of hashed elements that prove `N` is
   "fingerprint-free" (has the intended two-prime shape) **without revealing the
   factorization**. `SQRT_SAMPLES` is sized so soundness ≥ `α_min = 2^50`.
3. **`Client(cert)`** (`Client.jl`) — verifies the cert (shape bounds, `N ≡ 3
   mod 4`, gcd conditions, `jacobi(g,N)==1`, enough valid sqrts) and picks a
   random secret twist `x₀` (`jacobi(x₀,N) == -1`).
4. **`hll_generate(client, class)`** — maps `(x₀, class)` to a *stable* element
   `x = x₀·g^h` where `h = H(x₀, class)`, then **re-randomizes** it into a token
   `y = w·xᵗ` with fresh randomness per call. Two tokens from the same
   client+class are unlinkable, yet all decode to the same value. The client
   cannot decode or bias its own sampled value.
5. **`hll_decode(ring, y)`** — uses the secret factorization to recover
   `(bucket, geometric)`: bucket via a discrete log in the `C_B` part mod `P`
   (precompute with `bucket_map`), geometric via the order in the `C_(2^m)` part
   mod `Q`.

Supporting modules: `PrimePairs.jl` (primes of the form `scale·p + 1`),
`Jacobi.jl` (Jacobi symbol, random twist, `mod4`/`mod8` with BigInt fast paths),
`Hashing.jl` (`hash_into_ring` → element of ℤ_N via SHA-512; `hash_resource_class` →
UInt128 via SHA-256). `RingCert.jl` also holds `modsqrt` (Tonelli–Shanks plus
`p≡3 mod 4` / `p≡5 mod 8` fast paths).

Per-class independence: because `h = H(x₀, class)` salts the class with the
client's secret, a client's value is an independent draw per class — values
across classes are uncorrelated by construction.

## Docs structure & conventions (`docs/`)

- **`docs/make.jl` is the source of truth** for page order and titles, not the
  filenames. Filenames are numbered `01`–`06`. Don't renumber files to match
  a reordering — edit the `pages` list in `make.jl` instead.
- **LaTeX macros** used across the writeup are registered with KaTeX in the
  `MATH_MACROS` dict in `make.jl` (mirrors the original `preamble.sty`). Add new
  macros there, not per-page.
- **KaTeX is stricter than the MathJax** that Obsidian/HackMD use. Codified
  gotchas (see `check-math.js`):
  - Do **not** override `\not` — KaTeX builds `≠`/`∉` on top of it; redefining
    it renders them as "¬=" / "¬∈". Use `\lnot`/`\neg` for negation.
  - Brace multi-letter operator subscripts: `B_{\max}`, not `B_\max`.
  - Display math is ` ```math ` fences with `aligned` (not bare `align`).
  - Inline math at the **start** of a list item or paragraph must use
    double-backtick `` ``…`` `` — Documenter promotes block-initial `$…$` to a
    display block.
  - Tables are centered via `docs/src/assets/custom.css`.

## Publish flow

`.github/workflows/Documentation.yml` builds and deploys to GitHub Pages on push
to `main` and on tags (`deploydocs` only publishes from `main`/tags). The live
site is custom-domained at <https://karpinski.org/HyperLogLogOverRSA.jl/>. With
no tagged release yet, docs deploy under `/dev/`. Structural changes are
reviewed before they go live; tag a release (e.g. `v0.1.0`) to publish a
versioned `/stable/` URL.
