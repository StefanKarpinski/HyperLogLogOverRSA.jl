# Deployment in Julia’s Pkg Client

The protocol above is general; here we pin down the concrete choices for the Julia Pkg client, the setting these ideas were developed for.

**Acceptance bounds.** A client accepts a ring only within $B_{\max} = 2^{12}$, $m_{\max} = 127$, $L_{\max} = 2^{20}$, and $\alpha_{\min} = 2^{112}$ — so a conforming certificate carries $n ≥ 166$ square roots, all of which the client verifies. The relatively high $m_{\max} = 127$ is possible because the client does its per-request arithmetic in arbitrary-precision integers; the client doesn’t care what kind of integers the server needs to use for its arithmetic.

**Certificate endpoint.** The server publishes its certificate as TOML at `$server/hll_rsa.toml`. Client and server both use Julia’s TOML implementation, which reads and writes arbitrary-precision integers, so $N$, $g$, and the square roots are bare integer literals — no quoting is needed even though they far exceed 64 bits:

```toml
B = 4094
m = 63
N = 1152665851984795538…
g = 2154516298683041933…
sqrts = [3524590212…, 4461971058…]   # n = 166 of them at α = 2^112
```

The client fetches this endpoint with a plain download that does not pass through the metadata-header path, so the request for the certificate never carries an HLL header. That breaks what would otherwise be a circular dependency: needing a verified ring in order to build the header that would be sent while fetching that very ring.

**Stored ring record.** Once a certificate verifies, the client saves the ring — but not the large, no-longer-needed square roots — together with a freshly generated master key $x_0 \in J_N^-$, to `~/.julia/servers/<host>/hll_rsa.toml` (again as bare integers):

```toml
B = 4094
m = 63
N = 1152665851984795538…
g = 2154516298683041933…
x0 = 1055559624789921343…
expiry = 1793491200                   # unix time of the next x₀ rotation
```

The file is written atomically — to a temporary file, then renamed into place — and is rewritten when the ring changes (a fresh $x_0$ is generated in the new ring) or when the master key reaches its expiry and is rotated in place (see **Key regeneration** below).

**Ring hashing.** The header carries a compact **ring-id** that identifies the modulus. Write $B$, $m$, $N$, and $g$ in decimal — no leading zeros, and, being nonnegative, no signs — and join them in that order with commas (in Julia, `"$B,$m,$N,$g"`); the ring-id is the first four bytes of the SHA-256 of that string, as lowercase hex (eight characters). Change detection does not use the ring-id, instead comparing parameters directly.

**Class hashing.** The per-class exponent $h = \hash(x_0, \text{class})$ is instantiated as HMAC-SHA-256 keyed by the SHA-256 hash of $x_0$ (its little-endian bytes), applied to the class string’s bytes, with the first 16 bytes of the digest read big-endian into a 128-bit integer. The ring holder decodes with the factorization and never recomputes $h$, so this is purely a client-side choice: any keyed hash of the class works, and keying it by $x_0$ is what makes each client’s per-class draw independent and impossible for the client to bias.

**Semisharding.** The Pkg client uses the [semisharded](04-hyperloglog-over-rsa.md#Semisharding) construction, deriving each class’s element with $f = g^{B/2}$ so a client keeps a common bucket pair across all resource classes. Pkg requests are the case the variant exists for: a project’s dependency tree is fetched as a bundle and then re-fetched, resolved, and CI-tested repeatedly, so a client’s request bundles recur with overlap over time. The [request-bundle analysis](05-security-analysis.md#Request-bundles-and-cross-class-correlation) explains why that turns the plain-sharded bucket into a cross-class fingerprint, and why fixing it is the safe default.

**Deterministic semigenerator.** The certificate’s $g$ is not a free server choice but a deterministic function of $N$ — the hash of $N$ mapped into $J_N^+$ by the ring’s fixed twist — so the client recomputes it and rejects any certificate whose $g$ differs. This denies a malicious server the use of $g$ as a per-client tag; half of candidate moduli admit a usable $g$ at the first hash index, so a server regenerates $N$ in the rare case its first candidate fails, and the client’s local check remains the same $\Jacobi_N(g) = 1$ it already performs (see [Malicious semigenerators](05-security-analysis.md#Malicious-semigenerators)).

**Refresh.** The client re-fetches the certificate once at the start of a session and once a day thereafter, keeping the last-check time in memory (this timestamp does not need to be saved to a file, only in memory). The client detects a changed ring by comparing its parameters against those of the stored ring; on a change it re-verifies, regenerates $x_0$, and replaces the record. Any failure — no endpoint, a certificate that fails verification, a network error — is non-fatal: the client falls back to its stored ring if it has one, and otherwise simply sends no header.

**Key regeneration.** The stored record carries an expiry timestamp. When the client first generates $x_0$ it sets the expiry to a time drawn uniformly within the next $T$ (one year); on each daily refresh, if the current time is past the expiry, it generates a fresh $x_0$ in the same ring and advances the expiry by whole periods of $T$ — so a client returning from a long dormancy rotates once and lands back on its own uniformly-random phase rather than resetting to the moment it noticed. Rotating $x_0$ yearly bounds long-term linkability to a single period; the reasoning, and the survival model that lets counts be corrected for the induced churn, are in [Key regeneration and long-term unlinkability](07-reporting-and-retention.md#Key-regeneration-and-long-term-unlinkability).

**Header.** On every request to the package server, an enabled client sends a `Julia-Pkg-HLL-RSA` header whose value takes one of three forms:

```
<ring-id>,<token>    # request in a resource class: the HLL token
noclass              # request in no counted resource class
error,<reason>       # the token could not be produced (see below)
```

For a counted request, `<ring-id>` is the ring-id defined above — enough for the server to look up the modulus, and kept short so it costs little next to the token — and `<token>` is the base64 encoding of $y$ written as big-endian bytes.

A request that maps to *no* resource class (server metadata, for instance) can carry no token, so the client sends the fixed value `noclass`. This is a normal, explicitly-handled outcome, not an error. Reporting it explicitly — rather than omitting the header — is deliberate: it makes an *absent* header mean one thing only, that the client has opted out.

**Error reporting.** Anything in the client’s HLL machinery can fail on a request that would otherwise carry a token — the certificate can’t be fetched or parsed, it fails verification, or something unexpected goes wrong. Rather than fall silent, the client then sends `error,<reason>`, with `<reason>` always one of exactly three fixed strings — `fetch`, `verify`, or `internal` — never anything derived from client state, so it leaks nothing about the client. `internal` means an unexpected client-side failure specifically, distinct from the `noclass` status above. This never disrupts the request, but it lets the server distinguish a client that hit an error from one that has opted out (which sends no header at all) — useful for spotting a broken certificate or a client-side bug across the fleet.

**CI and ephemeral clients.** The client detects continuous-integration and other ephemeral environments from environment variables (`CI`, `GITHUB_ACTIONS`, …) and flags them, so that automated runs — which regenerate a master key on essentially every invocation — can be separated from persistent human clients rather than inflating the counts (see [Ephemeral clients](07-reporting-and-retention.md#Ephemeral-clients)).

**Opt-out.** The header is sent by default. Setting `JULIA_PKG_SERVER_HLL_RSA` to a false value (`0`, `false`, `no`, `f`, …) disables it entirely — the one and only case in which no header is sent at all.
