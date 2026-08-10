# Deployment in Julia’s Pkg Client

The protocol above is general; here we pin down the concrete choices for the Julia Pkg client, the setting these ideas were developed for.

**Acceptance bounds.** A client accepts a ring only within $B_{\max} = 2^{12}$, $m_{\max} = 127$, $L_{\max} = 2^{20}$, and $\alpha_{\min} = 2^{112}$ — so a conforming certificate carries $n ≥ 166$ square roots, all of which the client verifies. The relatively high $m_{\max} = 127$ is possible because the client does its per-request arithmetic in arbitrary-precision integers; the client doesn’t care what kind of integers the server needs to use for its arithmetic.

**Certificate endpoint.** The server publishes its certificate as TOML at `$server/hll_rsa.toml`. Client and server both use Julia’s TOML implementation, which reads and writes arbitrary-precision integers, so $N$, $g$, and the square roots are bare integer literals — no quoting is needed even though they far exceed 64 bits:

```toml
B = 4095
m = 63
N = 1152665851984795538…
g = 2154516298683041933…
sqrts = [3524590212…, 4461971058…]   # n = 166 of them at α = 2^112
```

The client fetches this endpoint with a plain download that does not pass through the metadata-header path, so the request for the certificate never carries an HLL header. That breaks what would otherwise be a circular dependency: needing a verified ring in order to build the header that would be sent while fetching that very ring.

**Stored ring record.** Once a certificate verifies, the client saves the ring — but not the large, no-longer-needed square roots — together with a freshly generated master key $x_0 \in J_N^-$, to `~/.julia/servers/<host>/hll_rsa.toml` (again as bare integers):

```toml
B = 4095
m = 63
N = 1152665851984795538…
g = 2154516298683041933…
x0 = 1055559624789921343…
```

The file is written atomically — to a temporary file, then renamed into place — and is rewritten only when the ring actually changes, at which point a fresh $x_0$ is generated in the new ring.

**Ring hashing.** The header carries a compact **ring-id** that identifies the modulus. Write $B$, $m$, $N$, and $g$ in decimal — no leading zeros, and, being nonnegative, no signs — and join them in that order with commas (in Julia, `"$B,$m,$N,$g"`); the ring-id is the first four bytes of the SHA-256 of that string, as lowercase hex (eight characters). Change detection does not use the ring-id, instead comparing parameters directly.

**Class hashing.** The per-class exponent $h = \hash(x_0, \text{class})$ is instantiated as HMAC-SHA-256 keyed by the SHA-256 hash of $x_0$ (its little-endian bytes), applied to the class string’s bytes, with the first 16 bytes of the digest read big-endian into a 128-bit integer. The ring holder decodes with the factorization and never recomputes $h$, so this is purely a client-side choice: any keyed hash of the class works, and keying it by $x_0$ is what makes each client’s per-class draw independent and impossible for the client to bias.

**Refresh.** The client re-fetches the certificate once at the start of a session and once a day thereafter, keeping the last-check time in memory (this timestamp does not need to be saved to a file, only in memory). The client detects a changed ring by comparing its parameters against those of the stored ring; on a change it re-verifies, regenerates $x_0$, and replaces the record. Any failure — no endpoint, a certificate that fails verification, a network error — is non-fatal: the client falls back to its stored ring if it has one, and otherwise simply sends no header.

**Header.** On every request to the package server, an enabled client sends a `Julia-Pkg-HLL-RSA` header whose value takes one of three forms:

```
<ring-id>,<token>    # request in a resource class: the HLL token
noclass              # request in no counted resource class
error,<reason>       # the token could not be produced (see below)
```

For a counted request, `<ring-id>` is the ring-id defined above — enough for the server to look up the modulus, and kept short so it costs little next to the token — and `<token>` is the base64 encoding of $y$ written as big-endian bytes.

A request that maps to *no* resource class (server metadata, for instance) can carry no token, so the client sends the fixed value `noclass`. This is a normal, explicitly-handled outcome, not an error. Reporting it explicitly — rather than omitting the header — is deliberate: it makes an *absent* header mean one thing only, that the client has opted out.

**Error reporting.** Anything in the client’s HLL machinery can fail on a request that would otherwise carry a token — the certificate can’t be fetched or parsed, it fails verification, or something unexpected goes wrong. Rather than fall silent, the client then sends `error,<reason>`, with `<reason>` always one of exactly three fixed strings — `fetch`, `verify`, or `internal` — never anything derived from client state, so it leaks nothing about the client. `internal` means an unexpected client-side failure specifically, distinct from the `noclass` status above. This never disrupts the request, but it lets the server distinguish a client that hit an error from one that has opted out (which sends no header at all) — useful for spotting a broken certificate or a client-side bug across the fleet.

**Opt-out.** The header is sent by default. Setting `JULIA_PKG_SERVER_HLL_RSA` to a false value (`0`, `false`, `no`, `f`, …) disables it entirely — the one and only case in which no header is sent at all.
