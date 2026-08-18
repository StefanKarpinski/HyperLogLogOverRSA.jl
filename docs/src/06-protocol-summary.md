# Protocol Summary

Here we summarize the full HyperLogLog Over RSA protocol for anonymously estimating the number of unique clients making various kinds of requests to a server or set of servers. At a high level:

- Clients and servers agree on protocol parameters.
- The server generates a HyperLogLog RSA ring and publishes it, along with a semigenerator and a certificate that the modulus is fingerprint-free.
- The client downloads the ring certificate and checks that the parameters are acceptable and that the certificate is valid.
- The client generates and stores a master key ring element.
- For each request, the client generates and sends a freshly randomized encrypted HyperLogLog value, keyed on the master key and the resource class of the request.
- The server responds to requests and collects request logs, including encrypted HyperLogLog values, which may or may not be valid.
- After request logs have been aggregated, HyperLogLog values can be validated and decoded.
- Once individual HyperLogLog values have been decoded, client count estimation is just a matter of normal HyperLogLog aggregation and estimation.

A test implementation of the entire HyperLogLog Over RSA protocol in Julia can be found [here](https://github.com/StefanKarpinski/HyperLogLogOverRSA.jl).

## Server step 0: Parameter selection

The system operator chooses HyperLogLog parameters, RSA bit-length, and a certificate strength:

- ``B`` is the bucket count
    - It must be an odd number
    - Example: $B = 2^{12} - 1$
- ``m`` is the geometric value cap
    - Ranks run $0$ to $m-1$; a saturated bucket reports $m-1$
    - Example: $m = 64$
- ``L`` is the bit-length of $N$ values
    - Controls cryptographic strength of the RSA ring
    - Example: $L = 1024$
- ``\alpha`` is the certificate strength
    - Effort to generate a valid-looking malicious $N$ value
    - Should be at least $2^\lambda$ for the system security level $\lambda$, and never below $2^{80}$ (see [Malicious servers](05-security-analysis.md#Malicious-servers))
    - Example: $\alpha = 2^{112}$

## Client step 0: Parameter acceptance criteria

On behalf of clients, the client implementor should choose acceptance criteria for protocol parameters, including:

- ``B_{\max}`` — maximum bucket count
    - The simplest way to fingerprint clients is just to choose $B = 2^{128}$ and let the bucket be the fingerprint. This limit prevents that kind of “attack”.
    - Example: $B_{\max} = 2^{12}$
- ``m_{\max}`` — maximum geometric cap
    - Mostly a sanity check: extreme geometric samples are vanishingly rare, and we don’t want a malicious server forcing a client to work in an absurdly large geometric range. The real ceiling is the width of the hash used to derive per-class values, since the geometric coordinate is uniform only while $m$ is at most that width — $128$ bits in the reference derivation.
    - Example: $m_{\max} = 128$
- ``L_{\max}`` — maximum modulus bit-length
    - This is also mostly a sanity check to make sure that clients aren’t DoSed by being made to do arithmetic in some absurdly large modulus.
    - Example: $L_{\max} = 2^{20}$
- ``\alpha_{\min}`` — minimum certificate strength
    - This is the least number of modulus values a malicious server would expect to have to try in order to find one that passes certificate checks. The server can provide more square roots than this strength implies, but not fewer.
    - Example: $\alpha_{\min} = 2^{112}$

## Server step 1: Ring generation

The RSA ring is $\Z_N$ where

```math
\begin{aligned}
N = PQ = (4 B p + 1)(2^m q + 1)
\end{aligned}
```

such that $P$, $Q$, $p$, $q$, are all distinct odd primes coprime to $B$. The factor of $4$ makes $P \equiv 5 \bmod 8$, hence $N \equiv 5 \bmod 8$.

The server generates $N$ such that:

```math
\begin{aligned}
2^{L-1} ≤ N < 2^L
\end{aligned}
```

This can be accomplished by:

- Choose random $p$
- Check that $p$ is prime
- Check that $P = 4 B p + 1$ is prime
- Check that $B \notin \set{P, p}$
- Choose random $q$
- Check that $q$ is prime
- Check that $Q = 2^m q + 1$ is prime
- Check that $\set{B, P, p} \inter \set{Q, q} = \emptyset$

The ranges that $p$ and $q$ are chosen from should be computed such that $N = PQ$ falls into $[2^{L-1}, 2^L)$ as required. It’s also desirable for $P$ and $Q$ to have $\log_2(P) \approx L/2 \approx \log_2(Q)$ for some notion of approximation.

## Server step 2: Semigenerator selection

The server chooses a random semigenerator, $g \in \Z_N^*$.

This can be accomplished by:

- Choose random $g \in \Z_N$
- Check that $\fmod(g, P)$ is a generator in $\Z_P$
- Check that $\fmod(g, Q)$ is a generator in $\Z_Q$

## Server step 3: Square root computation

The server computes the number of certificate square roots required based on the target certificate strength, $\alpha$:

```math
\begin{aligned}
n = \ceil{\frac{\log_2(\alpha)}{\log_2(8)-\log_2(5)}} 
\end{aligned}
```

The server generates $n$ pairs of values, $\set{x, y} \subseteq J_N^+$, using an agreed-upon hashing scheme that includes the value $N$ as a hash input. The hashed values must land in $J_N^+$ (positive Jacobi symbol): since a hash naturally produces values of either Jacobi symbol, the scheme multiplies any hashed value of negative Jacobi symbol by $2$, which lies in $J_N^-$ because $N \equiv 5 \bmod 8$, mapping it into $J_N^+$. (A hash with Jacobi symbol $0$ shares a factor with $N$, exposing $N$ as factorable; the server rejects such an $N$ outright.) For each pair, $\set{x, y}$:

- If $x$ is a quadratic residue, add $r$ such that $r^2 = x \bmod N$ to the list of square roots
- Otherwise, if $y$ is a quadratic residue, add $r$ such that $r^2 = y \bmod N$ to the list of square roots
- Otherwise, if $xy$ is a quadratic residue, add $r$ such that $r^2 = xy \bmod N$ to the list of square roots

Since $N$ is a semiprime and $x, y \in J_N^+$, one of these three checks must succeed (see [Malicious servers](05-security-analysis.md#Malicious-servers)), so one square root value is added to the list for each pair generated.

## Server step 4: Ring certificate publication

The server publishes a ring certificate containing:

- ``B`` — the bucket count
- ``m`` — the geometric value cap
- ``N`` — the RSA ring modulus
- ``g`` — the ring semigenerator value
- ``\text{sqrts}`` — a list of square roots

## Client step 1: Ring certificate checking

The client downloads the latest ring certificate at an agreed upon location. It may or may not have a previously downloaded and verified ring certificate. If it has a downloaded ring and the new one is the same, then it can simply use the existing ring and twist element. If it doesn’t have an existing ring certificate or the new one is different, then the client should delete the old ring certificate and proceed with checking and saving the new certificate and twist element.

To check a ring certificate, the client should:

- Check that $B ≤ B_{\max}$
- Check that $B$ is odd
- Check that $3 ≤ m ≤ m_{\max}$
- Check that $\log_2(N) ≤ L_{\max}$
- Check that $N = 5 \bmod 8$
- Check that $\gcd(B, N) = 1$
- Check that $\gcd(B, N-1) = 1$
- Check that $\Jacobi_N(g) = 1$
- Check that $\alpha_{\min} ≤ (8/5)^n$ where $n$ is the length of $\text{sqrts}$
- Check that each square root is valid

If a certificate passes checks, the client should generate a random twist element:

- Choose random $x_0 \in \Z_N$
- Check that $\Jacobi_N(x_0) = -1$

Replace any existing ring record with a new ring record:

- ``B``, $m$, $N$, $g$, $x_0$

This must be stored persistently so that the client uses the same $x_0$ in different sessions. The client should check for a new ring certificate periodically.

## Client step 2: Request generation

When the client needs to send a request in a resource class, it computes:

```math
\begin{aligned}
x = x_0 g^h = x_0 g^{\hash(x_0,\,\text{class})}
\end{aligned}
```

It then generates a random white noise element and a random exponent value:

- Choose $z \in \set{1, 2, \dots, N-1}$ and a random sign, and let $w = \pm z^{B2^{m-1}} \bmod N$
- Choose $i \in \set{0, 1, \dots, 2^{m-1}-1}$ and let $t = 2Bi + 1$

Finally, the client computes:

```math
\begin{aligned}
y = \fmod(wx^t, N)
\end{aligned}
```

This value $y$ is what the client sends with the request, in a header that also carries a short identifier for the ring — a truncated hash of the ring parameters — so the server knows which modulus $y$ belongs to. A few bytes are enough to tell a server’s rings apart, so the identifier stays far smaller than $y$ itself.

## Server step 5: Request validation & decoding

The server should not attempt to validate request headers while responding to requests, it should simply serve requests and log the header information for later processing. This also means that the factorization of $N$ does not need to reside on public-facing servers—it only needs to be available for offline log processing. This significantly reduces the chances of the factorization being accidentally leaked or exfiltrated by an attacker.

When post-processing request logs, the server can validate the submitted $y$ values by checking that the $N$ value is known and that $\Jacobi_N(y) = -1$. Any request with an unknown ring modulus or that fails the Jacobi symbol test should be ignored for client counting purposes.

The server decodes the HLL value by computing:

```math
\begin{aligned}
\text{bucket} &= \fmod(y^{4p}, P) \\
\text{geometric} &= \min\!\big(m-1,\; m - \log_2(\ord(\fmod(y^q, Q)))\big) \\
\end{aligned}
```

The bucket value is in $\Z_P$ which is huge, but there are only $B$ possible values it takes which can be mapped back to $\set{0, \dots, B-1}$ by any consistent mapping scheme. The geometric values land in $\set{0, \dots, m-1}$ and are computed efficiently by repeatedly squaring $\fmod(y^q,Q)$ until reaching one — the number of squarings is $\log_2(\ord)$. The result is then capped at $m-1$.

## Server step 6: Count estimation

Once the bucket and geometric count have been decoded, estimating the unique client count within any set of logs for the same resource class is a matter of doing normal HyperLogLog aggregation and estimation: compute the maximum geometric sample size for each bucket and use the histogram of those per-bucket maxima to estimate the client count. Several good estimators exist. The reference implementation uses the *improved estimator* of Otmar Ertl, [“New Cardinality Estimation Methods for HyperLogLog Sketches”](https://arxiv.org/abs/1706.07290): a closed-form estimator that is effectively unbiased across the whole cardinality range—with no separate small- or large-range corrections—and has relative standard error of about $1.04/\sqrt{B}$. The same paper’s *maximum likelihood estimator* is marginally more accurate (it attains the theoretical variance bound) but requires an iterative solve rather than a closed form. Both supersede the original Flajolet *et al.* estimator, and both are more principled than the empirical bias-correction tables of HyperLogLog++ (Heule *et al.*), which is a separate method and not as accurate. The reported value is capped at the number of aggregated requests, $\min(\text{estimate}, R)$—there cannot be more unique clients than requests. This also makes using reported count estimates as a geometric coordinate oracle sufficiently expensive as to no longer be an effective attack.

A pleasant feature of this protocol is that once HyperLogLog values are validated and decoded, aggregation and estimation can be done by anyone, so long as they only try to aggregate within resource classes. In other words, the “over RSA” part of the protocol can be forgotten about entirely once HyperLogLog values are decoded. At that point, it’s as if the collection of logs had been recorded with normal HyperLogLog values keyed on client ID and resource class for each record—which is also a reminder of what must *not* be done with them, the subject of the next section.

## Reporting and publishing counts safely

Client count estimates are not of any use unless they get published. HyperLogLog allows estimating the number of unique clients in arbitrary slices of requests, but some care should be taken in what slices are actually reported.

It is tempting to think that, because each decoded record is only a small HyperLogLog value, the decoded logs are harmless and could simply be published as raw request-level records. They should not be. In any class small enough that clients are unlikely to share $(b, k)$ pairs, these decoded pairs uniquely identify clients, in much the same way that being a FreeBSD user effectively uniquely identifies someone (we *know* what you’re up to, Alex). Publishing decoded HyperLogLog pairs, would also give an attacker an oracle for decoding their HLL values, which, at the very least, would allow them to learn which pairs are rare.

What *is* safe to expose is **aggregate** output for large enough slices: per-slice cardinality, rather than the per-request $(b, k)$ values themselves. This is safe so long as a **minimum-cardinality floor** is enforced—report no slice whose estimated cardinality falls below a threshold $T$, folding everything beneath it into a single “other” bucket. How large $T$ must be has a clean answer: a published estimate keeps a client anonymous exactly when it cannot reveal whether that client is present, and we can measure how close it comes. Adding or removing one client shifts the expected estimate by $1$, against a standard error of about $1.04\,n/\sqrt{B}$, so a single client’s detectability is the ratio $\sqrt{B}/(1.04\,n)$. Once that ratio drops well below $1$, the estimator’s own sampling noise drowns out any individual contribution. Setting $n = \sqrt{B}$ ($\approx 64$ for $B = 2^{12}$) makes that ratio about $1$, so at a cardinality of $\sqrt{B}$ a lone client is still fully visible; bringing it down to a quarter or an eighth takes only a small multiple more, a floor of four to eight times $\sqrt{B}$, a few hundred clients, with a round $T = 1024$ at the conservative end. The bound is for a single published estimate, though; an adversary who differences overlapping slices or repeats queries averages the noise away, and must be limited separately.

A second, independent rule protects *count integrity* rather than anonymity, and the two should not be conflated: **cap every reported estimate at the number of requests aggregated into the slice**. There cannot be more unique clients in a set of requests than there are requests, so this is correct regardless of any attack—and it is exactly what bounds a malicious client’s inflation. Without it, an attacker who can read fine-grained per-slice counts turns them into a decoding oracle and curates a flat batch that pushes a slice of $B$ requests to an estimate of $\sim \frac{1}{\ln 2} B 2^k$ (see [Malicious clients](05-security-analysis.md#Malicious-clients)); the cap clips that back to the $B$ requests that produced it, so inflation never exceeds the forged request count (plus any honest duplication slack). The floor and the cap are independent: the floor suppresses a small slice the cap would still report, and the cap bounds a large fake estimate the floor would happily publish.

## Deployment in Julia’s Pkg client

The protocol above is general; here we pin down the concrete choices for the Julia Pkg client, the setting these ideas were developed for.

**Acceptance bounds.** A client accepts a ring only within $B_{\max} = 2^{12}$, $m_{\max} = 128$, $L_{\max} = 2^{20}$, and $\alpha_{\min} = 2^{112}$ — so a conforming certificate carries $n ≥ 166$ square roots, all of which the client verifies. The relatively high $m_{\max} = 128$ is possible because the client does its per-request arithmetic in arbitrary-precision integers; the client doesn’t care what kind of integers the server needs to use for its arithmetic.

**Certificate endpoint.** The server publishes its certificate as TOML at `$server/hll_rsa.toml`. Client and server both use Julia’s TOML implementation, which reads and writes arbitrary-precision integers, so $N$, $g$, and the square roots are bare integer literals — no quoting is needed even though they far exceed 64 bits:

```toml
B = 4095
m = 64
N = 1152665851984795538…
g = 2154516298683041933…
sqrts = [3524590212…, 4461971058…]   # n = 166 of them at α = 2^112
```

The client fetches this endpoint with a plain download that does not pass through the metadata-header path, so the request for the certificate never carries an HLL header. That breaks what would otherwise be a circular dependency: needing a verified ring in order to build the header that would be sent while fetching that very ring.

**Stored ring record.** Once a certificate verifies, the client saves the ring — but not the large, no-longer-needed square roots — together with a freshly generated master key $x_0 \in J_N^-$, to `~/.julia/servers/<host>/hll_rsa.toml` (again as bare integers):

```toml
B = 4095
m = 64
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
