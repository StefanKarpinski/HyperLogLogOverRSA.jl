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
    - It must be congruent to $2 \bmod 4$
    - Example: $B = 2^{12} - 2$
- ``m`` is the height of the geometric ladder
    - Reported geometric samples are capped at $m-1$
    - Maximum client estimate is around $2^{m-1}$
    - Example: $m = 63$
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
- ``m_{\max}`` — maximum geometric sample
    - Mostly a sanity check: extreme geometric samples are vanishingly rare, and we don’t want a malicious server forcing a client to work in an absurdly large geometric range. The real ceiling is the width of the hash used to derive per-class values, since the geometric coordinate is uniform only while $m$ is at most that width — $128$ bits in the reference derivation. Within that, a client that does its per-request exponent arithmetic in fixed-width integers should cap near $63$ (to stay in `Int64`/`Int128`), while one using arbitrary-precision integers can accept up to $127$.
    - Example: $m_{\max} = 127$ (arbitrary-precision client) or $63$ (fixed-width client)
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
N = PQ = (2 B p + 1)(2^m q + 1)
\end{aligned}
```

such that $P$, $Q$, $p$, $q$, are all distinct odd primes coprime to $\nfrac{B}{2}$.

The server generates $N$ such that:

```math
\begin{aligned}
2^{L-1} ≤ N < 2^L
\end{aligned}
```

This can be accomplished by:

- Choose random $p$
- Check that $p$ is prime
- Check that $P = 2 B p + 1$ is prime
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

The server generates $n$ pairs of values, $\set{x, y} \subseteq J_N^+$, using an agreed-upon hashing scheme that includes the value $N$ as a hash input. The hashed values must land in $J_N^+$ (positive Jacobi symbol): since a hash naturally produces values of either Jacobi symbol, the scheme fixes a twist element $\tau$ with $\Jacobi_N(\tau) = -1$ — chosen deterministically, e.g. the smallest such value — and multiplies any hashed value of negative Jacobi symbol by $\tau$, which maps it into $J_N^+$. For each pair, $\set{x, y}$:

- If $x$ is a quadratic residue, add $r$ such that $r^2 = x \bmod N$ to the list of square roots
- Otherwise, if $y$ is a quadratic residue, add $r$ such that $r^2 = y \bmod N$ to the list of square roots
- Otherwise, if $xy$ is a quadratic residue, add $r$ such that $r^2 = xy \bmod N$ to the list of square roots

Since $N$ is a semiprime and $x, y \in J_N^+$, one of these three checks must succeed (see [Malicious servers](05-security-analysis.md#Malicious-servers)), so one square root value is added to the list for each pair generated.

## Server step 4: Ring certificate publication

The server publishes a ring certificate containing:

- ``B`` — the bucket count
- ``m`` — the maximum geometric sample
- ``N`` — the RSA ring modulus
- ``g`` — the ring semigenerator value
- ``\text{sqrts}`` — a list of square roots

## Client step 1: Ring certificate checking

The client downloads the latest ring certificate at an agreed upon location. It may or may not have a previously downloaded and verified ring certificate. If it has a downloaded ring and the new one is the same, then it can simply use the existing ring and twist element. If it doesn’t have an existing ring certificate or the new one is different, then the client should delete the old ring certificate and proceed with checking and saving the new certificate and twist element.

To check a ring certificate, the client should:

- Check that $B ≤ B_{\max}$
- Check that $B = 2 \bmod 4$
- Check that $2 ≤ m ≤ m_{\max}$
- Check that $\log_2(N) ≤ L_{\max}$
- Check that $N = 5 \bmod 8$
- Check that $\gcd(B, N) = 1$
- Check that $\gcd(B, N-1) = 2$
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

- Choose $z \in \set{1, 2, \dots, N-1}$ and let $w = z^{B2^{m-1}} \bmod N$
- Choose $i \in \set{0, 1, \dots, 2^{m-1}-1}$ and let $t = Bi + 1$

Finally, the client computes:

```math
\begin{aligned}
y = \fmod(wx^t, N)
\end{aligned}
```

This value $y$ is what the client sends with the request, in a header that also carries a short identifier for the ring — a truncated hash of the ring parameters — so the server knows which modulus $y$ belongs to. A few bytes are enough to tell a server’s rings apart, so the identifier stays far smaller than $y$ itself.

To keep a client’s bucket *common* across resource classes rather than independent per class — the *semisharding* variant — the client derives $x$ with $f = g^{B/2}$ in place of $g$, so that $x = x_0 f^h = x_0 g^{(B/2)h}$. Because $g^{B/2}$ has a trivial $C_{B/2}$ component, the bucket becomes the client’s fixed master-key bucket while the geometric sample still varies per class. This defends against the cross-class correlation attack on clients that make consistent request bundles over time ([Request bundles and cross-class correlation](05-security-analysis.md#Request-bundles-and-cross-class-correlation)); it is the variant Julia’s Pkg client ships.

## Server step 5: Request validation & decoding

The server should not attempt to validate request headers while responding to requests, it should simply serve requests and log the header information for later processing. This also means that the factorization of $N$ does not need to reside on public-facing servers—it only needs to be available for offline log processing. This significantly reduces the chances of the factorization being accidentally leaked or exfiltrated by an attacker.

When post-processing request logs, the server can validate the submitted $y$ values by checking that the $N$ value is known and that $\Jacobi_N(y) = -1$. Any request with an unknown ring modulus or that fails the Jacobi symbol test should be ignored for client counting purposes.

The server decodes the HLL value by computing:

```math
\begin{aligned}
a &= \log(\fmod(y^p, P)) \\
k &= m - \log_2(\ord(\fmod(y^q, Q))) \\
\end{aligned}
```

Raising to $p$ annihilates the $C_p$ component and leaves the $C_{2B}$ one, whose logarithm $a$ is recovered by table lookup: the value is in $\Z_P$, which is huge, but it takes only $2B$ possible values, which any consistent mapping scheme sends to $\set{0, \dots, 2B-1}$. The order in the second line is computed efficiently by repeatedly squaring $\fmod(y^q,Q)$ until reaching one. The reported pair is then

```math
\begin{aligned}
\text{geometric} &= \min(k,\, m-1) &&\in \set{0, \dots, m-1} \\
\text{bucket} &= \left(a \bmod \tfrac{B}{2}\right) + \tfrac{B}{2}s &&\in \set{0, \dots, B-1}
\end{aligned}
```

where $s$ is the one bit of $2$-torsion that survives re-randomization. For $k ≤ m-2$ it is $1$ exactly when $(a \bmod 4)\,u = 2, 3 \bmod 4$, with $u$ read off the order-$4$ element the squaring chain passes through on its way to one; at the capped level it is $1$ exactly when $k = m$. Capping at $m-1$ is what makes the decoded pair range over an exact $B \times m$ rectangle — see [the anonymity theorem](05-security-analysis.md#The-anonymity-theorem).

## Server step 6: Count estimation

Once the bucket and geometric count have been decoded, estimating the unique client count within any set of logs for the same resource class is a matter of doing normal HyperLogLog aggregation and estimation: compute the maximum geometric sample size for each bucket and use the histogram of those per-bucket maxima to estimate the client count. Several good estimators exist. The reference implementation uses the *improved estimator* of Otmar Ertl, [“New Cardinality Estimation Methods for HyperLogLog Sketches”](https://arxiv.org/abs/1706.07290): a closed-form estimator that is effectively unbiased across the whole cardinality range—with no separate small- or large-range corrections—and has relative standard error of about $1.04/\sqrt{B}$. The same paper’s *maximum likelihood estimator* is marginally more accurate (it attains the theoretical variance bound) but requires an iterative solve rather than a closed form. Both supersede the original Flajolet *et al.* estimator, and both are more principled than the empirical bias-correction tables of HyperLogLog++ (Heule *et al.*), which is a separate method and not as accurate.

A pleasant feature of this protocol is that once HyperLogLog values are validated and decoded, aggregation and estimation can be done by anyone, so long as they only try to aggregate within resource classes. In other words, the “over RSA” part of the protocol can be forgotten about entirely once HyperLogLog values are decoded. At that point, it’s as if the collection of logs had been recorded with normal HyperLogLog values keyed on client ID and resource class for each record—which is also a reminder of what must *not* be done with them, the subject of [Reporting, Publishing & Retention](07-reporting-and-retention.md).
