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

**Server step 0: Parameter selection**

The system operator chooses HyperLogLog parameters, RSA bit-length, and a certificate strength:

- ``B`` is the bucket count
    - It must be an odd number
    - Example: $B = 2^{12} - 1$
- ``m`` is the maximum geometric sample value
    - Maximum client estimate is around $2^m$
    - Example: $m = 63$
- ``L`` is the bit-length of $N$ values
    - Controls cryptographic strength of the RSA ring
    - Example: $L = 1024$
- ``\alpha`` is the certificate strength
    - Effort to generate a valid-looking malicious $N$ value
    - Should be at least $2^\lambda$ for the system security level $\lambda$, and never below $2^{80}$ (see [Malicious Servers](12-malicious-servers.md))
    - Example: $\alpha = 2^{128}$

**Client step 0: Parameter acceptance criteria**

On behalf of clients, the client implementor should choose acceptance criteria for protocol parameters, including:

- ``B_{\max}`` — maximum bucket count
    - The simplest way to fingerprint clients is just to choose $B = 2^{128}$ and let the bucket be the fingerprint. This limit prevents that kind of “attack”.
    - Example: $B_{\max} = 2^{16}$
- ``m_{\max}`` — maximum geometric sample
    - This is mostly a sanity check since extreme geometric samples are so rare as to not matter. We just want to prevent a malicious server from being able to DoS a client by having it compute enormous integer values.
    - Example: $m_{\max} = 128$
- ``L_{\max}`` — maximum modulus bit-length
    - This is also mostly a sanity check to make sure that clients aren’t DoSed by being made to do arithmetic in some absurdly large modulus.
    - Example: $L_{\max} = 2^{20}$
- ``\alpha_{\min}`` — minimum certificate strength
    - This is the least number of modulus values a malicious server would expect to have to try in order to find one that passes certificate checks. The server can provide more square roots than this strength implies, but not fewer.
    - Example: $\alpha_{\min} = 2^{80}$

**Server step 1: Ring generation**

The RSA ring is $\Z_N$ where

```math
\begin{aligned}
N = PQ = (2 B p + 1)(2^m q + 1)
\end{aligned}
```

such that $P$, $Q$, $p$, $q$, are all distinct odd primes coprime to $B$.

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

**Server step 2: Semigenerator selection**

The server chooses a random semigenerator, $g \in \Z_N^*$.

This can be accomplished by:

- Choose random $g \in \Z_N$
- Check that $\fmod(g, P)$ is a generator in $\Z_P$
- Check that $\fmod(g, Q)$ is a generator in $\Z_Q$

**Server step 3: Square root computation**

The server computes the number of certificate square roots required based on the target certificate strength, $\alpha$:

```math
\begin{aligned}
n = \ceil{\frac{\log_2(\alpha)}{\log_2(8)-\log_2(5)}} 
\end{aligned}
```

The server generates $n$ pairs of values, $\set{x, y} \subseteq \Z_N$, using an agreed-upon hashing scheme that includes the value $N$ as a hash input. For each pair, $\set{x, y}$:

- If $x$ is a quadratic residue, add $r$ such that $r^2 = x \bmod N$ to the list of square roots
- Otherwise, if $y$ is a quadratic residue, add $r$ such that $r^2 = y \bmod N$ to the list of square roots
- Otherwise, if $xy$ is a quadratic residue, add $r$ such that $r^2 = xy \bmod N$ to the list of square roots

Since $N$ is a semiprime, one of these three checks must succeed, so one square root value is added to the list for each pair generated.

**Server step 4: Ring certificate publication**

The server publishes a ring certificate containing:

- ``B`` — the bucket count
- ``m`` — the maximum geometric sample
- ``N`` — the RSA ring modulus
- ``g`` — the ring semigenerator value
- ``\text{sqrts}`` — a list of square roots

**Client step 1: Ring certificate checking**

The client downloads the latest ring certificate at an agreed upon location. It may or may not have a previously downloaded and verified ring certificate. If it has a downloaded ring and the new one is the same, then it can simply use the existing ring and twist element. If it doesn’t have an existing ring certificate or the new one is different, then the client should delete the old ring certificate and proceed with checking and saving the new certificate and twist element.

To check a ring certificate, the client should:

- Check that $B ≤ B_{\max}$
- Check that $B$ is odd
- Check that $2 ≤ m ≤ m_{\max}$
- Check that $\log_2(N) ≤ L_{\max}$
- Check that $N = 3 \bmod 4$
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

**Client step 2: Request generation**

When the client needs to send a request in a resource class, it computes:

```math
\begin{aligned}
x = x_0 g^h = x_0 g^{\hash(x_0,\,\text{class})}
\end{aligned}
```

It then generates a random white noise element and a random exponent value:

- Choose $z \in \set{1, 2, \dots, N-1}$ and let $w = z^{B2^m} \bmod N$
- Choose $i \in \set{0, 1, \dots, 2^{m-1}-1}$ and let $t = 2Bi + 1$

Finally, the client computes:

```math
\begin{aligned}
y = \fmod(wx^t, N)
\end{aligned}
```

This is the value it sends along with the request in a header. The request should also include some hash of the ring parameters so that the server can know which ring the value belongs to. This could just be the modulus, or a cryptographic hash of the ring certificate. Note that $N$ will be large enough that a cryptographic hash is actually the more compact option.

**Server step 5: Request validation & decoding**

The server should not attempt to validate request headers while responding to requests, it should simply serve requests and log the header information for later processing. This also means that the factorization of $N$ does not need to reside on public-facing servers—it only needs to be available for offline log processing. This significantly reduces the chances of the factorization being accidentally leaked or exfiltrated by an attacker.

When post-processing request logs, the server can validate the submitted $y$ values by checking that the $N$ value is known and that $\Jacobi_N(y) = -1$. Any request with an unknown ring modulus or that fails the Jacobi symbol test should be ignored for client counting purposes.

The server decodes the HLL value by computing:

```math
\begin{aligned}
\text{bucket} &= \fmod(y^{2p}, P) \\
\text{geometric} &= m - \log_2(\ord(\fmod(y^q, Q))) \\
\end{aligned}
```

The bucket value is in $\Z_P$ which is huge, but there are only $B$ possible values it takes which can be mapped back to $\set{0, \dots, B-1}$ by any consistent mapping scheme. The geometric values are already in $\set{0, \dots, m}$ and the value can be computed efficiently by repeatedly squaring $\fmod(y^q,Q)$ until reaching one.

**Server step 6: Count estimation**

Once the bucket and geometric count have been decoded, estimating the unique client count within any set of logs for the same resource class is a matter of doing normal HyperLogLog aggregation and estimation: compute the maximum geometric sample size for each bucket and use the histogram of maximum sample counts to estimate the most likely client count. The optimal maximum likelihood estimator was derived in [“New Cardinality Estimation Methods for HyperLogLog Sketches”](https://arxiv.org/abs/1706.07290). This is the recommended estimator—it supersedes the original estimator and the improved estimator that is sometimes called “HyperLogLog++.”

A pleasant feature of this protocol is that once HyperLogLog values are validated and decoded, aggregation and estimation can be done by anyone, so long as they only try to aggregate within resource classes. In other words, the “over RSA” part of the protocol can be forgotten about entirely once HyperLogLog values are decoded. At that point, it’s as if the collection of logs had been recorded with normal HyperLogLog values keyed on client ID and resource class for each record—which is also a reminder of what must *not* be done with them, the subject of the next section.

**Reporting and publishing counts safely**

Client count estimates are not of any use unless they get published. HyperLogLog allows estimating the number of unique clients in arbitrary slices of requests, but some care should be taken in what slices are actually reported.

It is tempting to think that, because each decoded record is only a small HyperLogLog value, the decoded logs are harmless and could simply be published as raw request-level records. They should not be. In any class small enough that clients are unlikely to share $(b, k)$ pairs, these decoded pairs uniquely identify clients, in much the same way that being a FreeBSD user effectively uniquely identifies someone (we *know* what you're up to, Alex). Publishing decoded HyperLogLog pairs, would also give an attacker an oracle for decoding their HLL values, which, at the very least, would allow them to learn which pairs are rare.

What *is* safe to expose is **aggregate** output for large enough slices. Per-slice cardinality can be published, rather than the per-request $(b, k)$ values themselves. This is safe so long as a **minimum-cardinality floor** is also enforced: when publishing roll-ups of requests with counts and cardinality estimates, refuse to report any slice whose estimated cardinality is below a threshold, folding everything beneath that threshold into a single “other” bucket. The cost in utility is small since HyperLogLog is a poor estimator at small cardinalities anyway—the suppressed slices are the ones whose counts were least reliable. Only reporting slices with cardinality $≥ 1024$ would be reasonable.

These rules also close the inflation channel discussed in [Malicious Clients](11-malicious-clients.md): an attacker who can read fine-grained per-slice counts gains a decoding oracle for their own forged tokens, and denying counts on sub-floor slices removes it. So the cardinality floor pays for itself three times over—anonymity in small populations, safety of published output, and resistance to inflation attacks.
