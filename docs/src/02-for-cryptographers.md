# For Cryptographers

This is a compressed, self-contained tour of the construction for readers comfortable with finite abelian groups and RSA. Everything here is developed in full in the numbered sections; the goal of this page is to let you see the whole shape quickly and decide where to dig in.
Readers preferring a more gradual development from first principles should start with the [next section](03-counting-users.md).

The protocol decouples into two parts that can be understood independently:

- a **HyperLogLog (HLL) sketch**, which carries the *entire* privacy model, and
- an **RSA encoding** of that sketch, which carries the *entire* integrity model and is provably privacy-neutral.

The design target is a per-request token $y \in \Z_N^*$ such that

1. **it reveals exactly the client’s HLL sketch value** $(b, k)$ and nothing else;
2. **it can be freshly re-randomized on every request** without knowing the factorization of $N$, so repeated tokens from one client are unlinkable beyond sharing a sketch value; and
3. **its sketch value cannot be read or biased by the client** without that factorization, so clients cannot manufacture rare values to inflate counts.

Properties (1) and (2) are what make it private; (3) is what makes the counts trustworthy.

## What leaks: the HLL sketch

Fix public parameters: an odd bucket count $B$ and a rank parameter $m \ge 3$. Each client holds a persistent per-class sketch

```math
\hll = (b, k) \in \Z_B \times \set{0, 1, \dots, m-1},
```

with $b$ uniform and $k$ geometric, $\Pr[k = j] = 2^{-(j+1)}$ for $j < m-1$ and $\Pr[k = m-1] = 2^{-(m-1)}$ (the top level absorbs the tail). The server takes the per-bucket maximum of $k$ over any subset of requests and applies a standard HLL estimator to recover the number of distinct clients in that subset. Everything the server is *supposed* to learn about a client is contained in $(b, k)$; the privacy argument — on-average anonymity within an equal-sketch class, and statistical independence across classes via sharding — is made entirely at this level in [Anonymously Counting Users](03-counting-users.md). The remainder of this page is about encoding $(b, k)$ into an RSA group so that a client can neither read nor forge it, *without changing what it leaks*.

## The ring and its coordinates

Take

```math
N = PQ = (4Bp + 1)(2^m q + 1),
```

where $P, Q, p, q$ are distinct odd primes coprime to $B$, which is odd. Then $P \equiv 5$ and $Q \equiv 1 \pmod 8$, so $N \equiv 5 \pmod 8$, and because $P - 1 = 4Bp$ and $Q - 1 = 2^m q$ split into pairwise-coprime cyclic factors ($B$ odd is what keeps $C_B$ off the 2-part),

```math
\Z_N^* \;\cong\; C_4 \times C_B \times C_{2^m} \times C_{pq}.
```

Fix a **semigenerator** $g$: an element whose projection into each cyclic factor generates that factor. Every $x \in \Z_N^*$ then has well-defined coordinates

```math
\log_g x = (a, b, c, d) \in \Z_4 \times \Z_B \times \Z_{2^m} \times \Z_{pq},
```

multiplication adds coordinates, and exponentiation by $t$ scales them by $t$. The sketch is read off the middle two:

```math
\hll(x) = \big(b,\; \min(\tz(c),\, m-1)\big),
```

where $\tz(c) = v_2(c)$ is the 2-adic valuation (trailing-zero count, with $\tz(0) = m$). For uniform $x$, $b$ is uniform on $\Z_B$ and $\tz(c)$ is geometric — exactly the HLL distribution, capped at $m-1$ because the noise below washes out the top bit of $c$ (made rigorous in [Geometric RSA rings](04-hyperloglog-over-rsa.md#Geometric-RSA-rings)). The remaining coordinates are noise to be washed out: $a \in \Z_4$ contributes only its parity, which tracks the Jacobi symbol — its high bit is randomized away — and $d \in \Z_{pq}$ is pure per-client identity entropy.

## The token, and why it leaks only the sketch

A client draws a persistent secret $x_0 \in J_N^-$ (Jacobi symbol $-1$) and derives its per-class secret by hashing into the exponent, $x = x_0\, g^{H(x_0, \text{class})}$. Per request it samples

```math
w \in W = \pm(\Z_N^*)^{B2^{m-1}}, \qquad
t \in T = \set{\, 2Bi + 1 \st i \in [0,\, 2^{m-1}) \,},
```

— both samplable without the factorization — and sends $y = w\,x^{t} \bmod N$. The server discards any $y$ with $\Jacobi_N(y) \ne -1$.

The whole construction is in what $(w, t)$ do coordinate-by-coordinate:

| coord. | factor    | multiplier $w \in W$   | exponent $(\cdot)^t,\ t = 2Bi+1$ | net effect / exposed                                                         |
| ------ | --------- | ---------------------- | -------------------------------- | ---------------------------------------------------------------------------- |
| $a$    | $C_4$     | high bit $\to$ uniform | parity fixed ($t$ odd)           | only the parity survives — pinned by $J_N^-$; the high bit is noise          |
| $b$    | $C_B$     | fixed                  | fixed ($t \equiv 1 \bmod B$)     | **exposed: the bucket**                                                      |
| $c$    | $C_{2^m}$ | top bit $\to$ uniform  | $\times$ uniform odd unit        | randomized within its valuation class — **only $\tz(c)$, capped at $m-1$, survives: the rank** |
| $d$    | $C_{pq}$  | $\to$ uniform          | (irrelevant)                     | washed out                                                                   |

Two facts make the table work:

- ``W`` is $\pm$ the image of the $B2^{m-1}$-power map. The $B2^{m-1}$-power map annihilates $C_B$ and — since $m \ge 3$ — the whole $C_4$ factor, sends $C_{2^m}$ onto just its top bit $\set{0, 2^{m-1}}$, and, as $\gcd(pq,\, B2^{m-1}) = 1$, acts as an automorphism on $C_{pq}$; the $\pm$ sign then adjoins $-1 = (2, 0, 2^{m-1}, 0)$. So $W = \set{0, 2} \times \set{0} \times \set{0, 2^{m-1}} \times C_{pq}$: multiplying by uniform $w$ replaces $d$ with fresh randomness, randomizes the high bit of $a$ and the top bit of $c$, and touches nothing else.
- For $t = 2Bi + 1$: $t \equiv 1 \pmod B$ and $t$ is odd, so $b$ and the parity of $a$ are fixed; and as $i$ ranges over $[0, 2^{m-1})$, $t \bmod 2^m$ runs over **every odd residue exactly once**, so $t$ acts as a uniform unit on $c$ — preserving $v_2(c) = \tz(c)$ but uniformizing $c$ within that valuation class.

Hence the $(w, t)$-orbit of $x$ is exactly the fiber $\hll^{-1}(\hll(x)) \cap J_N^-$, and $y$ is **uniform** on that fiber, not merely supported on it. Two honest clients with the same sketch induce *identical* token distributions, and the server’s entire view is a function of $(b, \min(\tz(c), m-1))$. This is the anonymity theorem of [Proof of Anonymity](05-security-analysis.md); the only bookkeeping subtlety is the 2-torsion ($C_4$/Jacobi) accounting, which is exactly why tokens are confined to $J_N^-$. The proof is phrased in these semigenerator coordinates throughout — if there is a more standard way to present it, I’d welcome the pointer.

## Integrity I — malicious clients (unforgeable rank)

A client who could steer $\hll(x)$ toward large $k$ could inflate counts cheaply, so we need clients unable to do better than chance. Without the factorization a client knows neither $b$ nor $c$ for any value it can produce. It can walk an arithmetic progression of exponents, shifting $c$ by a chosen $h$, but it cannot *see* $\tz(c + h)$, so reaching rank $\ge k$ costs $\sim 2^k$ tries — exactly chance. Inflation is therefore linear at $\approx 1/\ln 2 \approx 1.44$ forged requests per unit of count (the unique-ID gold standard is $1$); see [Malicious clients](05-security-analysis.md#Malicious-clients), which also notes that this linear bound assumes the attacker has no per-slice count oracle, and that a published-count interface must therefore floor small slices.

This is the **one acknowledged gap**: there is no reduction from biasing the rank to factoring or another standard RSA assumption. The evidence is structural. Everything a client can write down *together with its logarithms* lies in $\gen{g, -1}$, which is exactly $J_N^+$; forging a rare rank instead needs a $J_N^-$ element whose $c$-coordinate is known, and the smallest such object is a square root of $-1$ — a primitive fourth root of unity $i$, with $\Jacobi_N(i) = -1$ and structural logarithms $(\text{odd}, 0, 2^{m-2}, 0)$. Producing $i \bmod N$ is not known to be any easier than factoring: the honest statement is a one-way chain — factor $N \Rightarrow$ compute $i \Rightarrow$ forge — with no reduction known in either reverse direction. But this is a conjecture. It is also the *less critical* half of the formalization: a broken integrity bound costs accuracy, while a broken privacy bound costs anonymity; only the latter is proven. A reduction proving that biasing the rank is hard would be the most valuable single addition to this work.

## Integrity II — malicious servers (fingerprint-freedom)

The server chooses $N$, and a maliciously structured $N$ can turn the token into a fingerprint. For instance, if, instead of the above structure, we have

```math
N = (2^m B p + 1)(2^m B q + 1)
```

this gives $\Z_N^*$ two $C_B$ and two $C_{2^m}$ factors; the same honest token then exposes $2\log_2 B + m \approx 87$ bits — total deanonymization. So a client must verify, without learning the factorization, that $N$ has benign structure. The exact criterion ([Malicious servers](05-security-analysis.md#Malicious-servers)) is

```math
\begin{gathered}
N \equiv 5 \!\!\pmod 8 \\[0.5em]
\gcd(B, N) = \gcd(B, N-1) = 1 \\[0.5em]
J_N^+ / W_N \text{ cyclic of order dividing } B2^m.
\end{gathered}
```

The third condition is what caps the leak at one sketch’s worth of bits, and — combined with the first two — it holds **iff $N$ has at most two distinct prime factors**.

That equivalence rests on a clean quadratic-residue criterion. For odd $N$ with $D$ distinct prime factors, every pair $\set{x, y} \subseteq J_N^+$ has at least one of $\set{x, y, xy}$ a quadratic residue **iff $D \le 2$**; and for $D \ge 3$, a uniform pair has this property with probability at most $5/8$. So the server proves $D \le 2$ by answering challenges: for random pairs $\set{x, y} \subseteq J_N^+$ it returns a square root of one of $x$, $y$, $xy$. A $D \ge 3$ modulus survives $n$ independent challenges with probability $\le (5/8)^n$; to force a forging server to try $\ge \alpha$ candidate moduli, take

```math
n = \ceil{\frac{\log_2 \alpha}{\log_2(8/5)}}.
```

Drawing the challenge pairs by hashing $N$ makes this a non-interactive certificate $(B, m, N, g, \text{square roots})$. Verification is one Jacobi symbol, two gcds, and $n$ squarings.

Two cautions on $\alpha$ (detailed in [Malicious servers](05-security-analysis.md#Malicious-servers)): a forger knows its *own* candidate’s factorization, so the QR checks are free for it and it fails fast, meaning $\alpha$ must be sized against attacker wall-clock rather than a bare trial count. Since $n$ grows only as $\log \alpha$, the certificate stays tiny regardless — so use $\alpha \ge 2^\lambda$ for the system security level $\lambda$, never below $2^{80}$, with $\alpha = 2^{112}$ ($n = 166$, ~21 KB) the comfortable default.

## What is proven, and what is not

| Property | Status | Where |
|---|---|---|
| $\hll$ induces the correct HLL distribution | elementary; empirically confirmed | [Geometric RSA rings](04-hyperloglog-over-rsa.md#Geometric-RSA-rings) |
| Anonymity: equal-sketch tokens are identically distributed (uniform on the fiber) | **proven** | [Proof of Anonymity](05-security-analysis.md) |
| Fingerprint-freedom + a certificate with explicit $(5/8)^n$ soundness | **proven** | [Malicious servers](05-security-analysis.md#Malicious-servers) |
| Rank-unforgeability by malicious clients | argued, **not reduced** to a standard assumption | [Malicious clients](05-security-analysis.md#Malicious-clients) |

Feedback is welcome everywhere but especially on these three points:
1. A more standard presentation of the semigenerator-coordinate argument in [Proof of Anonymity](05-security-analysis.md);
2. Any route to a hardness reduction for client-side rank-unforgeability;
3. Whether the geometric-distribution argument in [Geometric RSA rings](04-hyperloglog-over-rsa.md#Geometric-RSA-rings) wants more rigor than it currently has.
