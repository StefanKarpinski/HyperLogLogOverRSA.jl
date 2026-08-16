# Security Analysis (Proofs)

This section formalizes the security properties of the protocol. The first part proves the central anonymity result: when $N$ is fingerprint-free, two clients with the same HLL value produce indistinguishable token distributions, so the server learns nothing beyond that value. Rather than assuming the specific structure $N = PQ = (4Bp+1)(2^m q+1)$, we work with a general characterization of fingerprint-free moduli — a condition that will also be needed when we analyze how servers can certify their modulus. We assume throughout that $B$ is odd and $m ≥ 3$.

The second part analyzes the integrity side from the client’s perspective: how much can a malicious client inflate count estimates, and why is that effort linear rather than exponential? The third part turns to the server side: how can a server prove that its modulus is fingerprint-free without revealing its factorization, and what does it cost an attacker to forge such a proof?

## Anonymity

We build to the anonymity theorem in three steps: the algebraic definitions it relies on, the theorem and its proof, and an interpretation of what it guarantees.

### Definitions

**Definition.** We write the *“positive Jacobi subgroup”* in $\Z_N^*$ as:

```math
\begin{aligned}
J_N^+
= \Jacobi_N^{-1}(+1)
= \set{\, x \in \Z_N^* \st \Jacobi_N(x) = +1 \,}
\end{aligned}
```

We’ll write the *“negative Jacobi set”* (it’s not a subgroup) as:

```math
\begin{aligned}
J_N^-
= \Jacobi_N^{-1}(-1)
= \set{\, x \in \Z_N^* \st \Jacobi_N(x) = -1 \,}
\end{aligned}
```

Since the Jacobi symbol only takes $±1$ values in $\Z_N^*$, $J_N^- = \Z_N^* \setminus J_N^+$, but it’s convenient to have shorter notation. Note that the Jacobi symbol is only well-defined for odd $N$, so when we talk about $J_N^±$ it implicitly requires that $N$ is odd.

**Definition.** We define the *“white noise”* subgroup of $\Z_N^*$ as:

```math
\begin{aligned}
W_N &= \pm(\Z_N^*)^{B2^{m-1}}
= \set{\, \pm z^{B2^{m-1}} \st z \in \Z_N^* \,}
\subseteq J_N^+
\end{aligned}
```

Since the exponent $B2^{m-1}$ is even (as $m ≥ 3$), every $z^{B2^{m-1}}$ has positive Jacobi symbol, and $-1 \in J_N^+$ because $N ≡ 5 \bmod 8$, so the entire noise group has positive Jacobi symbol. This subgroup is where we sample “noise” values to randomize the parts of $x$ that don’t encode the HyperLogLog value. We previously described deriving individual $w$ values from random $z$ values; here we consider the entire group.

**Definition.**  We call a positive integer, $N$, *“fingerprint-free”* if it is odd and the quotient $J_N^+/W_N$ is cyclic with order dividing $B2^{m-1}$.

This is an intrinsic condition on the noise quotient — it is what the certification of later sections verifies without the factorization — and it is exactly what supplies the homomorphism the anonymity proof needs:

**Proposition.**  If $N$ is fingerprint-free then the embedding $J_N^+/W_N \hookrightarrow C_B \times C_{2^{m-1}}$ extends to a group homomorphism on all of $\Z_N^*$,

```math
\begin{aligned}
\phi : \Z_N^* \to C_B \times C_{2^{m-1}},
\end{aligned}
```

whose kernel meets $J_N^+$ in exactly $W_N$; that is, $\ker(\phi) \cap J_N^+ = W_N$.

**Proof.**  Cyclicity with order dividing $B2^{m-1}$ means $J_N^+/W_N$ embeds into $C_B \times C_{2^{m-1}} \cong C_{B2^{m-1}}$, giving a homomorphism $J_N^+ \to C_B \times C_{2^{m-1}}$ with kernel $W_N$. We extend it across $J_N^-$: for $\xi \in J_N^-$ we have $\xi^2 \in J_N^+$, and its image is a square in $C_B \times C_{2^{m-1}}$ (the odd factor $C_B$ has unique square roots, and the $C_{2^{m-1}}$ component of a square is even, so a root exists), so choose such a square root to be $\phi(\xi)$ and set $\phi(\xi h) = \phi(\xi)\,\phi(h)$ for $h \in J_N^+$. This $\phi$ is a homomorphism on all of $\Z_N^*$ extending the embedding, so restricting it to $J_N^+$ recovers the embedding and $\ker(\phi) \cap J_N^+ = W_N$. $\square$

For a correctly constructed ring we can explicitly write down such a $\phi$. For $x \in C_4 \times C_B \times C_{2^m} \times C_{pq}$ with $\log_g(x) = (a, b, c, d)$ this is:

```math
\begin{aligned}
\phi(x) = (b,\; c \bmod 2^{m-1}) \in C_B \times C_{2^{m-1}}.
\end{aligned}
```

Its kernel meets $J_N^+$ in exactly $W_N$: an element lies in $J_N^+$ when $a + c$ is even, and $\phi(x) = 1$ forces $b = 0$ and $c \in \set{0, 2^{m-1}}$; the latter makes $c$ even, so $a$ is even too, and the surviving elements are precisely $\pm(\Z_N^*)^{B2^{m-1}} = W_N$. (Over all of $\Z_N^*$, $\ker(\phi)$ is twice as large — it also contains the $a$-odd elements of $J_N^-$ — but only its intersection with $J_N^+$ is used below.) Reading $c$ only modulo $2^{m-1}$ is what collapses the two rarest rungs into one saturated level, and $\phi(-1) = 1$ with $-1 \in J_N^+$ confirms $-1 \in W_N$.

When $N$ is fingerprint-free — so that $\phi$ from the Proposition exists — we call $\bar{x} := \phi(x) \in C_B \times C_{2^{m-1}}$ the *essential version* of $x$. If $\bar{f}$ is a function on $C_B \times C_{2^{m-1}}$ we write $f(x) = \bar{f}(\bar{x})$ for the corresponding function on $\Z_N^*$. So $\bar{\triangle}$ is generally the “essential version” of $\triangle$, whether an element or a function.

**Definition.** We generalize our earlier definition of a *“semigenerator”* in a multiplicative group. Let $G = \prod_{i=1}^n C_{\alpha_i}$ be a product of cyclic groups and let $\pi_i: G \to C_{\alpha_i}$ be the canonical projection onto the $i$th component. An element $g \in G$ is called a semigenerator if the projection of $g$ onto each cyclic component is a generator for that component:

```math
\begin{aligned}
\A i: C_{\alpha_i} = \set{\, \pi_i(g)^k \st k \in \Z \,}
\end{aligned}
```

If the ${} \alpha_i$ are pairwise coprime, then $G$ is cyclic and $g$ is a true generator. Otherwise, $G$ is not cyclic and $g$ cannot generate it. We can, however, express any $x \in G$ in terms of its logarithm with respect to the projection of $g$ onto each cyclic component:

```math
\begin{aligned}
\log_g(x) = (\ell_1, \ell_2, \dots, \ell_n)
\in \prod_{i=1}^n \Z_{\alpha_i}
\end{aligned}
```

In what follows, let $\bar{g} \in C_B \times C_{2^{m-1}}$ be a fixed semigenerator.

**Definition.** The *“essential HyperLogLog function”* maps each value in $C_B \times C_{2^{m-1}}$ to its HyperLogLog sample value:

```math
\begin{gathered}
\bar{\hll}: C_B \times C_{2^{m-1}} \to \Z_B \times \set{0, \dots, m-1} \\[0.5em]
\bar{\hll}(\bar{x}) = (b, \tz(c)) \\[0.8em]
\text{where}~\log_{\bar{g}}(\bar{x}) = (b, c)
\end{gathered}
```

Here $c$ ranges over $\Z_{2^{m-1}}$, so $\tz(c)$ ranges over $\set{0, \dots, m-1}$ with $\tz(0) = m-1$—the saturated level into which the two rarest rungs have already collapsed. The first value, $b$, implicitly depends on the choice of semigenerator, $\bar{g}$, whereas the latter, $\tz(c)$, does not: $\tz(c)$ only depends on the multiplicative order of the $C_{2^{m-1}}$ part of $\bar{x}$, which is independent of $\bar{g}$. The higher bits of $c$ do depend on $\bar{g}$, but the position of the last bit does not.

The HyperLogLog function on $\Z_N^*$ is then $\hll(x) = \bar{\hll}(\bar{x})$. In terms of a raw logarithm $\log_g(x) = (a, b, c, d)$ this is $\hll(x) = (b, \min(\tz(c), m-1))$, since $\bar{x}$ reads the geometric coordinate modulo $2^{m-1}$. It depends on the choice of $\bar{g}$ for $\bar\hll$ and, through $\bar{x}$, on the semigenerator $g$; both choices merely permute the output bucket indices.

### The anonymity theorem

**Lemma.** For $\bar{x}, \bar{y} \in C_B \times C_{2^{m-1}}$ we have $\bar{\hll}(\bar{x}) = \bar{\hll}(\bar{y})$ if and only if there exists $t \in \Z$ with $t = 1 \bmod{2B}$ such that $\bar{x}^t = \bar{y}$.

**Proof.**  Denote the logarithms of $\bar{x}$ and $\bar{y}$ as:

```math
\begin{aligned}
\log_{\bar{g}}(\bar{x})
&= (b_1, c_1) \\
\log_{\bar{g}}(\bar{y})
&= (b_2, c_2) \\
\end{aligned}
```

Suppose $\bar{\hll}(\bar{x}) = \bar{\hll}(\bar{y})$. This gives us two equalities:

```math
\begin{aligned}
b_1 &= b_2 \bmod B &
\tz(c_1) &= \tz(c_2) \\
\end{aligned}
```

The second equality implies that there exists odd $u \in \Z_{2^{m-1}}$ such that

```math
\begin{aligned}
uc_1 &= c_2 && \pmod{2^{m-1}} \\
\end{aligned}
```

We can apply the Chinese Remainder Theorem to find $t$ with:

```math
\begin{aligned}
t &= 1 && \pmod{B} \\
t &= u && \pmod{2^{m-1}} \\
\end{aligned}
```

Since $u$ is odd and ${} m ≥ 3$ we know that $t$ is odd as well. This means that $t = 1 \bmod{2B}$ as required. Now check that $\bar{x}^t = \bar{y}$:

```math
\begin{aligned}
\log_{\bar{g}}(\bar{x}^t)
&= (t b_1, t c_1) \\
&= (b_1, u c_1) \\
&= (b_2, c_2)
= \log_{\bar{g}}(\bar{y})
\end{aligned}
```

Thus we have found $t = 1 \bmod{2B}$ such that $\bar{x}^t = \bar{y}$.

In the other direction, suppose we have $t = 1 \bmod{2B}$ such that $\bar{x}^t = \bar{y}$. This gives us two equalities:

```math
\begin{aligned}
b_1 =
t b_1 &= b_2 && \pmod B \\
t c_1 &= c_2 && \pmod{2^{m-1}} \\
\end{aligned}
```

Since $t$ is odd, the second equality implies that $\tz(c_1) = \tz(c_2)$ which means we’ve shown that $\bar{\hll}(\bar{x}) = \bar{\hll}(\bar{y})$. $\square$

**Theorem.** Let $N$ be fingerprint-free. For $x, y \in \Z_N^*$ with the same Jacobi symbol: $\hll(x) = \hll(y)$ if and only if there exists $w \in W_N$ and $t \in \Z$ with $t = 1 \bmod{2B}$ such that $wx^t = y \bmod N$.

**Proof.**  Fingerprint-freeness gives us, via the Proposition, the essential homomorphism $\phi : \Z_N^* \to C_B \times C_{2^{m-1}}$, $x \mapsto \bar{x}$, with $\ker(\phi) \cap J_N^+ = W_N$. The following conditions are equivalent for $x$ and $y$ with $\Jacobi_N(x) = \Jacobi_N(y)$:

1. ``\hll(x) = \hll(y)``
2. ``\bar{\hll}(\bar{x}) = \bar{\hll}(\bar{y})``
3. ``\E\, t = 1 \bmod{2B}`` such that $\bar{x}^t = \bar{y}$
4. ``{} \E\, t = 1 \bmod{2B},\, w \in W_N {}`` such that $wx^t = y$

Equivalence of (1) and (2) is just the definition $\hll(x) = \bar{\hll}(\bar{x})$. The lemma we just proved gives equivalence of (2) and (3). That leaves us to prove equivalence of (3) and (4). We will show that for any odd exponent, $t$, we have:

```math
\begin{aligned}
\bar{x}^t = \bar{y} ~\iff~
\E\, w \in W_N\!:~ wx^t = y
\end{aligned}
```

which means (3) and (4) are logically equivalent. Suppose that $\bar{x}^t = \bar{y}$. Let $w = y x^{-t}$; we show $w \in W_N$ by landing it in both $J_N^+$ and $\ker(\phi)$. Its Jacobi symbol is positive,

```math
\begin{aligned}
\Jacobi_N(w) = \Jacobi_N(y)\,\Jacobi_N(x)^{-t} = 1,
\end{aligned}
```

using $\Jacobi_N(x) = \Jacobi_N(y) \in \set{±1}$ and $t$ odd so that raising to $-t$ doesn’t change the sign — so $w \in J_N^+$. And its essential version is trivial,

```math
\begin{aligned}
\bar{w} = \bar{y}\,\bar{x}^{-t} = 1,
\end{aligned}
```

so $w \in \ker(\phi)$. Together these give $w \in \ker(\phi) \cap J_N^+ = W_N$. In the other direction, suppose $w \in W_N$ with $wx^t = y$; since $W_N \subseteq \ker(\phi)$ we have $\bar{w} = 1$, and

```math
\begin{aligned}
\bar{y} = \overline{wx^t} = \bar{w}\,\bar{x}^t = \bar{x}^t.
\end{aligned}
```

This proves equivalence of (3) and (4). $\square$

### Interpretation

How does this result prove that our scheme preserves client anonymity? Suppose there are two clients with $x_1$ and $x_2$ as their respective client secrets and assume that $\hll(x_1) = \hll(x_2)$. We also assume that the clients follow the protocol so that $\Jacobi_N(x_1) = \Jacobi_N(x_2) = -1$. The server receives $w_1x_1^{t_1}$ from the first client. Can the server tell that it got the message from the first client rather than the second? We can apply the theorem twice to show that it cannot. First, the theorem tells us that when the client chooses $w_1 \in W_N$ and $t_1$ with $t_1 = 1 \bmod{2B}$, we always have:

```math
\begin{aligned}
\hll(w_1 x_1^{t_1}) = \hll(x_1)
\end{aligned}
```

In other words, the value the client sends has the same $\hll$-value as its secret. Next, since $\hll(w_1 x_1^{t_1}) = \hll(x_1) = \hll(x_2)$, the theorem tells us that there exists ${} w_2 \in W_N$ and $t_2 \in \Z$ with $t_2 = 1 \bmod{2B}$ such that

```math
\begin{aligned}
w_2 x_2^{t_2} = w_1 x_1^{t_1}
\end{aligned}
```

This means that either client could have sent the same value—the server cannot distinguish between messages from correctly behaving clients with the same $\hll$ value.

There is also the question of distribution of values: it could be possible that either client _could_ send the same value but it’s much more likely for one of them to do so. This depends on how $w$ and $t$ values are chosen, which we haven’t addressed yet. Fortunately, so long as they are chosen uniformly from publicly known ranges, every possible value is equally likely for each client, so there is truly no way to distinguish between clients with the same HyperLogLog value.

With the anonymity theorem in hand, we turn to the integrity side of the protocol.

## Malicious clients

Plaintext HyperLogLog counting has a serious inflation problem: an attacker can send one request for each bucket, each carrying the maximum geometric value, and instantly push the unique client estimate to its ceiling. How does our protocol compare? Our gold standard is the naive unique client ID approach, where each forged request inflates the count by exactly one—attack effort is linear with coefficient one.

If a malicious client sends $y$ with $\Jacobi_N(y) ≠ -1$, the server will detect it. So we can focus on a malicious client sending $y$ with $\Jacobi_N(y) = -1$. As discussed in the previous section, the values it can build from a fixed twist $x_0$ (with $\Jacobi_N(x_0) = -1$) and the published semigenerator $g$ are those of the form $y = \pm x_0 g^h$; take $y = x_0 g^h$, the other case being identical. Here we use $h$ as just an arbitrary attacker-controlled exponent value, not necessarily the output of a hash function. The question is whether a malicious client can influence $y$ to have large geometric sample values.

The server doesn’t publish any $x_0$ value, so clients have to generate one for themselves. However, since they don’t know the factorization of $N$, they have no idea what its logarithms are. For our analysis, write $\log(x_0) = (a, b, c, d)$, but keep in mind that the attacker has no idea what these values are. The HyperLogLog sample that $y = x_0 g^h$ encodes is $(b + h, \tz(c + h))$. The attacker controls $h$ but doesn’t know $b$ or $c$. Since they don’t know $c$, they cannot force $\tz(c + h)$ to be large or know how large it is for any particular $h$. In order to hit $k ≥ \tz(c + h)$ the attacker needs to happen to choose $h$ whose last $k$ bits complement $c$ perfectly, which occurs with probability $1/2^k$ — exactly the probability of picking a value that good by chance. What they can do, however, is scan through consecutive $h$ values. If they scan $2^k$ consecutive values, they’re guaranteed to hit $h = -c \bmod{2^k}$ for one of those values and they may happen to hit something better.

In the unique ID scheme, each request with a freshly “forged” client ID inflates the client count by one. How does HLL over RSA compare to this when an attacker tries to inflate its estimates? The expected maximum value of $\tz(c + h)$ is $k+1$ when scanning a consecutive block of $2^k$ values. This is slightly better than sending random values. If an attacker sends $B2^k$ spoofed requests they can push each bucket up to an expected value of $k + 1$, which gives a client count estimate of:

```math
\begin{aligned}
\hat{n} = \frac{1}{2 \ln(2)} B 2^{k+1} = \frac{1}{\ln(2)} B2^k
\end{aligned}
```

This uses the original Flajolet *et al.* (2007) estimator for $B ≥ 128$. Better estimators exist for small samples, and should be used for client count estimates, but for this analysis the original asymptotic formula is fine. So, at scale, an attacker can expect to inflate the estimate by $\ln(2)^{-1} \approx 1.44$ per malicious request. This is a rather good result: attack effort is linear and the coefficient is just a bit larger than one—only slightly worse than the unique ID gold standard. In short, whereas the naive unique ID scheme has perfect accuracy, resists inflation well, but has awful privacy, our HLL over RSA scheme is approximate (with tunable accuracy), recovers nearly the same inflation resistance, and provides provable client anonymity.

There is an important assumption hiding inside this linear bound: that the attacker cannot tell which of their forged tokens decoded to large geometric values. The bound holds precisely because, unable to decode, the attacker can’t cherry-pick—every forged request is a fresh blind sample, so pushing the count up takes linear effort. But the whole purpose of the system is to *publish counts*, and a published count is exactly the oracle the attacker otherwise lacks. Inject a token into an otherwise-sparse slice, watch whether the reported cardinality jumps, and you learn whether that token happened to decode high—without breaking any encryption. With that feedback an attacker can curate a collection of high-$k$ tokens, one per bucket, and a curated set is no longer linear: roughly $B$ requests can then push every bucket to level $k$, inflating the estimate to $\sim \frac{1}{\ln 2} B 2^k$, the same exponential leverage as the plaintext scheme. So the linear-resistance result holds specifically against an attacker who *cannot observe per-slice counts*. The mitigation is the minimum-cardinality floor introduced for anonymity ([Resource class sharding](03-counting-users.md#Resource-class-sharding) and the [Protocol Summary](06-protocol-summary.md)): if sub-threshold slices are never reported, the attacker has no sparse slice to probe and no decoding oracle to build. Pleasingly, the same mechanism that protects anonymity in small populations also closes this inflation channel—it is one more reason to enforce the floor, not a separate defense.

It’s worth flagging a stronger mitigation that we have deliberately *not* adopted, since it’s a natural question. One could try to *bind each token to its request*—mixing something request-specific (a server-issued nonce, or the full request path) into the value the client sends—so that a single high-value forgery can’t be reused across many requests. The difficulty is that honest counting needs the opposite: a client’s value must be *stable* across its requests within a class, or we’d be counting requests instead of unique clients. A binding that preserved per-(client, class) stability while denying an attacker reuse of a curated value is not obviously achievable, and in any case the cardinality floor already removes the count oracle that the curated-set attack depends on. So we rely on the floor and leave token-request binding as an open question—to revisit only if some deployment leaks count feedback through a channel the floor can’t cover.

While this is quite a good result, our analysis of resistance to malicious clients rests on a hardness assumption we cannot discharge. The break condition is sharp: an attacker who learns the $c$ coordinate—indeed the whole logarithm—of even one element of $J_N^-$ with respect to $g$ can forge arbitrarily rare samples, choosing $h$ so that $\tz(c + h)$ lands wherever it likes and $h \bmod B$ so as to hit any bucket. So the question is whether an attacker can produce a $J_N^-$ element *together with its logarithm*.

It is essential to be precise about “produce,” because merely *naming* an element of $J_N^-$ is trivial. For our moduli $\Jacobi_N(2) = -1$ always (this is what $N ≡ 5 \bmod 8$ buys), so $2 \in J_N^-$ on the nose; and the fingerprint-freedom certificate openly publishes a list of square roots, roughly half of which land in $J_N^-$. None of these are dangerous, because their logarithms are unknown. What an attacker needs is not a name but a logarithm.

The elements whose logarithms an attacker *can* write down are exactly those built from known powers of $g$ and $-1$. Since $\Jacobi_N(g) = 1$ and $\Jacobi_N(-1) = 1$, that subgroup sits inside $J_N^+$—and in fact fills it:

```math
\begin{aligned}
\gen{g, -1} = J_N^+, \qquad J_N^+ = \gen{g} \sqcup -\gen{g}
\end{aligned}
```

Every product of known powers of $g$ and $-1$ therefore has positive Jacobi symbol. Reaching $J_N^-$ with a known logarithm means exhibiting one more object: a square root of $-1$. A primitive fourth root of unity $i$ has logarithm $(\text{odd}, 0, 2^{m-2}, 0)$ and $\Jacobi_N(i) = -1$, so $i$ is exactly a known-logarithm element of $J_N^-$—and given $i$, the forgery is immediate. Computing such an $i \bmod N$ is not known to be any easier than factoring $N$: we have a one-way chain, factor $N \Rightarrow$ compute $i \Rightarrow$ forge, with no reduction known in either reverse direction. (Indeed a *mixed* square root of $-1$, $\xi = \mathrm{CRT}(-1 \bmod P,\, 1 \bmod Q)$, has $\Jacobi_N(\xi) = 1$ and $\gcd(\xi - 1, N) = Q$, so it simply *is* the factorization; the server plainly cannot publish it.) The gap between what a client can build with known logarithms and what a forgery requires is thus exactly one square root of $-1$.

This is a vaguer argument than I would like to make in a write up full of rigorous proofs. Fortunately, proving that client anonymity is preserved is the much more important result, and for that we do have a solid proof. Inability to provably guarantee resistance to malicious client attacks is more acceptable: that’s a risk that we, as the operators of the system, can choose to take on. If our estimates suddenly start looking ridiculously large, we can begin to wonder if someone has cracked the protocol.

## Malicious servers

So far we’ve worried about malicious clients but have assumed that servers behave as intended. Servers are supposed to generate $N$ with the following arithmetic structure:

```math
\begin{aligned}
N = P Q = (4 B p + 1)(2^m q + 1)
\end{aligned}
```

where $B$ and $m$ are published parameters and $P$, $Q$, $p$ and $q$ are secret primes. If a client knew $P$ and $Q$ they could check the construction of $N$, but the whole point, of course, is that they don’t know the factorization. Otherwise they could easily forge rare HyperLogLog values. What if a server generates $N$ with a different structure? Can that allow them to “fingerprint” and track individual clients? As it turns out, constructing $N$ with different structure than intended *can* allow a server to fingerprint clients. Because of this, our protocol needs a mechanism for servers to convince clients that they have not hidden any fingerprints in the construction of $N$. This only needs to be checked when a client talks to a server initially and downloads new protocol parameters (or when the server issues new protocol parameters).

### The fingerprinting threat

First, let’s get a feel for how a malicious server could fingerprint clients. Here’s an alternative structure that would allow a server to uniquely identify clients:

```math
\begin{aligned}
N = P Q = (2^m B p + 1)(2^m B q + 1)
\end{aligned}
```

In other words, $P-1$ is divisible by $2^m$ instead of the intended $4$ and $Q-1$ is divisible by $B$ when it shouldn’t be. With this alternative structure, $\Z_N^*$ has the following richer multiplicative structure:

```math
\begin{aligned}
\Z_N^* \cong
C_{2^m} \times C_B \times C_p \times
C_{2^m} \times C_B \times C_q
\end{aligned}
```

When the client follows the protocol and picks a random $x_0 \in J_N^-$, it has two random $C_B$ components and two random $C_{2^m}$ components instead of one of each. How much extra identifying information about each client does the malicious server get from this?

First, consider the two $C_B$ components. When the client raises their $x$ value to $t = 1 \bmod{2B}$ both $C_B$ components are preserved. Instead of getting $\log_2(B)$ bits of identifying information from them, the server gets $2 \log_2(B)$ bits. When $B \approx 2^{12}$ that’s 24 bits of fingerprint. Multiplying $x^t$ by a white noise value, $w \in W = \pm(\Z_N^*)^{B2^{m-1}}$, also doesn’t touch either of the two $C_B$ and $C_{2^m}$ components—the definition of $W_N$ is specifically crafted to leave those components intact while randomizing everything else.

Next, consider the two $C_{2^m}$ components. When there is only one $C_{2^m}$ component, taking $x^t$ with random odd $t$ destroys all information in it except for the position of the last bit, which is only $\log_2(m)$ bits of information. When there are two $C_{2^m}$ components, however, the information provided doesn’t just double: since the two $C_{2^m}$ components are scaled in lock-step, their _ratio_ remains fixed, which carries significantly more information. It’s not hard to see that two $C_{2^m}$ components raised to random $t$, convey a full $m$ bits of client fingerprint. Again, multiplying by $w$ doesn’t affect this, by design. So the malicious server gets an additional $m$ bits of fingerprint from the two $C_{2^m}$ components.

With this structure of $N$ then, the malicious server gets a total of $2 \log_2(B) + m$ bits of client fingerprint. For our usual choice of $B$ and $m$, that’s $2\cdot12 + 63 = 87$ bits, which is enough to uniquely identify billions of clients with random fingerprints without significant chance of collisions. If we left this unmitigated, it would be a serious defect in our protocol.

Can we convince a client that we aren’t smuggling fingerprint bits in the structure of $N$ without giving away its factorization? In this particular case, the unusual structure of $N$ is actually quite easy to detect:

- The correct structure has
  - ``P = 5 \bmod 8 \and Q = 1 \bmod 8 \implies N = PQ = 5 \bmod 8``
  - ``P = 1 \bmod B \and Q ≠ 1 \bmod B \implies N = PQ ≠ 1 \bmod B``
- This incorrect structure has
  - ``N = P = Q = 1 \bmod 8``
  - ``N = P = Q = 1 \bmod B``

Unfortunately, not all possible malicious structures are so easy to detect. For example, if $N = PQR$ where $P = Q = 1 \bmod 8$ and $P = Q = 1 \bmod B$ and $R = 5 \bmod 8$ and $R ≠ 1 \bmod B$, then you’d have $N = PQR = 5 \bmod 8$ and $N = PQR ≠ 1 \bmod B$ so $N$ looks normal from simple modular criteria, yet $PQ$ carries the $2 \log_2(B) + m$ bits of client fingerprint from before. To guarantee that a server cannot fingerprint clients, more evidence about the structure of $N$ needs to be provided.

Zero-knowledge proofs (ZKPs) are a popular solution to this kind of problem. I spent a good bit of time going down this rabbit hole. It’s definitely doable. However, every time I sat down to implement zero-knowledge proofs of semiprimality, I found myself getting bogged down in complex and fussy details. This isn’t just a matter of laziness—if the code is that hard to implement, I find it hard to convince myself that it’s fully correct. And if we’re not confident in the correctness of the code that checks whether $N$ has the right structure, then we haven’t really proven anything.

I found myself really wishing for a simpler way to demonstrate the structure of $\Z_N^*$. One that would be easier to have confidence in, both in the sense of having a simpler implementation to review, but also in the sense of having mathematics that’s easier for more people (including myself) to understand. After some research, I found some straightforward evidence that the server can provide to convince clients a server-provided modulus cannot be used to fingerprint clients.

### Criteria for fingerprint-freedom

Recall that an odd integer, $N$, is *fingerprint-free* when $J_N^+/W_N$ is cyclic of order dividing $B2^{m-1}$ — the definition from the anonymity section, which is what guarantees the map $\phi$ that makes two same-sketch clients indistinguishable. These subgroups of $\Z_N^*$ have the following definitions:

```math
\begin{aligned}
J_N^+
&= \set{\, x \in \Z_N^* \st \Jacobi_N(x) = 1 \,}
\\[0.5em]
W_N
&= \pm(\Z_N^*)^{B2^{m-1}}
= \set{\, \pm x^{B2^{m-1}} \st x \in \Z_N^* \,}
\end{aligned}
```

Our proof in that section justified the term “fingerprint-free” by showing that as long as $N$ satisfies this property, a server cannot distinguish between two correctly behaving clients with the same HyperLogLog values. Our task then, is to come up with a way for the server to convince clients that the value of $N$ it sends them is actually fingerprint-free, without revealing its factorization.

We first present a few supporting results about the arithmetic structure of $N$. For all of the following results, we assume that $N$ is positive and odd with prime factorization given by:

```math
\begin{aligned}
N = \prod_{i=1}^D p_i^{n_i}
\end{aligned}
```

where $p_i$ are distinct odd primes and $n_i ≥ 1$. This prime factorization is, of course, uniquely determined by $N$. The ring structure of $\Z_N$ is given by:

```math
\begin{aligned}
\Z_N \cong \prod_{i=1}^D \Z_{p_i^{n_i}}
\end{aligned}
```

The multiplicative structure of each component is:

```math
\begin{aligned}
\Z_{p_i^{n_i}}^* \cong C_{(p_i-1)\,p_i^{n_i-1}}
\end{aligned}
```

Let $g$ be a semigenerator for $\Z_N^*$ so we can write the logarithm of a generic element $x \in \Z_N^*$ as:

```math
\begin{aligned}
\log_g(x) = (\ell_1, \ell_2, \dots, \ell_D)
\in \prod_{i=1}^D \Z_{(p_i-1) \, p_i^{n_i-1}}
\end{aligned}
```

We will often prove results in terms of logarithmic coordinates. These kinds of results are not usually proven with this approach, but since we’ve already worked extensively with this viewpoint, it will hopefully be more illuminating in this context.

**Fact.**  $\gcd_i(p_i-1) \divides \gcd(N-1)$.

**Proof.**  Let $d = \gcd_i(p_i-1)$. Since $d \divides p_i-1$ we have $p_i = 1 \bmod d$, which implies:

```math
\begin{aligned}
N
= \prod_i p_i^{n_i}
= \prod_i 1^{n_i}
= 1 \bmod{d}
\end{aligned}
```

This means that $d \divides N-1$ as required. $\square$

**Definition** (standard)**.** An element $x \in \Z_N$ is called a “quadratic residue” if there exists $r \in \Z_N$ such that $r^2 = x \bmod N$. The set of quadratic residues may be denoted as $(\Z_N)^2$.

**Lemma.** Let $N$ be an odd, positive integer. $N$ has at most two distinct prime factors if and only if for all $\set{x, y} \subseteq J_N^+$ at least one of $\set{x, y, xy}$ is a quadratic residue.

**Proof.**  Since $(p_i-1)\,p_i^{n_i-1}$ is always even, each of the log-coordinates has a well-defined notion of parity. We can understand both the Jacobi symbol and quadratic residues in terms of parity of log-coordinates:

- ``\Jacobi_N(x) = 1`` iff an even number of the $\ell_i$ are odd
- ``x`` is a quadratic residue iff none of the $\ell_i$ are odd

In the $D = 1$ case, $\Jacobi_N(x) = 1$ if and only if $x$ is a quadratic residue, so all three of $x$, $y$ and $xy$ are quadratic residues.
In the $D = 2$ case, let’s write coordinates for $x$, $y$ and $xy$:

```math
\begin{aligned}
\log_g(x) &= (a_x, b_x) \\
\log_g(y) &= (a_y, b_y) \\
\log_g(xy) &= (a_x + a_y, b_x + b_y) \\
\end{aligned}
```

Since $\Jacobi_N(x) = \Jacobi_N(y) = 1$ we know that

```math
\begin{aligned}
a_x &= b_x && \pmod 2 \\
a_y &= b_y && \pmod 2 \\
a_x + a_y &= b_x + b_y && \pmod 2 \\
\end{aligned}
```

One of the following has to hold:

```math
\begin{aligned}
a_x = b_x &= 0 && \pmod 2 \\
a_y = b_y &= 0 && \pmod 2 \\
a_x + a_y = b_x + b_y &= 0 && \pmod 2 \\
\end{aligned}
```

Each condition means one of $x$, $y$ or $xy$ is a quadratic residue.
In the $D = 3$ case, we have a counterexample:

```math
\begin{aligned}
\log_g(x) &= (1, 1, 0) \\
\log_g(y) &= (0, 1, 1) \\
\log_g(xy) &= (1, 2, 1) \\
\end{aligned}
```

Each value has positive Jacobi symbol, yet none of them is a quadratic residue. The same essential counterexample works for $D > 3$ as well, letting trailing coordinates be zeros. $\square$

Combining these results, we get the following claim.

**Claim.** If the following conditions are satisfied:

- ``N = 5 \bmod 8``
- ``\gcd(B, N) = \gcd(B, N-1) = 1``
- For all $x, y \in J_N^+$ one of $\set{x, y, xy}$ is a quadratic residue

then $N$ is fingerprint-free.

**Proof.**  From the prior lemma the last condition is equivalent to $N$ having at most two distinct prime factors. We compute the noise quotient $\Z_N^*/W_N$ and then read off its index-2 Jacobi-kernel $J_N^+/W_N$, showing that kernel is cyclic of order dividing $B2^{m-1}$ — the definition of fingerprint-freedom.

Two general tools. First, for a cyclic group $C_M$ and any exponent $E$, the power map $x \mapsto x^E$ has image of index $\gcd(M, E)$, so $C_M/(C_M)^E \cong C_{\gcd(M, E)}$. Applied factor by factor to $\Z_N^* \cong \prod_i C_{M_i}$ with $E = B2^{m-1}$,

```math
\begin{aligned}
\Z_N^* / (\Z_N^*)^{B2^{m-1}} \cong \prod_i C_{\gcd(M_i,\, B2^{m-1})}.
\end{aligned}
```

Second, $W_N = \pm(\Z_N^*)^{B2^{m-1}}$ adjoins the single order-2 element $-1$, so the noise quotient is one further step:

```math
\begin{aligned}
\Z_N^*/W_N \cong \Big( \prod_i C_{\gcd(M_i,\, B2^{m-1})} \Big) \Big/ \gen{\,\overline{-1}\,},
\end{aligned}
```

where $\overline{-1}$ is the image of $-1$. We take the odd and 2-power parts of this in turn.

**One prime factor.**  $N = P^j$, so $\Z_N^*$ is cyclic of order $(P-1)P^{j-1}$. A prime power $P^j \equiv 5 \bmod 8$ forces $P \equiv 5 \bmod 8$ (the residues $1, 3, 7 \bmod 8$ generate cyclic subgroups of $(\Z/8)^*$ that never contain $5$), so $\tz(P-1) = 2$ and $\Z_N^* \cong C_4 \times C_U$ with $U$ odd. Then $\gcd(4, B2^{m-1}) = 4$ (since $m ≥ 3$) and $\gcd(U, B2^{m-1}) = \gcd(U, B)$, so

```math
\begin{aligned}
\Z_N^* / (\Z_N^*)^{B2^{m-1}} \cong C_4 \times C_{\gcd(U, B)}.
\end{aligned}
```

The gcd condition pins the odd factor: $\gcd(B, U) \divides \gcd(B, P-1)$, which divides $\gcd(B, N-1) = 1$ (as $N - 1 = (P-1)\,(P^{j-1} + \dots + 1)$ is a multiple of $P-1$), so $\gcd(U, B) = 1$. That leaves $C_4$; and $\overline{-1}$ is its order-2 element, so $\Z_N^*/W_N \cong C_2$. The Jacobi symbol is injective on this $C_2$ (its generator is the image of a coordinate-$a = 1$ element, with $\Jacobi_N = -1$), so the Jacobi-kernel $J_N^+/W_N$ is trivial — cyclic of order dividing $B2^{m-1}$. Fingerprint-free.

**Two prime factors.**  $N = P^j Q^k$. Writing $s = \tz(P-1)$ and $n = \tz(Q-1)$,

```math
\begin{aligned}
\Z_N^* \cong C_{2^s} \times C_U \times C_{2^n} \times C_V,
\quad
U = \tfrac{P-1}{2^s}\,P^{j-1},
\quad
V = \tfrac{Q-1}{2^n}\,Q^{k-1},
\end{aligned}
```

with $U, V$ odd. Put $B_U = \gcd(B, U)$ and $B_V = \gcd(B, V)$; these are the bucket-carrying factors. Since $B_U \divides P-1$ and $B_V \divides Q-1$, we have $\gcd(B_U, B_V) \divides \gcd(B, P-1, Q-1)$, which divides $\gcd(B, N-1) = 1$ (any common divisor of $P-1$ and $Q-1$ divides $N - 1 = (P-1)Q + (Q-1)$); so $B_U, B_V$ are coprime, and being coprime divisors of $B$, $B_U B_V \divides B$. The **odd part** of $\Z_N^*/W_N$ is therefore $C_{B_U} \times C_{B_V} \cong C_{B_U B_V}$, which embeds into $C_B$; $\overline{-1}$, being 2-power, does not touch it.

For the **2-part**, $N \equiv 5 \bmod 8$ pins the 2-Sylow. On coordinates $(a, b, c, d)$ the Jacobi symbol is $\Jacobi_N = (-1)^{a+c}$, and $-1 = (2^{s-1}, 0, 2^{n-1}, 0)$, so $\Jacobi_N(-1) = (-1)^{[s=1] + [n=1]}$. Now $N \equiv 5 \bmod 8$ gives $N \equiv 1 \bmod 4$, forcing $\Jacobi_N(-1) = +1$ and hence $[s=1] = [n=1]$; and $N \not\equiv 1 \bmod 8$ rules out $s, n \ge 3$ together (that would give $N \equiv 1 \bmod 8$). So either $s = n = 1$, or — relabelling the primes so $P$ carries the smaller 2-part — $s = 2$ and $n \ge 2$. In every case $s \le 2 \le m-1$, so $\gcd(2^s, B2^{m-1}) = 2^s$, while $\gcd(2^n, B2^{m-1}) = 2^r$ with $r = \min(n, m-1)$. The 2-part of $\Z_N^*/(\Z_N^*)^{B2^{m-1}}$ is thus $C_{2^s} \times C_{2^r}$, and $\overline{-1}$ has 2-part $(2^{s-1},\, 2^{n-1} \bmod 2^r)$, whose first component $2^{s-1}$ is the order-2 element of $C_{2^s}$. Quotienting by it:

- **$s = n = 1$:** the 2-part is $C_2 \times C_2$ and $\overline{-1} = (1, 1)$ is the diagonal, so the quotient is $C_2$.
- **$s = 2$, $n \ge m$:** then $r = m-1$ and $2^{n-1} \equiv 0 \bmod 2^{m-1}$, so $\overline{-1} = (2, 0)$ and the quotient is $C_2 \times C_{2^{m-1}}$ — the honest case, Jacobi bit times a geometric $C_{2^{m-1}}$.
- **$s = 2$, $2 \le n \le m-1$:** then $r = n$ and $2^{n-1} \bmod 2^n = 2^{n-1}$, so $\overline{-1} = (2, 2^{n-1})$ is a diagonal order-2 element and the quotient is cyclic of order $2^{n+1}$, with $n + 1 \le m$.

Now take the Jacobi-kernel. The Jacobi symbol is a character of $\Z_N^*/W_N$ whose kernel is exactly $J_N^+/W_N$, an index-2 subgroup; its 2-part is the index-2 subgroup of the 2-parts above, namely — respectively — the trivial group ($s = n = 1$), the cyclic $C_{2^{m-1}}$ (the diagonal kernel $\gen{(1,1)}$ inside $C_2 \times C_{2^{m-1}}$, since here $\Jacobi_N$ is the *sum* of the two coordinates' parities), and the cyclic $C_{2^n}$ (the index-2 subgroup of $C_{2^{n+1}}$, with $n \le m-1$). Each is cyclic of order dividing $2^{m-1}$. Together with the odd part $C_{B_U} \times C_{B_V} \cong C_{B_U B_V}$ (order dividing $B$), $J_N^+/W_N$ is cyclic of order dividing $B2^{m-1}$, so $N$ is fingerprint-free. $\square$

This gives us a concrete set of criteria on $N$, which, taken together, imply that $N$ is fingerprint-free. Of course, the obvious question is how can a client be convinced that for *every* pair $\set{x, y} \subseteq J_N^+$ one of $\set{x, y, xy}$ has a square root? This obviously cannot be checked exhaustively by client or server, since $J_N^+$ is huge for realistic $N$. The next section gives results that allow us to design a protocol that lets a server convince clients that ${} N$ is overwhelmingly likely to be fingerprint-free.

### Upper bound on quadratic residues

The server cannot provide a square root for one of $\set{x, y, xy}$ for every pair $\set{x, y} \subseteq J_N^+$. But if the client can challenge the server to produce square roots for arbitrary pairs, then the server can interactively convince the client that either every pair has this property or the server has gotten very lucky. And if the client can challenge the server to provide a square root for as many pairs as it wants, then it can require the server to be implausibly lucky. Moreover, an interactive protocol like this can be converted into non-interactive certificates in a standard way by using cryptographic hashing to produce random elements instead of having the client choose them. In order to know how many challenges a server should be required to answer, however, we need an upper bound on the fraction of pairs that have square roots when $N$ has incorrect structure. The following theorem provides this upper bound.

**Theorem.** Let $N$ be an odd, positive integer. If $N$ has ${} D$ distinct prime factors, then the fraction of pairs ${} \set{x, y} \subseteq J_N^+ {}$ with at least one of $\set{x, y, xy}$ being a quadratic residue is

```math
\begin{aligned}
\frac{3}{2^{D-1}} - \frac{2}{2^{2(D-1)}}
\end{aligned}
```

In particular, for $D ≥ 3$ it is at most $5/8$.

**Proof.**  It’s straightforward to see that the probability of $x \in J_N^+$ being a quadratic residue is $1/2^{D-1}$:

- There are $D$ even log-coordinates in $\log_g(x) = (\ell_1, \ell_2, \dots, \ell_D)$ which gives $D$ independent parity bits in $\Z_N^*$
- Requiring $x \in J_N^+$ forces $\sum_i \ell_i = 0 \bmod 2$, removing a single independent parity bit, leaving $D-1$ parity bits
- ``x`` is a quadratic residue if and only if all the parity bits are even, which has a probability in $J_N^+$ of $1/2^{D-1}$

Now, define three events:

- ``X``: $x$ is a quadratic residue
- ``Y``: $y$ is a quadratic residue
- ``Z``: $xy$ is a quadratic residue

Note that any two of these events are independent, but any two of them imply the third. Thus, we have:

```math
\begin{aligned}
\Pr(X \cap Y \cap Z)
&= \Pr(X \cap Y) \\
&= \Pr(X \cap Z) \\
&= \Pr(Y \cap Z)
= 1/2^{2(D-1)}
\end{aligned}
```

Apply the [inclusion-exclusion principle](https://en.wikipedia.org/wiki/Inclusion%E2%80%93exclusion_principle) to get:

```math
\begin{aligned}
\Pr(X \cup Y \cup Z)
&= \Pr(X) + \Pr(Y) + \Pr(Z)\\
&-\Pr(X \cap Y) - \Pr(X \cap Z) - \Pr(Y \cap Z) \\
&+\Pr(X \cap Y \cap Z) \\
&= \frac{3}{2^{D-1}} - \frac{2}{2^{2(D-1)}}
\end{aligned}
```

Reassuringly, for $D ≤ 2$ this is one—anything else would contradict our prior theorem. For $D = 3$ it is $5/8$ and it decreases as $D$ increases, so $5/8$ is the maximum for odd $N$ with more than two prime factors. $\square$

This bound allows a protocol whereby a server can convince a client that $N$ is fingerprint-free by providing enough quadratic residues.

### Certifying a good modulus

Based on these results, we can design a protocol for a server to convince clients that $N$ is fingerprint-free. First, the client checks that $N = 5 \bmod 8$, that $\gcd(B, N) = 1$, and that $\gcd(B, N-1) = 1$. These are simple numerical checks.

Why exactly $N \equiv 5 \bmod 8$? Two independent requirements fix it. First, *integrity*: $-1$ must not be a free $J_N^-$ forgery, i.e. $\Jacobi_N(-1) = (-1)^{(N-1)/2} = +1$, which holds iff $N \equiv 1 \bmod 4$ — this is what rules out the $N \equiv 3 \bmod 4$ of a would-be modulus with $-1 \in J_N^-$. Second, *anonymity*: at most one prime may carry a large 2-Sylow factor, since two of them would give the token a second geometric-like coordinate — a fingerprint. That is governed by $\tz(N-1)$: $N \equiv 5 \bmod 8$ forces $\tz(N-1) = 2$, pinning the smaller prime’s 2-part to $C_4$ (or both primes' to $C_2$), whereas $N \equiv 1 \bmod 8$ (that is $\tz(N-1) \ge 3$) would permit two large 2-Sylow factors that the client cannot rule out without the factorization. Requiring $N \equiv 1 \bmod 4$ while forbidding $N \equiv 1 \bmod 8$ is exactly $N \equiv 5 \bmod 8$. As a bonus, $5 \equiv -3 \bmod 8$ gives $\Jacobi_N(2) = -1$, so $2$ is a ready-made, publicly nameable element of $J_N^-$.

The client is then ready to be convinced that $N$ has at most two prime divisors. The interactive version is:

> The client picks $n$ random pairs $\set{x, y} \subseteq J_N^+$ and challenges the server to produce $r \in \Z_N^*$ for each pair such that $r^2 \in \set{x, y, xy} \bmod N$.

This convinces the client that there’s at most a $(5/8)^{n}$ probability that $N$ has more than two distinct prime factors. If the client wants to be convinced to a probability of $\nfrac{1}{\alpha}$ for some large $\alpha$, they should choose $n$ such that:

```math
\begin{alignedat}{2}
\left(\frac{5}{8}\right)^{n}\! &≤ \frac{1}{\alpha}
~~&\iff~
n &≥ \frac{\log_2(\alpha)}{\log_2(8)-\log_2(5)}
\end{alignedat}
```

We can turn this into a non-interactive protocol by picking a cryptographic hashing scheme to generate pairs of values in such a way that the server cannot influence them. The hash should include $N$ as an input so that the set of hash-generated values is not fixed and must be computed after choosing a candidate modulus. This prevents generating $n$ pairs and then solving for an $N$ that happens to have quadratic residues for those pairs. In order to generate an incorrectly structured $N$ that passes this test, a malicious server would have to try about $\alpha$ candidate $N$ values before expecting to find one that passes. Since $n$ grows only as $\log_2(\alpha)$, we can demand an enormous amount of attacker work while keeping the certificate small—but we have to choose $\alpha$ with the real cost of a forgery attempt in mind, and here it pays to be pessimistic, because a successful forgery is catastrophic. The fingerprinting structures from the start of this section yield on the order of $2\log_2(B) + m \approx 87$ bits of fingerprint—enough to deanonymize every client of the system at once—so $\alpha$ is guarding against total failure, not a marginal leak.

It’s tempting to reason about the per-attempt cost in terms of the $n$ quadratic-residue checks, but that reasoning cuts the wrong way: the attacker is forging *their own* candidate moduli, so they know each candidate’s factorization, which makes the QR checks essentially free for them. Worse for the defender, the attacker fails fast—each challenge pair passes a bad modulus with probability at most $5/8$, so on a doomed candidate they bail after two or three pairs and pay the full $n$ square roots only on the eventual success. The real per-attempt cost is dominated by generating a constrained candidate modulus, and the search is embarrassingly parallel. We should not lean on QR-testing being expensive; for this attacker, it isn’t.

To put numbers on it, take a deliberately attacker-favorable estimate of $10^{11}$ modular exponentiations per second of aggregate compute, and pretend each attempt costs only a single modexp (generating a constrained prime is much more expensive, so this *understates* the attacker’s work, which is the safe direction for a security margin):

| $\alpha$ | candidates | wall-clock | $n$ | cert size |
|:---:|:---:|:---:|:---:|:---:|
| $2^{50}$ | $1.1\times10^{15}$ | ~3 hours | 74 | ~9.5 KB |
| $2^{64}$ | $1.8\times10^{19}$ | ~6 years | 95 | ~12 KB |
| $2^{80}$ | $1.2\times10^{24}$ | ~$4\times10^5$ years | 118 | ~15 KB |
| $2^{112}$ | infeasible | infeasible | 166 | ~21 KB |
| $2^{128}$ | infeasible | infeasible | 189 | ~24 KB |

Because $n$ grows only logarithmically in $\alpha$, buying a gigantic safety margin is nearly free: the certificate stays a few dozen kilobytes and the server’s square-root computation stays sub-second. The recommendation is therefore to set $\alpha \ge 2^\lambda$ for the system’s overall security level $\lambda$, and never below $2^{80}$. The simplest safe choice is $\alpha = 2^{112}$ ($n = 166$), which makes finding a malicious $N$ that passes the test so unlikely that it is no longer a concern.

The non-interactive version of this protocol serves as a certificate of fingerprint-freedom for a published $N$ value. The certificate structure contains:

- ``B`` — the number of buckets
- ``m`` — the geometric value cap
- ``N`` — the ring modulus
- ``g`` — a server-selected semigenerator for $\Z_N^*$
- ``\text{sqrts}`` — a list of square roots

When downloading a new ring structure, a client checks the following requirements based on the data in this certificate:

- ``N = 5 \bmod 8``
- ``\gcd(B, N) = 1``
- ``\gcd(B, N-1) = 1``
- That enough square roots are provided
- That all the square roots are valid

That’s it. Once the client has done this, it can safely use $N$ and proceed with the protocol for generating and sending $y = xw^t$ values with each request, confident that it has provable anonymity.
