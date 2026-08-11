# Security Analysis (Proofs)

This section formalizes the security properties of the protocol. The first part proves the central anonymity result: when $N$ is fingerprint-free, two clients with the same HLL value produce indistinguishable token distributions, so the server learns nothing beyond that value. Rather than assuming the specific structure $N = PQ = (2Bp+1)(2^m q+1)$, we work with a general characterization of fingerprint-free moduli — a condition that will also be needed when we analyze how servers can certify their modulus. We assume throughout that $B = 2 \bmod 4$ and $m ≥ 3$.

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
W_N &= (\Z_N^*)^{(B/2)2^m}
= \set{\, z^{(B/2)2^m} \st z \in \Z_N^* \,} ≤ J_N^+
\end{aligned}
```

Since $(B/2)2^m$ is even, every element has positive Jacobi symbol, so this is a subgroup of $J_N^+$. This subgroup is where we sample “noise” values to randomize the parts of $x$ that don’t encode the HyperLogLog value. We previously described deriving individual $w$ values from random $z$ values; here we consider the entire group.

The exponent takes the *odd* part of $B$, and it matters that it does. Raising to a power collapses each cyclic factor only as far as that power’s own $2$-adic valuation reaches. An exponent of $B2^m$ has valuation $m+1$ once $B$ is even, which would leave the $2$-part of the quotient truncated one level too high, and a server could then choose $2^{m+1} \divides Q-1$ to widen what it can distinguish while passing every check a client makes. Since $\nfrac{B}{2}$ is odd, $(B/2)2^m$ has valuation exactly $m$ and no choice of $Q$ reaches past the intended ladder.

**Definition.**  We call a positive integer, $N$, *“fingerprint-free”* if it is odd and the quotient by the white noise subgroup embeds into the HyperLogLog value group:

```math
\begin{aligned}
\Z_N^*/W_N \hookrightarrow C_{2B} \times C_{2^m}
\end{aligned}
```

We write $\phi: \Z_N^* \to C_{2B} \times C_{2^m}$ for the quotient map followed by such an embedding, so that $\ker(\phi) = W_N$ exactly.

The definition bounds the quotient from *above*, and only from above: it asks that dividing out the noise leave no more structure than the protocol means to reveal. A modulus whose quotient is *smaller* than $C_{2B} \times C_{2^m}$ — a degenerate ring that collapses part of the value space — satisfies it. Such a ring counts badly, but it cannot be used to tell two clients apart, which is the property the term names.

It may seem strange that a definition about anonymity says nothing about $J_N^±$, since a client’s secret is required to lie in $J_N^-$ and the noise lies in $J_N^+$. One could instead ask for a homomorphism onto $C_B \times C_{2^m}$ whose kernel meets the two Jacobi classes in $+W_N$ and $-W_N$ respectively. That formulation has to name an element of $J_N^-$ to relate the two halves — and the only one available without the factorization is $-1$, which is in $J_N^-$ exactly when $N = 3 \bmod 4$. Since a forger who can name an element of $J_N^-$ can also forge rare HyperLogLog values (see [Malicious clients](05-security-analysis.md#Malicious-clients)), that is the last thing we want the modulus to supply. The $±$ construction turns out to be an artifact of a target group one factor of two too small: keeping the whole $C_{2B}$ coordinate rather than dropping its $2$-part removes the need for it.

The quotient is also insensitive to which of the two admissible shapes $N$ has. With $B$ odd, $C_2 \times C_B \cong C_{2B}$; with $B = 2 \bmod 4$, $C_4 \times C_{B/2} \cong C_{2B}$. Only $B \bmod 4$ separates the two, and that is precisely the dial deciding whether $-1$ lands in $J_N^-$.

Our convention when $N$ is fingerprint-free will be that if $x \in \Z_N^*$ we’ll write $\bar{x} = \phi(x) \in C_{2B} \times C_{2^m}$ and if $\bar{f}$ is a function on $C_{2B} \times C_{2^m}$, we’ll write $f = \bar{f}\phi$ for the composition whose domain is $\Z_N^*$. So you can generally think of $\bar{\triangle}$ as the “essential version” of $\triangle$, whether $\triangle$ is an element or a function.

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

In what follows, let $\bar{g} \in C_{2B} \times C_{2^m}$ be a fixed semigenerator.

Because the $C_p$ and $C_q$ components of an element leave its Jacobi symbol alone, $\Jacobi_N(x) = (-1)^{a+c}$ where $\log_{\bar{g}}(\phi(x)) = (a, c)$. So $\phi$ carries $J_N^-$ onto the *odd-parity subset* of the value group:

```math
\begin{aligned}
S = \set{\, \bar{x} \in C_{2B} \times C_{2^m}
\st \log_{\bar{g}}(\bar{x}) = (a, c) ~\text{with}~ a + c ~\text{odd} \,}
\end{aligned}
```

Client secrets live in $J_N^-$, so $S$ is the only part of the value group the anonymity theorem has to speak about.

**Definition.** The *“essential HyperLogLog function”* maps each value in $S$ to its HyperLogLog sample value. Write $\log_{\bar{g}}(\bar{x}) = (a, c)$; split $a$ by the Chinese Remainder Theorem into $(\alpha, \beta) \in \Z_4 \times \Z_{B/2}$, which is possible because $2B = 4 \cdot \nfrac{B}{2}$ with $\nfrac{B}{2}$ odd; and when $\tz(c) ≤ m-2$ write $c = 2^{\tz(c)}u$ with $u$ odd. Then

```math
\begin{gathered}
\bar{\hll}: S \to B \times m \\[0.5em]
\bar{\hll}(\bar{x}) = \left(\beta + \tfrac{B}{2}s,~ \tz(c)\right) \\[0.8em]
s = \begin{cases}
0 & \tz(c) ≥ m-1 \\
1 & \alpha u = 2, 3 \bmod 4 \\
0 & \text{otherwise}
\end{cases}
\end{gathered}
```

The bucket index splits into two halves. The first, $\beta$, is the client’s position in the odd part of $C_{2B}$ — the direct descendant of the old bucket coordinate. The second, $s$, is a single bit read off the $2$-part. Neither factor of $s$ means anything alone, since $\alpha$ and $u$ are each pinned only up to an odd scalar; their product mod $4$ is another matter, because an odd $t$ has $t^2 = 1 \bmod 4$. That is what lets $s$ survive re-randomization, and it is why the bit must be *declared* part of the HyperLogLog value: the ring holder can read it whether or not we name it, so naming it is exactly what keeps it from being a fingerprint. The two top geometric levels are the exception — there the two orbits merge, $s$ is pinned to $0$, and only $\nfrac{B}{2}$ buckets occur, which affects a $2^{-(m-1)}$ share of clients.

As before, $\beta$ and $s$ depend on the choice of semigenerator $\bar{g}$ whereas $\tz(c)$ does not, and as before that choice only permutes bucket labels.

The HyperLogLog function for fingerprint-free $N$ on $\Z_N^*$ is defined as $\hll = \bar{\hll}\phi$, *i.e.* the composition of the $\phi$ whose existence is guaranteed by fingerprint-freeness with the essential HyperLogLog function. This depends on the choice of $\bar{g}$ for $\bar\hll$ and on which particular $\phi$ is chosen. For the purposes of the main proof, we can just assume that some fixed $\phi$ is chosen and used. The choice of $\phi$ doesn’t actually introduce any more ambiguity than already introduced by the choice of $\bar{g}$ — both choices merely permute the output bucket indices.

### The anonymity theorem

**Lemma.** For $\bar{x}, \bar{y} \in S$ we have $\bar{\hll}(\bar{x}) = \bar{\hll}(\bar{y})$ if and only if there exists $t \in \Z$ with $t = 1 \bmod{B}$ such that $\bar{x}^t = \bar{y}$.

**Proof.**  Write the logarithms of $\bar{x}$ and $\bar{y}$ as $(a_1, c_1)$ and $(a_2, c_2)$, with $a_i$ split into $(\alpha_i, \beta_i) \in \Z_4 \times \Z_{B/2}$ and, where it is defined, $c_i = 2^ku_i$ with $u_i$ odd.

Suppose first that $\bar{x}^t = \bar{y}$ for some $t = 1 \bmod B$. Since $B$ is even, $t$ is odd; since $\nfrac{B}{2}$ divides $B$, also $t = 1 \bmod \nfrac{B}{2}$. So $\beta_2 = t\beta_1 = \beta_1$, and $\tz(c_2) = \tz(tc_1) = \tz(c_1)$ because $t$ is odd. Write $k$ for that common value. If $k ≤ m-2$ then $\alpha_2 = t\alpha_1$ and $u_2 = tu_1$, so

```math
\begin{aligned}
\alpha_2u_2 = t^2\alpha_1u_1 = \alpha_1u_1 \pmod 4
\end{aligned}
```

since $t^2 = 1 \bmod 4$ for odd $t$. Hence $s$ agrees too, and $\bar{\hll}(\bar{x}) = \bar{\hll}(\bar{y})$.

Now suppose instead that $\bar{\hll}(\bar{x}) = \bar{\hll}(\bar{y})$, which gives $\beta_1 = \beta_2$, a common $k = \tz(c_1) = \tz(c_2)$, and $s_1 = s_2$. We build a suitable $t$ in two steps: first an odd $v$ acting correctly on the $2$-part, then $t$ itself by the Chinese Remainder Theorem.

When $k ≤ m-2$, take $v$ to be any odd lift of $u_2u_1^{-1} \bmod{2^{m-k}}$, so that $vc_1 = 2^k(vu_1) = 2^ku_2 = c_2$. Since $m-k ≥ 2$ this also fixes $v = u_2u_1^{-1} \bmod 4$. From $s_1 = s_2$ we have $\alpha_1u_1 = \alpha_2u_2 \bmod 4$, hence $\alpha_2 = \alpha_1u_1u_2^{-1}$, while $v\alpha_1 = \alpha_1u_2u_1^{-1}$. These agree because $u_1^2 = u_2^2 = 1 \bmod 4$ for odd $u_i$, so $v\alpha_1 = \alpha_2 \bmod 4$.

When $k ≥ m-1$ the situation is simpler: the only elements of $\Z_{2^m}$ with $m-1$ or $m$ trailing zeros are $2^{m-1}$ and $0$, so $c_1 = c_2$ outright, and $vc_1 = c_2$ for *any* odd $v$. Both $\alpha_i$ are odd, since $c_i$ is even and $a_i + c_i$ is odd, so $v = \alpha_2\alpha_1^{-1} \bmod 4$ is odd and gives $v\alpha_1 = \alpha_2$.

Either way we have an odd $v$ with $vc_1 = c_2 \bmod{2^m}$ and $v\alpha_1 = \alpha_2 \bmod 4$. Apply the Chinese Remainder Theorem — the moduli are coprime because $\nfrac{B}{2}$ is odd — to find $t$ with:

```math
\begin{aligned}
t &= 1 && \pmod{\nfrac{B}{2}} \\
t &= v && \pmod{2^m} \\
\end{aligned}
```

Then $t$ is odd, and being $1$ modulo both $2$ and $\nfrac{B}{2}$ it is $1 \bmod B$ as required. Since $m ≥ 2$ we have $t = v \bmod 4$, so $t\alpha_1 = \alpha_2$; likewise $tc_1 = vc_1 = c_2$ and $t\beta_1 = \beta_1 = \beta_2$. Reassembling $a$ from $(\alpha, \beta)$ gives $\bar{x}^t = \bar{y}$. $\square$

**Theorem.** Let $N$ be fingerprint-free. For $x, y \in J_N^-$: $\hll(x) = \hll(y)$ if and only if there exists $w \in W_N$ and $t \in \Z$ with $t = 1 \bmod{B}$ such that $wx^t = y \bmod N$.

**Proof.**  Since $N$ is fingerprint-free there exists $\phi: \Z_N^* \to C_{2B} \times C_{2^m}$ with $\ker(\phi) = W_N$, and $\phi$ carries $J_N^-$ into $S$. The following conditions are equivalent:

1. ``\hll(x) = \hll(y)``
2. ``\bar{\hll}(\phi(x)) = \bar{\hll}(\phi(y))``
3. ``\E\, t = 1 \bmod{B}`` such that $\phi(x)^t = \phi(y)$
4. ``{} \E\, t = 1 \bmod{B},\, w \in W_N {}`` such that $wx^t = y$

Equivalence of (1) and (2) is just the definition of $\hll = \bar{\hll}\phi$. The lemma we just proved gives equivalence of (2) and (3), and applies because $\phi(x)$ and $\phi(y)$ both lie in $S$. For (3) implying (4), set $w = yx^{-t}$ and observe

```math
\begin{aligned}
\phi(w) = \phi(y x^{-t}) = \phi(y)\phi(x)^{-t} = 1
\end{aligned}
```

so that $w \in \ker(\phi) = W_N$. For (4) implying (3), any $w \in W_N = \ker(\phi)$ gives

```math
\begin{aligned}
\phi(y) = \phi(wx^t) = \phi(w)\phi(x)^t = \phi(x)^t
\end{aligned}
```

This proves equivalence of (3) and (4). $\square$

Because $\ker(\phi)$ is $W_N$ on the nose, rather than $W_N$ only after intersecting with $J_N^+$, the last step needs no bookkeeping about Jacobi symbols: membership in the kernel *is* membership in the noise subgroup. That is the concrete dividend of defining fingerprint-freedom on $\Z_N^*$ rather than on $J_N^+$.

### Interpretation

How does this result prove that our scheme preserves client anonymity? Suppose there are two clients with $x_1$ and $x_2$ as their respective client secrets and assume that $\hll(x_1) = \hll(x_2)$. We also assume that the clients follow the protocol so that $\Jacobi_N(x_1) = \Jacobi_N(x_2) = -1$. The server receives $w_1x_1^{t_1}$ from the first client. Can the server tell that it got the message from the first client rather than the second? We can apply the theorem twice to show that it cannot. First, the theorem tells us that when the client chooses $w_1 \in W_N$ and $t_1$ with $t_1 = 1 \bmod{B}$, we always have:

```math
\begin{aligned}
\hll(w_1 x_1^{t_1}) = \hll(x_1)
\end{aligned}
```

In other words, the value the client sends has the same $\hll$-value as its secret. Next, since $\hll(w_1 x_1^{t_1}) = \hll(x_1) = \hll(x_2)$, the theorem tells us that there exists ${} w_2 \in W_N$ and $t_2 \in \Z$ with $t_2 = 1 \bmod{B}$ such that

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

If a malicious client sends $y$ with $\Jacobi_N(y) ≠ -1$, the server will detect it. So we can focus on a malicious client sending $y$ with $\Jacobi_N(y) = -1$. As discussed in the previous section, every element $y \in J_N^-$ is of the form $y = x_0 g^h$ where $g$ is our chosen semigenerator and $x_0$ is an arbitrary twist element with $\Jacobi_N(x_0) = -1$. Here we use $h$ as just an arbitrary attacker-controlled exponent value, not necessarily the output of a hash function. The question is whether a malicious client can influence $y$ to have large geometric sample values.

The server doesn’t publish any $x_0$ value, so clients have to generate one for themselves. However, since they don’t know the factorization of $N$, they have no idea what its logarithms are. That last claim is doing a great deal of work, and it is a property the shape of $N$ has to earn rather than one that comes for free; we return to it at the end of this section. For our analysis, write $\log(x_0) = (a, b, c, d)$, but keep in mind that the attacker has no idea what these values are. The HyperLogLog sample that $y = x_0 g^h$ encodes is $(b + h, \tz(c + h))$. The attacker controls $h$ but doesn’t know $b$ or $c$. Since they don’t know $c$, they cannot force $\tz(c + h)$ to be large or know how large it is for any particular $h$. In order to hit $k ≥ \tz(c + h)$ the attacker needs to happen to choose $h$ whose last $k$ bits complement $c$ perfectly, which occurs with probability $1/2^k$ — exactly the probability of picking a value that good by chance. What they can do, however, is scan through consecutive $h$ values. If they scan $2^k$ consecutive values, they’re guaranteed to hit $h = -c \bmod{2^k}$ for one of those values and they may happen to hit something better.

In the unique ID scheme, each request with a freshly “forged” client ID inflates the client count by one. How does HLL over RSA compare to this when an attacker tries to inflate its estimates? The expected maximum value of $\tz(c + h)$ is $k+1$ when scanning a consecutive block of $2^k$ values. This is slightly better than sending random values. If an attacker sends $B2^k$ spoofed requests they can push each bucket up to an expected value of $k + 1$, which gives a client count estimate of:

```math
\begin{aligned}
\hat{n} = \frac{1}{2 \ln(2)} B 2^{k+1} = \frac{1}{\ln(2)} B2^k
\end{aligned}
```

This uses the original Flajolet *et al.* (2007) estimator for $B ≥ 128$. Better estimators exist for small samples, and should be used for client count estimates, but for this analysis the original asymptotic formula is fine. So, at scale, an attacker can expect to inflate the estimate by $\ln(2)^{-1} \approx 1.44$ per malicious request. This is a rather good result: attack effort is linear and the coefficient is just a bit larger than one—only slightly worse than the unique ID gold standard. In short, whereas the naive unique ID scheme has perfect accuracy, resists inflation well, but has awful privacy, our HLL over RSA scheme is approximate (with tunable accuracy), recovers nearly the same inflation resistance, and provides provable client anonymity.

There is an important assumption hiding inside this linear bound: that the attacker cannot tell which of their forged tokens decoded to large geometric values. The bound holds precisely because, unable to decode, the attacker can’t cherry-pick—every forged request is a fresh blind sample, so pushing the count up takes linear effort. But the whole purpose of the system is to *publish counts*, and a published count is exactly the oracle the attacker otherwise lacks. Inject a token into an otherwise-sparse slice, watch whether the reported cardinality jumps, and you learn whether that token happened to decode high—without breaking any encryption. With that feedback an attacker can curate a collection of high-$k$ tokens, one per bucket, and a curated set is no longer linear: roughly $B$ requests can then push every bucket to level $k$, inflating the estimate to $\sim \frac{1}{\ln 2} B 2^k$, the same exponential leverage as the plaintext scheme. So the linear-resistance result holds specifically against an attacker who *cannot observe per-slice counts*. The mitigation is the minimum-cardinality floor introduced for anonymity ([Resource class sharding](03-counting-users.md#Resource-class-sharding) and the [Protocol Summary](06-protocol-summary.md)): if sub-threshold slices are never reported, the attacker has no sparse slice to probe and no decoding oracle to build. Pleasingly, the same mechanism that protects anonymity in small populations also closes this inflation channel—it is one more reason to enforce the floor, not a separate defense.

It’s worth flagging a stronger mitigation that we have deliberately *not* adopted, since it’s a natural question. One could try to *bind each token to its request*—mixing something request-specific (a server-issued nonce, or the full request path) into the value the client sends—so that a single high-value forgery can’t be reused across many requests. The difficulty is that honest counting needs the opposite: a client’s value must be *stable* across its requests within a class, or we’d be counting requests instead of unique clients. A binding that preserved per-(client, class) stability while denying an attacker reuse of a curated value is not obviously achievable, and in any case the cardinality floor already removes the count oracle that the curated-set attack depends on. So we rely on the floor and leave token-request binding as an open question—to revisit only if some deployment leaks count feedback through a channel the floor can’t cover.

While this is quite a good result, the whole bound rests on a single structural assumption, and it is worth stating that assumption sharply because the shape of $N$ has to earn it.

Everything above depends on the attacker not being able to *name* an element of $J_N^-$ — to exhibit some $z$ with $\Jacobi_N(z) = -1$ together with its $C_{2^m}$ logarithm relative to $g$. Given such a $z$ with logarithm $c$, forging is immediate rather than lucky: send $zg^h$ with $h$ chosen to cancel $c$ to whatever depth is wanted, and sweep $h$ across residues mod $\nfrac{B}{2}$ to cover every bucket. So the integrity of the scheme is exactly the claim that $J_N^-$ contains nothing anyone can name.

With $N = 3 \bmod 4$ that claim is false, and it is false for the most public element there is. It is a standard identity that $\Jacobi_N(-1) = (-1)^{(N-1)/2}$, so $N = 3 \bmod 4$ puts $-1$ in $J_N^-$; and the logarithms of $-1$ are not merely learnable but structural, since $-1$ is the unique element of order two: trivial in the bucket component and $2^{m-1}$ in $C_{2^m}$. Worse, exploiting it requires knowing nothing about $g$. Because $2^{m-1}u = 2^{m-1} \bmod{2^m}$ for every odd $u$, taking $h = 2^{m-1} \bmod{2^m}$ zeroes the geometric coordinate of $-g^h$ whatever $g$'s own logarithm happens to be. Each forgery costs one modular exponentiation, and a few thousand of them peg every bucket at the maximum geometric value. Nor does the cardinality floor help: the attacker needs no oracle, because the value is chosen open-loop rather than discovered by probing. That is the full plaintext-HyperLogLog attack, and it means the linear bound above holds only against an attacker who declines to use $-1$.

This is why $B$ is required to be $2 \bmod 4$ rather than odd. It makes $P = 2Bp+1 = 5 \bmod 8$ and hence $N = 5 \bmod 8$, so $\Jacobi_N(-1) = +1$. The elements whose logarithms are public are exactly $±\gen{g}$ — generated by the two things everybody has, the semigenerator and $-1$ — and since $\Jacobi_N(-1) = \Jacobi_N(g) = +1$, every one of them lands in $J_N^+$, where the server’s existing Jacobi check already refuses it. The bootstrapping problem the attacker faces is then real rather than assumed.

It is not, however, a proof, and the residual assumption should be named. Since $N = 5 \bmod 8$ makes both $P$ and $Q$ congruent to $1 \bmod 4$, $-1$ is a quadratic residue mod $N$ and fourth roots of unity exist. Such an $i$ has logarithms that are known outright — trivial in the bucket component and $2^{m-2}$ in $C_{2^m}$ — and $\Jacobi_N(i) = -1$ whenever $m ≥ 3$. So $i$ *is* a nameable element of $J_N^-$, and anyone who can compute one has broken the scheme. No method better than factoring $N$ is known for computing it. But finding a square root of $-1$ modulo $N$ is not known to be *equivalent* to factoring: an oracle handing back one of the four roots does not obviously yield a second, and it is the difference of two roots that factors. The honest statement is therefore a chain,

```math
\begin{aligned}
\text{factor}~N
\implies \text{compute}~i
\implies \text{forge arbitrary HyperLogLog values}
\end{aligned}
```

with no reduction known in either reverse direction, leaving the protocol’s integrity somewhere below both. Two consequences are worth drawing out. First, a single $i$, once found, works forever and for every client of that ring, which is an argument for rotating $N$ on a schedule rather than treating a modulus as permanent. Second, the certificate publishes square roots of hash-derived elements of $J_N^+$, and $-1$ now lies in $J_N^+$; if two rooted challenges were ever negatives of one another, dividing their published roots would hand out $i$. Finding such a relation among random challenges means computing the relation lattice of the rooted values, which is a discrete logarithm problem mod $N$ — the same barrier as before, and the reason the challenges must be hash-derived rather than chosen.

This is a vaguer argument than I would like to make in a write up full of rigorous proofs, and it is worth being precise about what improved and what did not. What improved is not the *form* of the argument — it is still an unproven hardness assumption — but its content: the object an attacker needs has gone from one every client already possesses to one nobody knows how to compute. Fortunately, proving that client anonymity is preserved is the much more important result, and for that we do have a solid proof. Inability to provably guarantee resistance to malicious client attacks is more acceptable: that’s a risk that we, as the operators of the system, can choose to take on. If our estimates suddenly start looking ridiculously large, we can begin to wonder if someone has cracked the protocol.

## Malicious servers

So far we’ve worried about malicious clients but have assumed that servers behave as intended. Servers are supposed to generate $N$ with the following arithmetic structure:

```math
\begin{aligned}
N = P Q = (2 B p + 1)(2^m q + 1)
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

In other words, $P-1$ is divisible by $2^m$ instead of just $2$ and $Q-1$ is divisible by $B$ when it shouldn’t be. With this alternative structure, $\Z_N^*$ has the following richer multiplicative structure:

```math
\begin{aligned}
\Z_N^* \cong
C_{2^m} \times C_B \times C_p \times
C_{2^m} \times C_B \times C_q
\end{aligned}
```

When the client follows the protocol and picks a random $x_0 \in J_N^-$, it has two random $C_B$ components and two random $C_{2^m}$ components instead of one of each. How much extra identifying information about each client does the malicious server get from this?

First, consider the two $C_B$ components. When the client raises their $x$ value to $t = 1 \bmod{B}$ both $C_B$ components are preserved. Instead of getting $\log_2(B)$ bits of identifying information from them, the server gets $2 \log_2(B)$ bits. When $B \approx 2^{12}$ that’s 24 bits of fingerprint. Multiplying $x^t$ by a white noise value, $w \in W = (\Z_N^*)^{(B/2)2^m}$, also doesn’t touch either of the two $C_B$ and $C_{2^m}$ components—the definition of $W_N$ is specifically crafted to leave those components intact while randomizing everything else.

Next, consider the two $C_{2^m}$ components. When there is only one $C_{2^m}$ component, taking $x^t$ with random odd $t$ destroys all information in it except for the position of the last bit, which is only $\log_2(m)$ bits of information. When there are two $C_{2^m}$ components, however, the information provided doesn’t just double: since the two $C_{2^m}$ components are scaled in lock-step, their _ratio_ remains fixed, which carries significantly more information. It’s not hard to see that two $C_{2^m}$ components raised to random $t$, convey a full $m$ bits of client fingerprint. Again, multiplying by $w$ doesn’t affect this, by design. So the malicious server gets an additional $m$ bits of fingerprint from the two $C_{2^m}$ components.

With this structure of $N$ then, the malicious server gets a total of $2 \log_2(B) + m$ bits of client fingerprint. For our usual choice of $B$ and $m$, that’s $2\cdot12 + 63 = 87$ bits, which is enough to uniquely identify billions of clients with random fingerprints without significant chance of collisions. If we left this unmitigated, it would be a serious defect in our protocol.

Can we convince a client that we aren’t smuggling fingerprint bits in the structure of $N$ without giving away its factorization? In this particular case, the unusual structure of $N$ is actually quite easy to detect:

- The correct structure has
  - ``P = 5 \bmod 8 \and Q = 1 \bmod 8 \implies N = PQ = 5 \bmod 8``
  - ``P = 1 \bmod B \and \gcd(\nfrac{B}{2}, Q-1) = 1 \implies \gcd(B, N-1) = 2``
- This incorrect structure has
  - ``N = P = Q = 1 \bmod 8``
  - ``\gcd(B, N-1) = B``

Unfortunately, not all possible malicious structures are so easy to detect. For example, if $N = PQR$ where $P = Q = 1 \bmod 8$ and $P = Q = 1 \bmod B$ and $R = 5 \bmod 8$ and $\gcd(B, R-1) = 2$, then you’d have $N = PQR = 5 \bmod 8$ and $\gcd(B, N-1) = \gcd(B, R-1) = 2$ so $N$ looks normal from simple modular criteria, yet $PQ$ carries the $2 \log_2(B) + m$ bits of client fingerprint from before. To guarantee that a server cannot fingerprint clients, more evidence about the structure of $N$ needs to be provided.

Zero-knowledge proofs (ZKPs) are a popular solution to this kind of problem. I spent a good bit of time going down this rabbit hole. It’s definitely doable. However, every time I sat down to implement zero-knowledge proofs of semiprimality, I found myself getting bogged down in complex and fussy details. This isn’t just a matter of laziness—if the code is that hard to implement, I find it hard to convince myself that it’s fully correct. And if we’re not confident in the correctness of the code that checks whether $N$ has the right structure, then we haven’t really proven anything.

I found myself really wishing for a simpler way to demonstrate the structure of $\Z_N^*$. One that would be easier to have confidence in, both in the sense of having a simpler implementation to review, but also in the sense of having mathematics that’s easier for more people (including myself) to understand. After some research, I found some straightforward evidence that the server can provide to convince clients a server-provided modulus cannot be used to fingerprint clients.

### Criteria for fingerprint-freedom

Recall from the section with the formal anonymity proof that an odd integer, $N$, is fingerprint-free when the quotient $\Z_N^*/W_N$ embeds into $C_{2B} \times C_{2^m}$. The subgroups of $\Z_N^*$ involved have the following definitions:

```math
\begin{aligned}
J_N^+
&= \set{\, x \in \Z_N^* \st \Jacobi_N(x) = 1 \,}
\\[0.5em]
W_N
&= (\Z_N^*)^{(B/2)2^m}
= \set{\, x^{(B/2)2^m} \st x \in \Z_N^* \,}
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
- ``\gcd(B, N) = 1 \and \gcd(B, N-1) = 2``
- For all $x, y \in J_N^+$ one of $\set{x, y, xy}$ is a quadratic residue

then $N$ is fingerprint-free.

**Proof.**  From the prior lemma we know that the last condition is equivalent to $N$ having at most two distinct prime factors. Write $N = P^jQ^k$ with $P$ and $Q$ distinct odd primes, $j ≥ 1$ and $k ≥ 0$, where $k = 0$ covers the single-prime case by making the second factor below trivial. Put

```math
\begin{aligned}
d_P = (P-1)P^{j-1}
\hspace{4em}
d_Q = (Q-1)Q^{k-1}
\end{aligned}
```

so that $\Z_N^* \cong C_{d_P} \times C_{d_Q}$, each factor cyclic because $P$ and $Q$ are odd primes. Raising a cyclic group to a power leaves a cyclic quotient whose order is a greatest common divisor:

```math
\begin{aligned}
C_M/(C_M)^E \cong C_{\gcd(M, E)}
\end{aligned}
```

Applying this componentwise with $E = (B/2)2^m$ gives

```math
\begin{aligned}
\Z_N^*/W_N \cong C_{\gcd(d_P, E)} \times C_{\gcd(d_Q, E)}
\end{aligned}
```

and it remains to show that this embeds into $C_{2B} \times C_{2^m}$. Take the $2$-parts and the odd parts in turn.

For the $2$-parts, write $e_P = \tz(P-1)$ and $e_Q = \tz(Q-1)$; these are the $2$-adic valuations of $d_P$ and $d_Q$, since $P^{j-1}$ and $Q^{k-1}$ are odd. As $\nfrac{B}{2}$ is odd, $E$ has $2$-adic valuation exactly $m$, so the $2$-part of the quotient is

```math
\begin{aligned}
C_{2^{\min(e_P,\, m)}} \times C_{2^{\min(e_Q,\, m)}}
\end{aligned}
```

If $e_P ≥ 3$ and $e_Q ≥ 3$ both held, then $P = Q = 1 \bmod 8$ and hence $N = 1 \bmod 8$, contradicting $N = 5 \bmod 8$. So at least one of them — say $e_P$, without loss of generality — is at most $2$. Then $C_{2^{\min(e_P, m)}}$ embeds into $C_4$ and $C_{2^{\min(e_Q, m)}}$ embeds into $C_{2^m}$, so the $2$-part embeds into $C_4 \times C_{2^m}$.

This is the step that the modular condition on $N$ exists to supply, and it is worth noting how little has changed from the version of this protocol with odd $B$. There, $N = 3 \bmod 4$ forced one prime to be $3 \bmod 4$ and so pinned its $2$-part to exactly $C_2$. Here $N = 5 \bmod 8$ forces one prime to be $5$ or $3$ or $7 \bmod 8$ and so pins its $2$-part to $C_4$ or smaller. Either way a single congruence on $N$, free for a client to check, is what prevents two tall $2$-parts from coexisting.

For the odd parts, let $U$ and $V$ be the odd parts of $d_P$ and $d_Q$; since the odd part of $E$ is $\nfrac{B}{2}$, the odd part of the quotient is $C_{\gcd(U,\, B/2)} \times C_{\gcd(V,\, B/2)}$. Suppose some odd prime $r$ divided both gcds. Then $r$ divides $\nfrac{B}{2}$ and also both $U$ and $V$; since $\gcd(B, N) = 1$ keeps $r$ from dividing $P$ or $Q$, it must divide both $P-1$ and $Q-1$, and hence $N-1$ by the Fact above. But then $r$ divides $\gcd(B, N-1) = 2$, which no odd prime does. So the two gcds are coprime, their product divides $\nfrac{B}{2}$, and the odd part embeds into $C_{B/2}$.

A finite abelian group is the direct product of its $2$-part and its odd part, so putting the two halves together:

```math
\begin{aligned}
\Z_N^*/W_N
\hookrightarrow C_4 \times C_{2^m} \times C_{B/2}
\cong C_{2B} \times C_{2^m}
\end{aligned}
```

where the last isomorphism is $C_4 \times C_{B/2} \cong C_{4 \cdot B/2} = C_{2B}$, valid because $\nfrac{B}{2}$ is odd. Thus $N$ is fingerprint-free. $\square$

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

Based on these results, we can design a protocol for a server to convince clients that $N$ is fingerprint-free. First, the client checks that $N = 5 \bmod 8$, that $\gcd(B, N) = 1$, and that $\gcd(B, N-1) = 2$ — the last being $\gcd(\nfrac{B}{2}, N-1) = 1$ written without the halving, since $N-1$ is even and $\nfrac{B}{2}$ is odd. These are simple numerical checks. The client is then ready to be convinced that $N$ has at most two prime divisors. The interactive version is:

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
- ``m`` — the maximum geometric sample value
- ``N`` — the ring modulus
- ``g`` — a server-selected semigenerator for $\Z_N^*$
- ``\text{sqrts}`` — a list of square roots

When downloading a new ring structure, a client checks the following requirements based on the data in this certificate:

- ``N = 5 \bmod 8``
- ``\gcd(B, N) = 1``
- ``\gcd(B, N-1) = 2``
- That enough square roots are provided
- That all the square roots are valid

That’s it. Once the client has done this, it can safely use $N$ and proceed with the protocol for generating and sending $y = xw^t$ values with each request, confident that it has provable anonymity.

### Request bundles and cross-class correlation

Sharding stops a log-holder from following a client from one resource class to another: the per-class HLL values are independent, and nothing in a client’s value for package A predicts its value for package B. That guarantee is real, but it has a boundary worth stating plainly, because requests do not arrive one at a time.

A single `Pkg.add` or `Pkg.update` resolves into a *bundle* of requests — the registry, the target package, and every dependency in its tree — issued together, from one client, in a tight window. An observer with request logs can group such a bundle by timing (and, before it is dropped, by source address) with no help from the protocol. And because a project’s dependency set is stable, the *same* client re-issues an *overlapping* bundle every time it resolves, updates, or runs CI against that project.

Once a bundle is grouped, the sharding no longer hides much. The client’s decoded values across the bundle — $(b_1, k_1), (b_2, k_2), \dots, (b_n, k_n)$ for the $n$ co-requested classes — form a tuple that is *reproducible* (a deterministic function of the client’s master key and the class set) and, for even a modest $n$, *unique*. The bucket carries the weight: it is close to uniform on $\Z_B$, worth $\log_2 B \approx 12$ bits per class, against the geometric sample’s mere $2$ bits (its distribution $P(k) = 2^{-(k+1)}$ has exactly two bits of entropy, concentrated near zero). Six co-requested packages already give $6 \cdot (12 + 2) = 84$ bits — enough to single out one client among every human on earth, many times over. Match that tuple against a later bundle and the client is re-identified across time, despite the sharding.

Two things must line up for this to enable *tracking* rather than a one-shot fingerprint. The bundle must be groupable — timing and dependency structure supply that. And the client’s bundles must *recur with overlap* over time — a client that issued a fresh, disjoint set of classes on every occasion could not be matched to itself, since its tuples would share no coordinates to compare. This is the criterion that decides whether the attack touches a given deployment: **it bites exactly when the same client re-issues overlapping bundles of sharded requests over time.** Julia’s Pkg client meets it squarely — a project’s dependency tree is requested again and again — which is why Julia does something about it. An application whose sharding guarantees that a client’s bundles do not consistently overlap over time is free of it and can leave the bucket sharded.

### Semisharding

The bundle fingerprint is dominated, overwhelmingly, by the bucket term — $n \cdot \log_2 B$ against the geometric’s $2n$ — so the natural defense is to stop the bucket from varying across classes at all. That is exactly what [semisharding](04-hyperloglog-over-rsa.md#Semisharding) does: deriving each class’s element with $f = g^B$ in place of $g$, a client keeps one fixed bucket $b_0$ in every class while its geometric sample still varies per class. The construction is a one-line change and costs nothing (see the linked section); here we account for what it buys and what it leaves.

It collapses the dominant term. The bucket’s contribution to the cross-class tuple drops from $n \cdot \log_2 B$ — reproduced afresh in each class — to a single shared $\log_2 B \approx 12$ bits, the same $b_0$ everywhere. A tuple that was unique at six classes ($84$ bits) becomes $\log_2 B + 2n$ bits: $24$ at six classes, growing only two bits per additional class rather than fourteen.

What remains is the geometric profile $(k_1, \dots, k_n)$, still reproducible and independent per class. It contributes $\sim 2n$ bits, but weakly: the distribution is concentrated near zero, so a typical client’s profile is mostly zeros (min-entropy about $n$ bits, not $2n$), and two clients whose samples are all zero are simply indistinguishable. A large enough bundle can still accumulate identifying power from the occasional high sample, so semisharding *weakens* the bundle fingerprint dramatically rather than eliminating it — which is why it is one layer among several, backed by the eager aggregation and short retention of [Reporting, Publishing & Retention](07-reporting-and-retention.md).

The trade has one cost worth naming: fixing the bucket makes $b_0$ a *stable* value that the decoder sees on every one of a client’s requests, in every class. That sounds like a pseudonym — but at the scale of a real deployment it is nothing of the kind, as the next section works out.

### The bucket as a cohort

A stable 12-bit value per client sounds alarming until you count. With $B = 2^{12}$ buckets and a population of $N$ clients, each bucket holds about $N/B$ clients, and $b_0$ says only *which* of $B$ cohorts a client falls in — never which client. It is the same construction that privacy-preserving cohort proposals such as Google’s FLoC, and its Topics successor, adopted deliberately as the anonymized alternative to per-user tracking: a coarse cohort label, not an identity. Our version is strictly stronger than those, because $b_0$ is never exposed in the clear — it rides *encrypted* inside the token and can be recovered only by the offline ring holder with the factorization, whereas FLoC’s cohort was readable by every website. The only observer who can even see $b_0$ is the one trusted party the whole protocol is built around.

How big is the cohort? Julia has on the order of $2$ million persistent clients ($N \approx 2^{21}$): the VS Code extension alone reports over a million installs, and if perhaps half of users are on VS Code, the total is around $2$M. (That is *human* clients; ephemeral CI runners are separable — the Pkg client can detect CI from environment variables — and in any case churn through a fresh master key every run, so they inflate raw request volume but form no stable cohort and are the wrong denominator here.) At $N \approx 2^{21}$ each bucket is a cohort of $2^{21}/2^{12} = 2^{9} = 512$ clients. The crowd scales with the population, so the demanding case is a *small* population, not a large one — the opposite of the intuition that more users means more exposure.

Narrowing within a cohort means stacking other stable attributes. Operating system is roughly stable and worth about a bit for a common OS — much more for a rare one (the FreeBSD user again). The geometric profile across a bundle adds a handful more. So a *typical* client — common OS, low geometric samples — sits in a slice of hundreds; the *tail* — a rare OS with a high geometric profile — can be narrowed much further, toward a handful, at this population. Two things bound even the tail. Only the *stable* attributes accumulate for cross-time tracking: OS persists, but version, region, and the geometric profile drift or reset, and $b_0$ itself resets at the yearly key regeneration described in [Reporting, Publishing & Retention](07-reporting-and-retention.md). And assembling the tuple at all still requires bundling a client’s requests, which dropping source addresses before decode, and aggregating logs eagerly, are designed to prevent. The residual is the familiar rare-cell problem — now with $b_0$ as one more, cohort-scale coordinate in the stack — handled where every rare-cell risk is handled: by the minimum-cardinality floor on what may be published, not by the cohort crowd alone. Indeed $512$ is already thinner than FLoC’s few-thousand-client anonymity target, so it would be a mistake to lean on the bucket crowd by itself; it is a useful layer, not the guarantee.

### Malicious semigenerators

We have assumed the server’s semigenerator $g$ is honest. A client cannot verify that it is — checking that $g$ generates $\Z_P^*$ and $\Z_Q^*$ needs the factorization — so we should ask what a malicious $g$ buys. Under semisharding the answer is: essentially nothing, for the same reason $f = g^B$ fixes the bucket.

The decoded token exposes only two coordinates, the bucket and the geometric sample. With $f = g^{B/2}$ the odd part of the bucket is $b_0$ for *any* $g$ — that component is annihilated by the $\nfrac{B}{2}$-th power whatever $g$ was — so the server cannot touch it. The only lever left is the geometric coordinate, $c_0 + B\,g_{2^m} h$: it covers $\Z_{2^m}$ uniformly exactly when $g_{2^m}$ is odd, and a malicious server can only make coverage *worse* (a $g$ with even $g_{2^m}$ shrinks it, degrading the geometric distribution). But degrading the geometric collapses it toward a constant, which *reduces* the per-client fingerprint while corrupting the count — it helps no attacker, and a server has no reason to sabotage the very counts it wants. And since a single published $g$ is common to all clients, it cannot target anyone. So a malicious common $g$ gains nothing against a semisharded client. As a bonus, the same $f = g^B$ neutralizes an attack that *does* work in the plain-sharded construction: a $g$ whose bucket component fails to generate $C_B$ would confine a client’s buckets to a coset, exposing a cross-class-stable residue — a real fingerprint that fixing the bucket removes.

The one thing $g$ can still do is act as a *tag*, if the server hands out a *different* $g$ to each client: the ring-id in the header is a hash of $(B, m, N, g)$, so a per-client $g$ yields a per-client ring-id. But this is parasitic. To serve a chosen $g$ to a chosen client, the server must already identify that client at certificate-fetch time — and if it can do that, it already holds an identifier at least as good as the tag it would plant. The tag adds real power in exactly one case: an active server that identifies clients live at fetch time and *then* strips addresses from its stored logs, using the persistent ring-id to re-link the de-identified token stream back to the fetch identity.

Two measures close even that. First, derive $g$ *deterministically* from $N$ by hashing — take the hash of $N$ mapped into $J_N^+$ by the same twist the certificate already uses — so a client can recompute the canonical $g$ and reject any certificate whose $g$ does not match. Because a hashed candidate in $J_N^+$ needs only its $C_{2^m}$ component to be a generator (the $C_B$ component is discarded by $f = g^B$, and the noise components are re-randomized per token), *half* of all hashed values qualify — independent of $B$ — so the server simply regenerates $N$ in the rare case its first candidate does not, yielding a canonical $g$ at no cost and with no room to choose. That removes $g$ as an independent tagging knob, leaving only a per-client *$N$* — the general “distinct certificate per client” problem, closed the way distinct-$N$ always is: by certificate transparency, clients confirming they all see the same ring. What a client *can* check locally — that $\Jacobi_N(g) = 1$ — is the anchor that makes the canonical derivation verifiable.
