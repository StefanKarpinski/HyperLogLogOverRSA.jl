# Malicious Servers

So far we've worried about malicious clients but have assumed that servers behave as intended. Servers are supposed to generate $N$ with the following arithmetic structure:

```math
\begin{aligned}
N = P Q = (2 B p + 1)(2^m q + 1)
\end{aligned}
```

where $B$ and $m$ are published parameters and $P$, $Q$, $p$ and $q$ are secret primes. If a client knew $P$ and $Q$ they could check the construction of $N$, but the whole point, of course, is that they don't know the factorization. Otherwise they could easily forge rare HyperLogLog values. What if a server generates $N$ with a different structure? Can that allow them to "fingerprint" and track individual clients? In fact, constructing $N$ with different structure than intended *can* allow a server to fingerprint clients. For completeness, our protocol needs a mechanism for servers to convince clients that they have not hidden any fingerprints in the construction of $N$.

First, let's get a feel for how a malicious server could fingerprint clients. Here's an alternative structure that would allow a server to uniquely identify clients:

```math
\begin{aligned}
N = P Q = (2^m B p + 1)(2^m B q + 1)
\end{aligned}
```

In other words, $P-1$ is divisible by $2^m$ instead of just $2$ and $Q-1$ is divisible by $B$ when it shouldn't be. With this alternative structure, $\Z_N^*$ has this multiplicative structure:

```math
\begin{aligned}
\Z_N^* \cong
C_{2^m} \times C_B \times C_p \times
C_{2^m} \times C_B \times C_q
\end{aligned}
```

When the client follows the protocol and picks a random $x_0 \in J_N^-$, it has two random $C_B$ components and two random $C_{2^m}$ components instead of one of each. How much extra identifying information about each client does the malicious server get from this?

First we'll consider the two $C_B$ components. When the client raises their $x$ value to $t = 1 \bmod{2B}$ both $C_B$ components are preserved. Instead of getting $\log_2(B)$ bits of identifying information from them, the server gets $2 \log_2(B)$ bits. When $B \approx 2^{12}$ that's 24 bits of fingerprint. Multiplying $x^t$ by a white noise value, $w \in W = (\Z_N^*)^{B2^m}$, also doesn't touch either of the two $C_B$ and $C_{2^m}$ components—the definition of $W_N$ is specifically crafted to leave those components intact while randomizing everything else.

Next, we consider the two $C_{2^m}$ components. When there is only one $C_{2^m}$ component, taking $x^t$ with random odd $t$ destroys all information in it except for the position of the last bit, which is only $\log_2(m)$ bits of information. When there are two $C_{2^m}$ components, however, the information provided doesn't just double: since the two $C_{2^m}$ components are scaled in lock-step, their _ratio_ remains fixed, which carries significantly more information. It's not hard to see that two $C_{2^m}$ components raised to random $t$, convey a full $m$ bits of client fingerprint. Again, multiplying by $w$ doesn't affect this, by design. So the malicious server gets an additional $m$ bits of fingerprint from the two $C_{2^m}$ components.

With this structure of $N$, the malicious server gets a total of $2 \log_2(B) + m$ bits of client fingerprint. For our usual choice of $B$ and $m$, that's $2\cdot12 + 63 = 87$ bits, which is enough to uniquely identify billions of clients with random fingerprints without significant chance of collisions. If we left this unmitigated, it would be a serious defect in our protocol.

Can we convince a client that we aren't smuggling fingerprint bits in the structure of $N$ without giving away its factorization? In this particular case, the unusual structure of $N$ is actually quite easy to detect:

- The correct structure has
  - ``P = 3 \bmod 4 \and Q = 1 \bmod 4 \implies N = PQ = 3 \bmod 4``
  - ``P = 1 \bmod B \and Q ≠ 1 \bmod B \implies N = PQ ≠ 1 \bmod B``
- This incorrect structure has
  - ``N = P = Q = 1 \bmod 4``
  - ``N = P = Q = 1 \bmod B``

Unfortunately, not all possible malicious structures are so easy to detect. For example, if $N = PQR$ where $P = Q = 1 \bmod 4$ and $P = Q = 1 \bmod B$ and $R = 3 \bmod 4$ and $R ≠ 1 \bmod B$, then you'd have $N = PQR = 3 \bmod 4$ and $N = PQR ≠ 1 \bmod B$ so $N$ looks normal from simple modular criteria, yet $PQ$ carries the $2 \log_2(B) + m$ bits of client fingerprint from before. To guarantee that a server cannot fingerprint clients, more evidence about the structure of $N$ needs to be provided.

Zero-knowledge proofs (ZKPs) are a popular solution to this kind of problem. I spent a good bit of time going down this rabbit hole. It's definitely doable. However, every time I sat down to implement zero-knowledge proofs of semiprimality, I found myself getting bogged down in complex and fussy details. This isn't just a matter of laziness—if the code is that hard to implement, I find it hard to convince myself that it's fully correct. And if we're not confident in the correctness of the code that checks whether $N$ has the right structure, then we haven't really proven anything.

I found myself really wishing for a simpler way to demonstrate the structure of $\Z_N^*$. One that would be easier to have confidence in, both in the sense of having a simpler implementation to review, but also in the sense of having mathematics that's easier for more people (including myself) to understand. After some research, I found some straightforward evidence that the server can provide to convince clients a server-provided modulus cannot be used to fingerprint clients.

### Criteria for fingerprint-freedom

Recall from the section with the formal anonymity proof that an odd integer, $N$, is shown to be fingerprint-free if $N = 3 \bmod 4$ and $J_N^+/W_N$ is cyclic with order dividing $B2^m$. These subgroups of $\Z_N^*$ have the following definitions:

```math
\begin{aligned}
J_N^+
&= \set{\, x \in \Z_N^* \st \Jacobi_N(x) = 1 \,}
\\[0.5em]
W_N
&= (\Z_N^*)^{B2^m}
= \set{\, x^{B2^m} \st x \in \Z_N^* \,}
\end{aligned}
```

Our proof in that section justified the term "fingerprint-free" by showing that as long as $N$ satisfies this property, a server cannot distinguish between two correctly behaving clients with the same HyperLogLog values. Our task then, is to come up with a way for the server to convince clients that the value of $N$ it sends them is actually fingerprint-free, without revealing its factorization.

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

We will often prove results in terms of logarithmic coordinates. These kinds of results are not usually proven with this approach, but since we've already worked extensively with this viewpoint, it will hopefully be more illuminating in this context.

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

**Definition** (standard)**.** An element $x \in \Z_N$ is called a "quadratic residue" if there exists $r \in \Z_N$ such that $r^2 = x \bmod N$. The set of quadratic residues may be denoted as $(\Z_N)^2$.

**Lemma.** Let $N$ be an odd, positive integer. $N$ has at most two distinct prime factors if and only if for all $\set{x, y} \subseteq J_N^+$ at least one of $\set{x, y, xy}$ is a quadratic residue.

**Proof.**  Since $(p_i-1)\,p_i^{n_i-1}$ is always even, each of the log-coordinates has a well-defined notion of parity. We can understand both the Jacobi symbol and quadratic residues in terms of parity of log-coordinates:

- ``\Jacobi_N(x) = 1`` iff an even number of the $\ell_i$ are odd
- ``x`` is a quadratic residue iff none of the $\ell_i$ are odd

In the $D = 1$ case, $\Jacobi_N(x) = 1$ if and only if $x$ is a quadratic residue, so all three of $x$, $y$ and $xy$ are quadratic residues.
In the $D = 2$ case, let's write coordinates for $x$, $y$ and $xy$:

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

- ``N = 3 \bmod 4``
- ``\gcd(B, N) = \gcd(B, N-1) = 1``
- For all $x, y \in J_N^+$ one of $\set{x, y, xy}$ is a quadratic residue

then $N$ is fingerprint-free.

**Proof.**  From the prior lemma we know that the last condition is equivalent to $N$ having at most two distinct prime factors.

First we'll consider the case of a single prime factor: $N = P^j$ where $P$ is an odd prime and ${} j ≥ 1 {}$. If $P = 1 \bmod 4$ then we'd have $N = 1 \bmod 4$ as well, so we know that $P = 3 \bmod 4$. Let $d = \gcd(B, P-1)$. Since $d \divides P-1$ we have $d \divides N-1$ and therefore $d \divides \gcd(B, N-1) = 1$. So $B$ and $P-1$ must be coprime. Thus, the structure of $\Z_N^*$ is:

```math
\begin{aligned}
\Z_N^* \cong C_2 \times C_U
~~\text{where}~~
U = \frac{P-1}{2} \, P^{j-1}
\end{aligned}
```

We know that $U$ is odd and coprime to $B$ since $(P-1)/2$ and $P$ are. To prove that $N$ is fingerprint-free we need a homomorphism from $J_N^+$ to a subset of $C_B \times C_{2^m}$ whose kernel is a subset of $W_N ≤ J_N^+$. In this case, however, $W_N = J_N^+$ so the constant homomorphism, which has all of $J_N^+$ as its kernel, works. To see this, consider generic $x \in \Z_N^*$:

```math
\begin{aligned}
\log_g(x) = (a, b) \in \Z_2 \times \Z_U
\end{aligned}
```

``\Jacobi_N(x) = (-1)^a`` since $U$ is odd, so $x \in J_N^+$ if and only if $a = 0$. On the other hand, $w \in W$ if and only if $w = z^{B2^m}$ for some $z \in \Z_N^*$ or:

```math
\begin{aligned}
\log_g(w)
= (a_w, b_w)
= (a_z B2^m, b_z B2^m)
= B2^m \log_2(z)
\end{aligned}
```

This forces $a_w = 0$ since $m ≥ 2$. Does it impose any restriction on $b_w$? No, since for any $b_w$ we can let

```math
\begin{aligned}
b_z = b_w(B2^m)^{-1} \pmod U
\end{aligned}
```

This inverse exists since $B2^m$ is coprime to $U$. This shows that both $J_N^+$ and $W_N$ are precisely the subset where $a = 0 \bmod 2$ and that they are equal, and therefore the constant homomorphism witnesses that $N$ is fingerprint-free.

Next, consider the case of two distinct prime factors: $N = P^j Q^k$ where $P$ and $Q$ are distinct odd primes and $j ≥ 1$ and $k ≥ 1$. Since $N = 3 \bmod 4$, we can assume without loss of generality that $P^j = P = 3 \bmod 4$, which implies that $1 = \tz(P-1)$. Let $n = \tz(Q-1)$. Then the structure of $\Z_N^*$ is:

```math
\begin{gathered}
\Z_N^* \cong C_2 \times C_U \times C_{2^n} \times C_V \\[0.5em]
\text{where} \\
U = \frac{P-1}{2} \, P^{j-1}
\hspace{4em}
V = \frac{Q-1}{2^n} \, Q^{k-1}
\end{gathered}
```

``U`` and $V$ are both odd and coprime to $B$. They are odd since $P$, $Q$, $(P-1)/2$ and $(Q-1)/2^n$ are all odd. We know that $P$ and $Q$ are coprime to $B$ since $\gcd(B, PQ) = 1$. Let $d = \gcd(B, P-1, Q-1)$. Since $d$ divides both $P-1$ and $Q-1$ it must also divide $N-1$, but since $\gcd(B, N-1) = 1$ the only option is $d = 1$. Thus, $U$ and $V$ are also coprime to $B$.

Let $B_U = \gcd(B, U)$ and $B_V = \gcd(B, V)$. These are coprime since $\gcd(B, U, V) = 1$. Let $\bar{m} = \min(m, n) ≥ 1$ and pick a semigenerator $\bar{g} \in C_{B_U} \times C_{2^{\bar{m}}} \times C_{B_V}$. For $x \in \Z_N^*$ with $\log_g(x) = (a, b, c, d)$ define

```math
\begin{aligned}
\phi: \Z_N^* \to
C_{B_U} \times C_{2^{\bar{m}}} \times
C_{B_V}
\end{aligned}
```

```math
\begin{aligned}
\log_{\bar{g}}(\phi(x)) = (b, c, d)
\in \Z_{B_U}
\times \Z_{2^{\bar{m}}}
\times \Z_{B_V}
\end{aligned}
```

Note that in addition to dropping the first coordinate, the remaining three coordinates are modularly reduced. In order to show that $\phi$ witnesses that $N$ is fingerprint-free, we need to show that $\ker(\phi|_{J_N^+}) \subseteq W_N$. Let

```math
\begin{aligned}
w \in \ker(\phi|_{J_N^+}) = J_N^+ \cap \ker(\phi)
\end{aligned}
```

Since $w \in J_N^+$ we have:

```math
\begin{aligned}
a &= c \pmod 2 \\
\end{aligned}
```

Since $w \in \ker(\phi)$ we have:

```math
\begin{aligned}
b &= 0 \pmod{B_U} \\
c &= 0 \pmod{2^{\bar{m}}} \\
d &= 0 \pmod{B_V} \\
\end{aligned}
```

Thus, there exist $b'$, $c'$ and $d'$ such that:

```math
\begin{aligned}
b &= b' B_U &
c &= c' 2^{\bar{m}} &
d &= d' B_V
\end{aligned}
```

To show that $w \in W$ we need to find $i \in \Z$ such that:

```math
\begin{aligned}
i B 2^m &= 0 && \pmod 2 \\
i B 2^m &= b = b' B_U && \pmod U \\
i B 2^m &= c = c' 2^{\bar{m}} && \pmod{2^n} \\
i B 2^m &= d = d' B_V && \pmod V \\
\end{aligned}
```

The first equation is automatically satisfied since $m ≥ 1$. The other three equations are equivalent to:

```math
\begin{aligned}
i (B/B_U) 2^m &= b' && \pmod{U/B_U} \\
i B 2^{m-\bar{m}} &= c' && \pmod{2^{n-\bar{m}}} \\
i (B/B_V) 2^m &= d' && \pmod{V/B_V} \\
\end{aligned}
```

Here we have divided common factors—$B_U$, $2^{\bar{m}}$ and $B_V$, respectively—out of each equation and modulus. These equations are in turn equivalent to:

```math
\begin{aligned}
i &= b' (B/B_U 2^m)^{-1} && \pmod{U/B_U} \\
i &= c' B^{-1}     && \pmod{2^{n-\bar{m}}} \\
i &= d' (B/B_V 2^m)^{-1} && \pmod{V/B_V} \\
\end{aligned}
```

This set of equations can be solved via the Chinese Remainder Theorem since the moduli are pairwise coprime. Let $x = g^i$, which gives:

```math
\begin{aligned}
\log_g(x^{B2^m})
&= \log_g(g^{iB2^m}) \\
&= (iB2^m, iB2^m, iB2^m, iB2^m) \\
&= (0, b' B_U, c' 2^{\bar{m}}, d' B_V) \\
&= (0, b, c, d) \\
&= \log_g(w)
\end{aligned}
```

This shows that $w \in W_N$ as claimed, so $\ker(\phi|_{J_N^+}) \subseteq W_N$ and thus $\phi$ witnesses that $N$ is fingerprint-free in the two prime factor case. $\square$

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

**Proof.**  It's straightforward to see that the probability of $x \in J_N^+$ being a quadratic residue is $1/2^{D-1}$:

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

Based on these results, we can design a protocol for a server to convince clients that $N$ is fingerprint-free. First, the client checks that $N = 3 \bmod 4$, that $\gcd(B, N) = 1$, and that $\gcd(B, N-1) = 1$. These are simple numerical checks. The client is then ready to be convinced that $N$ has at most two prime divisors. The interactive version is:

> The client picks $n$ random pairs $\set{x, y} \subseteq J_N^+$ and challenges the server to produce $r \in \Z_N^*$ for each pair such that $r^2 \in \set{x, y, xy} \bmod N$.

This convinces the client that there's at most a $(5/8)^{n}$ probability that $N$ has more than two distinct prime factors. If the client wants to be convinced to a probability of $\nfrac{1}{\alpha}$ for some large $\alpha$, they should choose $n$ such that:

```math
\begin{alignedat}{2}
\left(\frac{5}{8}\right)^{n}\! &≤ \frac{1}{\alpha}
~~&\iff~
n &≥ \frac{\log_2(\alpha)}{\log_2(8)-\log_2(5)}
\end{alignedat}
```

We can turn this into a non-interactive protocol by picking a cryptographic hashing scheme to generate pairs of values in such a way that the server cannot influence them. The hash should include $N$ as an input so that the set of hash-generated values is not fixed and must be computed after choosing a candidate modulus. This prevents generating $n$ pairs and then solving for an $N$ that happens to have quadratic residues for those pairs. In order to generate an incorrectly structured $N$ that passes this test, a malicious server would have to try about $\alpha$ candidate $N$ values before expecting to find one that passes. Since $n$ grows only as $\log_2(\alpha)$, we can demand an enormous amount of attacker work while keeping the certificate small—but we have to choose $\alpha$ with the real cost of a forgery attempt in mind, and here it pays to be pessimistic, because a successful forgery is catastrophic. The fingerprinting structures from the start of this section yield on the order of $2\log_2(B) + m \approx 87$ bits of fingerprint—enough to deanonymize every client of the system at once—so $\alpha$ is guarding against total failure, not a marginal leak.

It's tempting to reason about the per-attempt cost in terms of the $n$ quadratic-residue checks, but that reasoning cuts the wrong way: the attacker is forging *their own* candidate moduli, so they know each candidate's factorization, which makes the QR checks essentially free for them. Worse for the defender, the attacker fails fast—each challenge pair passes a bad modulus with probability at most $5/8$, so on a doomed candidate they bail after two or three pairs and pay the full $n$ square roots only on the eventual success. The real per-attempt cost is dominated by generating a constrained candidate modulus, and the search is embarrassingly parallel. We should not lean on QR-testing being expensive; for this attacker, it isn't.

To put numbers on it, take a deliberately attacker-favorable estimate of $\sim 10^{11}$ modular exponentiations per second of aggregate compute, and pretend each attempt costs only a single modexp (generating a constrained prime is much more expensive, so this *understates* the attacker's work, which is the safe direction for a security margin):

| $\alpha$ | candidates | wall-clock @ $10^{11}$/s | $n = \ceil{\log_2\alpha / \log_2(8/5)}$ | cert size @ $L=1024$ |
|---|---|---|---|---|
| $2^{50}$ | $1.1\times10^{15}$ | ~3 hours | 74 | ~9.5 KB |
| $2^{64}$ | $1.8\times10^{19}$ | ~6 years | 95 | ~12 KB |
| $2^{80}$ | $1.2\times10^{24}$ | ~$4\times10^5$ years | 118 | ~15 KB |
| $2^{112}$ | infeasible | infeasible | 166 | ~21 KB |
| $2^{128}$ | infeasible | infeasible | 189 | ~24 KB |

The lesson is that the $\alpha = 2^{50}$ I used in early test code is *not* "completely infeasible": at a few hours of attacker-favorable wall-clock it is squarely within reach of a motivated operator—which is exactly the threat we care about, a platform with both the means and an incentive to deanonymize its own users. Because $n$ grows only logarithmically in $\alpha$, buying a gigantic safety margin is nearly free: the certificate stays a few dozen kilobytes and the server's square-root computation stays sub-second. The recommendation is therefore to set $\alpha \ge 2^\lambda$ for the system's overall security level $\lambda$, and never below $2^{80}$; the simplest safe choice is $\alpha = 2^{128}$ ($n = 189$, a ~24 KB certificate), at which point the parameter can simply be forgotten about.

The non-interactive version of this protocol serves as a certificate of fingerprint-freedom for a published $N$ value. The certificate structure contains:

- ``B`` — the number of buckets
- ``m`` — the maximum geometric sample value
- ``N`` — the ring modulus
- ``g`` — a server-selected semigenerator for $\Z_N^*$
- ``\text{sqrts}`` — a list of square roots

When downloading a new ring structure, a client checks the following requirements based on the data in this certificate:

- ``N = 3 \bmod 4``
- ``\gcd(B, N) = 1``
- ``\gcd(B, N-1) = 1``
- That enough square roots are provided
- That all the square roots are valid

That's it. Once the client has done this, it can safely use $N$ and proceed with the protocol for generating and sending $xw^t$ values with each request, confident that it has provable anonymity.
