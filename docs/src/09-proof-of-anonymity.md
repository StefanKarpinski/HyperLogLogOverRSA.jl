# Proof of Anonymity

In this section we present a formal proof that client anonymity is preserved by the protocol that we have outlined. Instead of assuming that $N = PQ = (2Bp + 1)(2^m q + 1)$, we generalize the setting somewhat. This weakens what needs to be shown to guarantee that clients cannot be fingerprinted by the server, which complicates the proof somewhat, but we will need this generality later on. We assume throughout that $B$ is odd and $m ≥ 2$. First, we'll present some definitions, then the main theorem and its proof, and finally conclude with why it guarantees client anonymity.

**Definition.** We write the *"positive Jacobi subgroup"* in $\Z_N^*$ as:

```math
\begin{aligned}
J_N^+
= \Jacobi_N^{-1}(+1)
= \set{\, x \in \Z_N^* \st \Jacobi_N(x) = +1 \,}
\end{aligned}
```

We'll write the *"negative Jacobi set"* (it's not a subgroup) as:

```math
\begin{aligned}
J_N^-
= \Jacobi_N^{-1}(-1)
= \set{\, x \in \Z_N^* \st \Jacobi_N(x) = -1 \,}
\end{aligned}
```

Since the Jacobi symbol only takes $±1$ values in $\Z_N^*$, $J_N^- = \Z_N^* \setminus J_N^+$, but it's convenient to have shorter notation. Note that the Jacobi symbol is only well-defined for odd $N$, so when we talk about $J_N^±$ it implicitly requires that $N$ is odd.

**Definition.** We define the *"white noise"* subgroup of $\Z_N^*$ as:

```math
\begin{aligned}
W_N &= (\Z_N^*)^{B2^m}
= \set{\, z^{B2^m} \st z \in \Z_N^* \,} ≤ J_N^+
\end{aligned}
```

Since $B2^m$ is even ($m ≥ 2$), every element has positive Jacobi symbol, so this is a subgroup of $J_N^+$. This subgroup is where we sample "noise" values to randomize the parts of $x$ that don't encode the HyperLogLog value. We previously described deriving individual $w$ values from random $z$ values; here we consider the entire group.

**Definition.**  We call a positive integer, $N$, *"fingerprint-free"* if it is odd and there exists a group homomorphism

```math
\begin{aligned}
\phi: \Z_N^* \to C_B \times C_{2^m}
\end{aligned}
```

such that

```math
\begin{aligned}
\ker(\phi) \cap J_N^+ &= +W_N \\
\ker(\phi) \cap J_N^- &= -W_N. \\
\end{aligned}
```

**Proposition.**  $N$ is fingerprint-free if and only if $N = 3 \bmod 4$ and there exists $\phi: J_N^+ \to C_B \times C_{2^m}$ such that $\ker(\phi) \subseteq W_N$.

**Proof.**  The forward direction is easy: restricting $\phi$ to $J_N^+$ already gives most of what is required, we only need to show that $N = 3 \bmod 4$. It's a standard identity that

```math
\begin{aligned}
\Jacobi_N(-1) = (-1)^{\frac{N-1}{2}}
\end{aligned}
```

which means that $N = 3 \bmod 4$ iff $-1 \in J_N^-$. Since $1 \in W_N$ we know that $-1 \in -W_N \subseteq J_N^-$ which implies that $N = 3 \bmod 4$.

The reverse direction is harder. First, we'll show that if $\phi$ exists with $\ker(\phi) \subseteq W_N$ then we can also find $\phi'$ with $\ker(\phi') = W_N$. Let $K = \ker(\phi)$. The existence of $\phi: J_N^+ \to C_B \times C_{2^m}$ with $K = \ker(\phi) ≤ W_N$ tells us that $J_N^+/K$ is cyclic with order dividing $B2^m$. It also gives the following subgroup chain:

```math
\begin{aligned}
K ≤ W_N ≤ J_N^+
\end{aligned}
```

This implies an inclusion of quotient groups:

```math
\begin{aligned}
W_N/K ≤ J_N^+/K
\end{aligned}
```

Moreover, there is a canonical isomorphism:

```math
\begin{aligned}
J_N^+/W_N \cong (J_N^+/K)/(W_N/K)
\end{aligned}
```

So $J_N^+/W_N$ is a quotient of a cyclic group with order dividing $B2^m$ which means it must also match that description. This means there exists a homomorphism, $\phi': J_N^+ \to C_B \times C_{2^m}$ with $\ker(\phi') = W_N$ exactly. In what follows, we'll just assume that we had $\ker(\phi) = W_N$ in the first place.

Recall that $N = 3 \bmod 4$ implies that $-1 \in J_N^-$. This allows us to extend $\phi$ to all of ${} \Z_N^*$ by $\phi(x) = \phi(-x)$ for $x \in J_N^-$. We need to check four identities to verify that this is a homomorphism:

1. ``x \in J_N^-``: $\phi(x) \phi(x^{-1}) = \phi(-x) \phi(-x^{-1}) = \phi(1) = 1$
2. ``x \in J_N^-, y \in J_N^-``: $\phi(x)\phi(y) = \phi(-x)\phi(-y) = \phi(xy)$
3. ``x \in J_N^-, y \in J_N^+``: $\phi(x)\phi(y) = \phi(-x)\phi(y) = \phi(-xy) = \phi(xy)$
4. ``x \in J_N^+, y \in J_N^-``: $\phi(x)\phi(y) = \phi(x)\phi(-y) = \phi(-xy) = \phi(xy)$

We already know that

```math
\begin{aligned}
\ker(\phi) \cap J_N^+ = \ker(\phi |_{J_N^+}) = W_N
\end{aligned}
```

Suppose $x \in -W_N$. Since $-x \in W_N \subseteq \ker(\phi)$, we have $\phi(x) = \phi(-x) = 1$ so $x$ is in the kernel of $\phi$ as well. This gives $-W_N \subseteq \ker(\phi)$. Now suppose $x \in \ker(\phi) \cap J_N^-$. Since $x \in J_N^-$ we have $\phi(x) = \phi(-x)$ by definition, and $x \in \ker(\phi)$ means $\phi(x) = 1$. Thus, $-x \in \ker(\phi) \cap J_N^+ = W_N$, so $x \in -W_N$. This gives the other inclusion, and together we have:

```math
\begin{aligned}
\ker(\phi) \cap J_N^- &= -W_N \\
\end{aligned}
```

This shows that our extension's kernel has the necessary intersections with $±W_N$. $\square$

Our convention when $N$ is fingerprint-free will be that if $x \in \Z_N^*$ we'll write $\bar{x} = \phi(x) \in C_B \times C_{2^m}$ and if $\bar{f}$ is a function on $C_B \times C_{2^m}$, we'll write $f = \bar{f}\phi$ for the composition whose domain is $\Z_N^*$. So you can generally think of $\bar{\triangle}$ as the "essential version" of $\triangle$, whether $\triangle$ is an element or a function.

**Definition.** We generalize our earlier definition of a *"semigenerator"* in a multiplicative group. Let $G = \prod_{i=1}^n C_{\alpha_i}$ be a product of cyclic groups and let $\pi_i: G \to C_{\alpha_i}$ be the canonical projection onto the $i$th component. An element $g \in G$ is called a semigenerator if the projection of $g$ onto each cyclic component is a generator for that component:

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

In what follows, let $\bar{g} \in C_B \times C_{2^m}$ be a fixed semigenerator.

**Definition.** The *"essential HyperLogLog function"* maps each value in $C_B \times C_{2^m}$ to its HyperLogLog sample value:

```math
\begin{gathered}
\bar{\hll}: C_B \times C_{2^m} \to B \times m \\[0.5em]
\bar{\hll}(\bar{x}) = (b, \tz(c)) \\[0.8em]
\text{where}~\log_{\bar{g}}(\bar{x}) = (b, c)
\end{gathered}
```

The first value, $b$, implicitly depends on the choice of semigenerator, $\bar{g}$, whereas the latter, $\tz(c)$, does not: $\tz(c)$ only depends on the multiplicative order of the $C_{2^m}$ part of $\bar{x}$, which is independent of $\bar{g}$. The higher bits of $c$ do depend on $\bar{g}$, but the position of the last bit does not.

The HyperLogLog function for fingerprint-free $N$ on $\Z_N^*$ is defined as $\hll = \bar{\hll}\phi$, *i.e.* the composition of the $\phi$ whose existence is guaranteed by fingerprint-freeness with the essential HyperLogLog function. This depends on the choice of $\bar{g}$ for $\bar\hll$ and on which particular $\phi$ is chosen. For the purposes of the main proof, we can just assume that some fixed $\phi$ is chosen and used. The choice of $\phi$ doesn't actually introduce any more ambiguity than already introduced by the choice of $\bar{g}$ — both choices merely permute the output bucket indices.

**Lemma.** For $\bar{x}, \bar{y} \in C_B \times C_{2^m}$ we have $\bar{\hll}(\bar{x}) = \bar{\hll}(\bar{y})$ if and only if there exists $t \in \Z$ with $t = 1 \bmod{2B}$ such that $\bar{x}^t = \bar{y}$.

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

The second equality implies that there exists odd $u \in \Z_{2^m}$ such that

```math
\begin{aligned}
uc_1 &= c_2 && \pmod{2^m} \\
\end{aligned}
```

We can apply the Chinese Remainder Theorem to find $t$ with:

```math
\begin{aligned}
t &= 1 && \pmod{B} \\
t &= u && \pmod{2^m} \\
\end{aligned}
```

Since $u$ is odd and ${} m ≥ 1$ we know that $t$ is odd as well. This means that $t = 1 \bmod{2B}$ as required. Now check that $\bar{x}^t = \bar{y}$:

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
t c_1 &= c_2 && \pmod{2^m} \\
\end{aligned}
```

Since $t$ is odd, the second equality implies that $\tz(c_1) = \tz(c_2)$ which means we've shown that $\bar{\hll}(\bar{x}) = \bar{\hll}(\bar{y})$. $\square$

**Theorem.** Let $N$ be fingerprint-free. For $x, y \in \Z_N^*$ with the same Jacobi symbol: $\hll(x) = \hll(y)$ if and only if there exists $w \in W_N$ and $t \in \Z$ with $t = 1 \bmod{2B}$ such that $wx^t = y \bmod N$.

**Proof.**  Since $N$ is fingerprint-free there exists $\phi: \Z_N^* \to C_B \times C_{2^m}$ with $\ker(\phi) = ±W_N$. The following conditions are equivalent for $x$ and $y$ with $\Jacobi_N(x) = \Jacobi_N(y)$:

1. ``\hll(x) = \hll(y)``
2. ``\bar{\hll}(\phi(x)) = \bar{\hll}(\phi(y))``
3. ``\E\, t = 1 \bmod{2B}`` such that $\phi(x)^t = \phi(y)$
4. ``{} \E\, t = 1 \bmod{2B},\, w \in W_N {}`` such that $wx^t = y$

Equivalence of (1) and (2) is just the definition of $\hll = \bar{\hll}\phi$. The lemma we just proved gives equivalence of (2) and (3). That leaves us to prove equivalence of (3) and (4). We will show that for any odd exponent, $t$, we have:

```math
\begin{aligned}
\phi(x)^t = \phi(y) ~\iff~
\E\, w \in W_N\!:~ wx^t = y
\end{aligned}
```

which means (3) and (4) are logically equivalent. Suppose that $\phi(x)^t = \phi(y)$. Let $w = y x^{-t}$ and check that $w \in \ker(\phi) \cap J_N^+ = W_N$:

```math
\begin{aligned}
\phi(w)
= \phi(y x^{-t})
= \phi(y) \phi(x)^{-t}
&= 1 \\
\Jacobi_N(w)
= \Jacobi_N(y) \Jacobi_N(x)^{-t}
&= 1
\end{aligned}
```

This requires $\Jacobi_N(x) = \Jacobi_N(y) \in \set{±1}$ and $t$ odd so that raising to $-t$ doesn't change the sign. In the other direction, suppose $w \in W_N \subseteq \ker(\phi)$ with $wx^t = y$. Check the required equality:

```math
\begin{aligned}
\phi(y)
= \phi(wx^t)
= \phi(w) \phi(x)^t
= \phi(x)^t
\end{aligned}
```

This proves equivalence of (3) and (4). $\square$

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

There is also the question of distribution of values: it could be possible that either client _could_ send the same value but it's much more likely for one of them to do so. This depends on how $w$ and $t$ values are chosen, which we haven't addressed yet. Fortunately, so long as they are chosen uniformly from publicly known ranges, every possible value is equally likely for each client, so there is truly no way to distinguish between clients with the same HyperLogLog value.
