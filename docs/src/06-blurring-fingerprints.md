# Blurring Fingerprints

We now have a construction where a client can sample a value with a geometric distribution without knowing the value they've sampled. We can't, however, have a client just generate a random $x \in \Z_N^*$ and send that $x$ every time—this would be a massive, uniquely identifying client fingerprint. We need some way of destroying all identifying information about $x$ _except_ for the geometric sample that we want. Fortunately, this turns out to be entirely doable.

We want to erase all information about the $C_{2^m}$ component except for its order, effectively preserving only $\tz(c)$. Fortunately, this is actually straightforward: for any odd $t$, $x^t$ has the same geometric sample value as $x$. We can see this easily from logarithms:

```math
\begin{aligned}
\log_g(x^t) = t (a, b, c, d) = (ta, tb, tc, td)
\end{aligned}
```

Since ${} \tz(tc) = \tz(t) + \tz(c) = \tz(c)$ we know that $x^t$ and $x$ have the same geometric value. This crucially depends on $t$ being odd. The higher bits of $tc$ are fully random: given $c$ and $c'$ with $k = \tz(c) = \tz(c')$ we can find $t$ that connects them, namely $t = (c'/2^k)(c/2^k)^{-1} \bmod{2^m}$. Instead of sending $x$, the client can choose a new random, odd $t$ for each request and send $x^t$. As long as $t$ can fall into any odd residue class modulo $2^m$, the leading bits of $tc \bmod{2^m}$ are completely arbitrary and only the position of the trailing bit is preserved.

Sending $x^t$ obscures $x$'s exact value of $c$, but what about the other components? The Chinese Remainder Theorem tells us that if $b$ and $d$ are non-zero, for any $b'$ and ${} d'$ there exists $t$ such that the following modular equalities all hold:

```math
\begin{aligned}
t &= 1        &&\pmod 2 \\
t &= b'b^{-1} &&\pmod p \\
t &= d'd^{-1} &&\pmod q \\
\end{aligned}
```

This is possible since $2$, $p$ and $q$ are pairwise coprime, and gives:

```math
\begin{aligned}
\log_g(x^t) = (ta, tb, tc, td) = (a, b', tc, d')
\end{aligned}
```

We have $ta = a \bmod 2$ since $t$ is odd. In other words, we can hit any combination of values in the $C_p$ and $C_q$ components for some value of $t \in \Z_{pq2^m}$.

Being able to reach any possible values in the $C_p$ and $C_q$ components via exponentiation is predicated on $b$ and $d$ being non-zero, however. In an RSA ring with our proposed structure, if $N$ is thousands of bits (as it would be in real usage), it is astronomically unlikely to choose $x$ with zero exponents in the $C_p$ or $C_q$ components. For the sake of sheer thoroughness, however, let's erase even that tiny bit of information. Choose $z \in \Z_N^*$ at random and multiply by $w = z^{2^m}$. If $\log_g(z) = (a_{z}, b_{z}, c_{z}, d_{z})$, we have:

```math
\begin{aligned}
\log_g(w)
&= 2^m \, (a_{z}, b_{z}, c_{z}, d_{z}) \\
&= (0, b_{z} 2^m, 0, d_{z} 2^m) \\
\end{aligned}
```

The $C_2$ and $C_{2^m}$ exponents of $w$ are both zero since $2^m = 0$ in both moduli. The $C_p$ and $C_q$ parts are random and arbitrary since $b_{z}$ and $d_{z}$ are random and arbitrary and multiplication by $2^m$ just permutes those already random values. Since $p$ and $q$ are coprime, the Chinese Remainder Theorem guarantees the existence of $e \in \Z$ such that:

```math
\begin{aligned}
e = b_{z} 2^m \pmod p \\
e = d_{z} 2^m \pmod q \\
\end{aligned}
```

This lets us write the logarithm of $w$ with a single unknown:

```math
\begin{aligned}
\log_g(w) = (0, e, 0, e)
\end{aligned}
```

We can multiply by $w$ before or after exponentiating by $t$:

```math
\begin{aligned}
\log_g(wx^t) &= (a, tb + e, tc, td + e) \\
\log_g((wx)^t) &= (a, t(b + e), tc, t(d + e)) \\
\end{aligned}
```

Either way, every possible value in the $C_p$ and $C_q$ components can be reached for some value of $w$ and $t$. We'll use $wx^t$; this version has a somwhat stronger guarantee: regardless of the other values, including $t$, there is some value of $w$ that can produce any pair of values in the $C_p$ and $C_q$ components.

The final component we need to consider is $C_2$: the value of $a$ in $x$ is unchanged by both $x \mapsto x^t$ and by multiplication by $w^{2^m}$. As it turns out, however, we actually need this parity bit to to be preserved in order to prevent clients from artificially inflating their geometric samples.
