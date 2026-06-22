# Fighting Inflation

Clients are supposed to raise $x$ to an odd power, $t$. What happens if they raise it to an even $t$ instead? For any values $t$ and ${} c$ we have $\tz(tc) = \tz(t) + \tz(c)$. When $t$ is odd we have $\tz(t) = 0$, so $\tz(tc) = \tz(c)$. But if $t$ is even then $\tz(t) > 0$ and $x^t$ inflates the geometric sample. If a client sends $x^{2^m}$, for example, then it is guaranteed to produce the maximal geometric sample value, which is supposed to be vanishingly rare, occurring with $1/2^m$ probability. If we require that $a = 1 \bmod 2$ in $x$, however, then $ta = t \bmod 2$, which lets us read the parity of $t$ form the first component of $wx^t$. Parity bit to the rescue! So we'd like to require clients to choose $x$ with odd $a$ and then check that $ta = 1$ in the $wx^t$ value that is sent. Any requests not satisfying this should be ignored for client count estimation purposes, since that indicates a malicious or malfunctioning client.

Readers may wonder how clients can check whether the $x$ that they've chosen has $a = 1 \bmod 2$. After all, the whole point of working in $\Z_N$ is that people who don't know the factorization of $N$ can't extract the components of $x \in \Z_N$. Someone who knows the factorization of $N$ can, of course, check this easily:

```math
\begin{aligned}
a = 0 ~~\iff~~ \fmod(x,P)^p = 1
\end{aligned}
```

But the client doesn't know $P$ and can't check this. It is, however, possible to efficiently compute the [Jacobi symbol](https://en.wikipedia.org/wiki/Jacobi_symbol), $\Jacobi_N(x)$, without knowing the factorization of $N$. For this particular ring shape, the value of the Jacobi symbol is $(-1)^{a + c}$. In general, the Jacobi symbol is the total parity of all the log-coordinates whose moduli are even. This doesn't tell us the value of $a$ by itself, but if $\Jacobi_N(x) = -1$ then we know that $a + c = 1 \bmod 2$ so we are in one of these two cases:

```math
\begin{alignedat}{3}
a &= 1 \bmod 2 &~~\wedge~~&& c &= 0 \bmod 2 \\
a &= 0 \bmod 2 &~~\wedge~~&& c &= 1 \bmod 2 \\
\end{alignedat}
```

This isn't exactly what we wanted, but it does give us what we need. Instead of requiring $a = 1$, we can require that $\Jacobi_N(wx^t) = -1$. This forces us into one of these two cases:

1. $ta = 1 \bmod 2$ and $tc = 0 \bmod 2$
2. $ta = 0 \bmod 2$ and $tc = 1 \bmod 2$

Whereas if $t$ is even, then $ta = tc = 0 \bmod 2$, so this check guarantees that odd $t$ was used.

Before proceeding with this, we should verify that conditioning on $\Jacobi_N(x) = -1$ does not change the relative probabilities of different $\tz(c)$ values. Fortunately it doesn't:

```math
\begin{aligned}
\mathcal{P}(\tz(c) = 0 ~|~ \Jacobi_N(x) = -1)
= \mathcal{P}(\tz(c) = 0)
= \tfrac{1}{2}
\end{aligned}
```

In summary, instead of requiring $x$ with $a = 1 \bmod 2$, which would work but cannot be checked by clients, we require $x$ with $\Jacobi_N(x) = -1$, which they can check. This implies that $\Jacobi_N(wx^t) = -1$ if and only if $t$ is odd. The server can check that this is the case and disregard any requests that don't satisfy this parity requirement.
