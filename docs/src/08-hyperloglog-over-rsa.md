# HyperLogLog Over RSA

We now have the ability for clients to blindly sample geometric values, but for full HyperLogLog sampling we also need a uniformly distributed bucket value. Can we modify our RSA ring to include a bucket value too? We certainly can:

```math
\begin{aligned}
N = P Q = (2 B p + 1)(2^m q + 1)
\end{aligned}
```

This is the full RSA ring shape for encrypted HyperLogLog sampling. As before, $P$, $Q$, $p$, and $q$ are distinct, odd primes, $m$ is the maximum geometric sample value, and $B$ is the number of HyperLogLog buckets, which must be odd and coprime to everything else. The multiplicative structure of a HyperLogLog RSA ring is:

```math
\begin{aligned}
\Z_N^* \cong (C_2 \times C_B \times C_p) \times (C_{2^m} \times C_q)
\end{aligned}
```

The $C_B$ component encodes the encrypted bucket value, the order of the $C_{2^m}$ component encodes the geometric sample value, and the $C_2$ component is used to ensure that clients don't inflate geometric samples by raising values to even powers. The $C_p$ and $C_q$ components are what make values hard to decode—we don't actually care about the values in these components. Since this is starting to be a lot of components and we don't actually care about $C_p$ and $C_q$, we'll combine these two into a single cyclic component with order $pq$:

```math
\begin{aligned}
\Z_N^* \cong C_2 \times C_B \times C_{2^m} \times C_{pq}
\end{aligned}
```

Note that $C_a \times C_b \cong C_{ab}$ only holds if $a$ and $b$ are coprime, which is the case here since $p$ and $q$ are distinct primes.

The client randomly chooses and saves a persistent $x \in \Z_N^*$ value. For each new request, it randomly chooses $w \in (\Z_N^*)^{B2^m}$ and $t = 1 \bmod{2B}$ and sends $wx^t$. It's not too hard to intuitively see that the only meaningful pieces of information conveyed about $x$ are:

1. The parity bit
2. The bucket value
3. The geometric sample

We can analyze it in terms of generator exponents. Choose a semigenerator, $g$, in this new ring and write $\log_g(x) = (a, b, c, d)$ and $\log_g(w) = (0, 0, 0, e)$. Then:

```math
\begin{aligned}
\log_g(wx^t)
&= (ta, tb, tc, td + e) \\
&= (a, b, tc, td + e) \\
\end{aligned}
```

Spelling these four components out:

- Parity bit: $ta = a \bmod 2$ since $t = 1 \bmod 2$
- Bucket index: $tb = b$ since $t = 1 \bmod B$
- Geometric sample: $\tz(tc) = \tz(t) + \tz(c) = \tz(c)$ since $\tz(t) = 0$
- The rest: $td + e \bmod{pq}$ is freshly random in each request

The Jacobi symbol works as before to guarantee that an odd exponent was used. If the above intuitive explanation of why only the Jacobi symbol and the HyperLogLog value can be recovered from $wx^t$ is good enough for you, you can skip the next section, but if you want to see how one can formally prove the claim of anonymity, read on.
