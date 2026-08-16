# HyperLogLog Over RSA

What alternative is there to signing HyperLogLog values? In the last section we encrypted HLL values. What if clients could sample encrypted HLL values for themselves? Of course, the way we discussed encrypting them previously, each distinct HLL value would have an equal chance of being chosen, which is entirely the wrong distribution—the whole point is that it’s geometric. But can we encrypt HLL values so that sampling encrypted values has the right distribution? This requires more common HLL values to appear in distinct encrypted forms many times: the most common values would have to appear $2^{m-1}$ times as often as the least common ones, where $m$ bounds the geometric range. This is certainly possible—we can just encrypt a 128-bit unsigned random value and derive the HLL value from it on the server side. But that’s a massive client fingerprint. We need another crucial ingredient: the client needs to be able to randomize what they send so that any client with the same HLL value could have sent that value. This seemingly impossible trick turns out to be possible in the mathematical setting of RSA rings.

## RSA ring refresher

An RSA ring is a modular integer ring, $\Z_N$, where $N$ is a large, balanced semiprime—*i.e.* a product of two similarly sized primes: $N = P Q$. When $N$ is large enough, this is a hard factoring problem—hard enough that it is currently impossible when $\log_2(N) ≥ 1024$. For good measure, people often use $2048$-bit and $4096$-bit semiprimes these days. Constructing such an $N$, on the other hand, is straightforward and takes no more than a fraction of a second: find two distinct primes, $P$ and $Q$, of appropriate size, and multiply them to get $N = PQ$. The party that generates $N$ can safely publish it and keep $P$ and $Q$ secret. This gives them knowledge of the structure of $\Z_N$ that no one else has, yet everyone can perform computations in $\Z_N$. This was one of the first public-key encryption schemes and is still widely used. If all you want to do is public-key cryptography, the additional structure of $\Z_N$ is actually a liability and has to be carefully tiptoed around. But for us the additional structure is crucial to making what we want to do possible.

The classic use of an RSA ring is for public-key cryptography. The party that generates $N$ chooses an exponent, $E$, and computes $D = E^{-1} \bmod{\lambda(N)}$. Here $\lambda(N)$ is the [Carmichael function](https://en.wikipedia.org/wiki/Carmichael_function), defined as the smallest exponent such that for all $x \in \Z_N$, $x^{\lambda(N)} = 1 \bmod N$. This function is crucial for understanding exponents in $\Z_N^*$ because exponents behave modularly as well, but with a modulus of $\lambda(N)$: *i.e.* $x^a = x^b \bmod N$ whenever $a = b \bmod{\lambda(N)}$. For semiprime $N = PQ$, we have $\lambda(N) = \lcm(P-1, Q-1)$, which is easy for someone who knows $P$ and $Q$ to compute, but infeasible for someone who only knows $N$. Since $DE = 1 \bmod{\lambda(N)}$ we have $x^{DE} = x^1 = x \bmod N$ for all $x \in \Z_N^*$. Thus, the exponents $E$ and $D$ form a public/private key pair:

- You can publish $E$ and anyone can encrypt a message $x \in \Z_N^*$ by exponentiating by $E$: $c = x^E \bmod N$.
- Only someone who knows $D$ can decrypt the message by exponentiating by $D$: $c^D = (x^E)^D = x^{DE} = x \bmod N$.

This system can also be used for signatures. To sign a message, $x$, the party who knows $D$ computes $s = x^D \bmod N$ and sends $(x, s)$ as a message/signature pair. Anyone with the public key can check the signature like this: 

```math
\begin{aligned}
s^E = (x^D)^E = x^{DE} = x \bmod N.
\end{aligned}
```

 If $s^E$ doesn’t match the message, the signature is invalid. Anyone can check a signature, but only someone who knows $E$ can generate a valid one, so a valid signed message proves that the message is authentically from the party who knows $E$.

We will not actually be using RSA rings for public-key cryptography or signatures. Instead, we’ll encode HyperLogLog values in $\Z_N^*$. Only the server who constructed $N$ and knows the factors $P$ and $Q$ can decode the HLL values. The structure of $\Z_N^*$ also allows the client to randomize what they send so that the HLL value and _only_ the HLL value is preserved by this randomization. This makes two clients with the same HLL value indistinguishable from the server’s perspective—the server learns the client’s HLL value and nothing else. Finally, this can all be done in such a way that the client only needs to generate and store a single “master key” value and can quickly and reproducibly generate resource-class-specific HLL values that are statistically independent from each other.

## Geometric RSA rings

The first thing we’ll tackle is encoding geometric samples in an RSA ring with the right frequency distribution. To do this, we define a special class of RSA rings with some additional structure. Consider the following special shape of “geometric RSA ring”:

```math
\begin{aligned}
N = P Q = (4 p + 1)(2^m q + 1)
\end{aligned}
```

Here $P$, $Q$, $p$ and $q$ are all primes and $m \ge 3$ is an integer. This ring’s structure is $\Z_N \cong \Z_P \times \Z_Q$, which cannot be broken down further as rings since $P$ and $Q$ are prime. But if we drop the additive structure and consider only the multiplicative structure, then we *can* break it down further into a product of cyclic groups:

```math
\begin{aligned}
\Z_N^* \cong (C_4 \times C_p) \times (C_{2^m} \times C_q)
\end{aligned}
```

Here $C_n$ denotes a cyclic group of order $n$, meaning that there is some generator element, $g \in C_n$, such that $g^{n-1} ≠ 1$ and $g^n = 1$ and every element of $C_n$ is of the form $g^k$ for some $k$. When $n$ is prime, every element besides $g^0 = 1$ is also a generator; when $n$ is composite, it’s more complicated. This structure connects back to the Carmichael function and the modular exponential structure of $\Z_N^*$.

Since $\Z_N^*$ is not cyclic overall — because the orders of $C_4$ and $C_{2^m}$ share a factor of four — no element can generate the entire multiplicative group. We can, however, choose a “semigenerator” element, $g \in \Z_N^*$, such that $\fmod(g,P)$ and $\fmod(g,Q)$ are generators in $\Z_P^*$ and $\Z_Q^*$, respectively. Each component of $g$ is a generator of that cyclic component, so if we take the elementwise logarithm of the cyclic components of $g$ with respect to itself then by definition we have:

```math
\begin{aligned}
\log_g(g) = (1, 1, 1, 1)
\in \Z_4 \times \Z_p \times \Z_{2^m} \times \Z_q
\end{aligned}
```

This is analogous to how $\log_x(x) = 1$ for any positive real $x$. An arbitrary unit, $x \in \Z_N^*$, can be expressed in terms of its logarithm vector with respect to $g$:

```math
\begin{aligned}
\log_g(x) = (a, b, c, d)
\in \Z_4 \times \Z_p \times \Z_{2^m} \times \Z_q
\end{aligned}
```

The behavior of $\Z_N^*$ can be entirely understood in terms of elementwise modular operations on these vectors of logarithms. Suppose we have $x_1$ and $x_2$ such that

```math
\begin{aligned}
\log_g(x_1) &= (a_1, b_1, c_1, d_1) \\
\log_g(x_2) &= (a_2, b_2, c_2, d_2) \\
\end{aligned}
```

Then we have:

```math
\begin{aligned}
\log_g(x_1 x_2)
&= (a_1, b_1, c_1, d_1) + (a_2, b_2, c_2, d_2) \\
&= (a_1 + a_2, b_1 + b_2, c_1 + c_2, d_1 + d_2) \\[0.5em]
\log_g(x_1^t)
&= t (a_1, b_1, c_1, d_1) \\
&= (t a_1, t b_1, t c_1, t d_1) \\
\end{aligned}
```

Products in $\Z_N^*$ turn into vector addition of exponents and exponentiation turns into vector scaling of exponents.

It’s worth emphasizing that the exponent in each component is modular with respect to its own modulus. This is especially significant when exponentiating by $t$, which manifests as scaling all exponents by $t$, since scaling by the same value can have very different effects on components based on their respective moduli. The Chinese Remainder Theorem guarantees that a single value $t$ can be chosen to have independent effects on coprime components, which is very powerful. Components whose moduli are not coprime, on the other hand, will always be scaled in lock step up to the greatest common divisor of their moduli. For example, if we want to scale $b \in \Z_p$ by $t_p$ and $d \in \Z_q$ by $t_q$ then we can find $t$ such that $t = t_p \bmod p$ and $t = t_q \bmod q$ and have precisely the desired effect. If we want to scale $a \in \Z_4$ and $c \in \Z_{2^m}$, on the other hand, since $4 \divides 2^m$ we cannot scale them independently: the scaling of $a$ is fully dictated by the scaling of $c$.

The key fact about our special RSA ring is that the multiplicative orders of elements of the $C_{2^m}$ component have precisely the geometric distribution we need for HyperLogLog. Here is a table of element orders, how many elements have that order, the probability of randomly choosing such an element, and the corresponding number of trailing zeros for a geometric sample value:

|   Order   |   Count   |     Probability      | Trailing zeros |
| :-------: | :-------: | :------------------: | :------------: |
|   $2^m$   | $2^{m-1}$ |    $\dfrac{1}{2}$    |      $0$       |
| $2^{m-1}$ | $2^{m-2}$ |    $\dfrac{1}{4}$    |      $1$       |
| $2^{m-2}$ | $2^{m-3}$ |    $\dfrac{1}{8}$    |      $2$       |
| $\vdots$  | $\vdots$  |       $\vdots$       |    $\vdots$    |
|    $4$    |    $2$    | $\dfrac{1}{2^{m-1}}$ |     $m-2$      |
|    $2$    |    $1$    |   $\dfrac{1}{2^m}$   |     $m-1$      |
|    $1$    |    $1$    |   $\dfrac{1}{2^m}$   |      $m$       |

This correspondence is more than just lucky coincidence. If $g$ is a generator for $C_{2^m}$, then every element of $C_{2^m}$ is of the form $g^k$ where ${} k \in \Z_{2^m}$. The order of ${} g^k {}$ in this group is determined by the largest $2^j$ that divides ${} k$—*i.e.* the number of trailing zeros when ${} k$ is written in binary:

- If ${} k$ is odd then $g^k$ is another generator with full order, $2^m$
- If $k$ has a single trailing zero then ${} g^k$ has order $2^{m-1}$
- If ${} k = 2^{m-1}$ then $g^k = -1$, the unique order $2$ value
- If ${} k = 2^m = 0 \pmod{2^m}$ then $g^k = 1$, the unique order $1$ value

In other words, sampling an element from $C_{2^m}$ and taking its order is fundamentally the same operation as sampling an $m$-bit integer and counting its trailing zeros—specifically the trailing zeros of its logarithm with respect to a generator (any generator—it doesn’t matter which one).

Suppose we take a random element of our ring, $x \in \Z_N$. If $x$ isn’t in $\Z_N^*$ then we can’t take its logarithm. It’s straightforward to test if $x$ is invertible, however: just check that $\gcd(x,N) = 1$. For large, semiprime $N$, there is an astronomically small chance that $x$ isn’t invertible. There are $PQ$ elements in $\Z_N$ and $(P-1)(Q-1) = PQ - P - Q + 1$ elements in $\Z_N^*$, so the probability of getting a non-invertible element at random is:

```math
\begin{aligned}
\frac{P + Q - 1}{PQ}
= \frac{1}{P} + \frac{1}{Q} - \frac{1}{PQ}
\end{aligned}
```

If we assume that $P \approx Q$ are large, which should be the case for RSA rings, then this is about $2/P$. When $P$ and $Q$ are large enough to be cryptographically secure, this is vanishingly small. Moreover, if someone happens to sample a non-zero $x$ that isn’t coprime with $N$, they have lucked into being able to factor $N$, which is a bigger problem. Of course, a non-malicious client would test for invertibility of $x$ and try again, kindly forgetting that they have accidentally cracked our RSA ring. But when $N$ has thousands of bits, a trillion clients could do this constantly for a trillion years and none of them would ever get lucky. In what follows, we simply assume that clients generating random elements of $\Z_N$ always get invertible ones.

Given $x \in \Z_N^*$, write its logarithm as before:

```math
\begin{aligned}
\log_g(x) = (a, b, c, d)
\end{aligned}
```

We’d like to find $k = \tz(c)$, where $\tz$ is the “trailing zeros” function, sometimes also called the 2-adic valuation and denoted as $\nu_2$, but we’ll use the more memorable function name, $\tz$, in this writeup. We can’t actually efficiently compute logarithms in $\Z_N^*$, however—even knowing the factorization, the size of $P$ and $Q$ make computing discrete logarithms in $\Z_P^*$ and $\Z_Q^*$ infeasible. Fortunately, we don’t have to explicitly compute logarithms. First, we can project onto $\Z_Q^*$ by computing $\fmod(x, Q) \in \Z_Q^*$ which has:

```math
\begin{aligned}
\log_g(\fmod(x,Q)) = (c,d) \in \Z_{2^m} \times \Z_q
\end{aligned}
```

The $C_q$ part of this can be annihilated by raising to the $q$ power: 

```math
\begin{aligned}
\log_g(\fmod(x,Q)^q) = q\,(c, d) = (qc, 0)
\in \Z_{2^m} \times \Z_q
\end{aligned}
```

 We could raise this to $q^{-1} \bmod{2^m}$ in order to undo the scaling by $q$, but we don’t actually need to: since $q$ is odd, we have $k = \tz(qc) = \tz(c)$ and $\ord(\fmod(x,Q)^q) = 2^k$. We can efficiently find $k$ by starting with $\fmod(x,Q)^q$ and squaring it in $\Z_Q^*$ until we reach one, counting how many squarings are needed.

This gives us a computationally efficient way for someone who knows the factorization of $N$ to turn a random, uniformly sampled element, $x \in \Z_N^*$, into a geometrically distributed sample value, $k = \tz(c)$. Crucially for our purposes, someone who doesn’t know the factorization of $N$ cannot compute this geometric sample value.

## Blurring fingerprints

We now have a construction where a client can sample a value with a geometric distribution without knowing the value they’ve sampled. We can’t, however, have a client just generate a random $x \in \Z_N^*$ and send that $x$ every time—this would be a massive, uniquely identifying client fingerprint. We need some way of destroying all identifying information about $x$ _except_ for the geometric sample that we want. Fortunately, this turns out to be entirely doable.

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

Being able to reach any possible values in the $C_p$ and $C_q$ components via exponentiation is predicated on $b$ and $d$ being non-zero, however. In an RSA ring with our proposed structure, if $N$ is thousands of bits (as it would be in real usage), it is astronomically unlikely to choose $x$ with zero exponents in the $C_p$ or $C_q$ components. For the sake of sheer thoroughness, let’s erase even that tiny bit of information. While we are at it, there is another bit we’d like to erase: the top bit of the $C_4$ exponent, $a$. (Why did we introduce that bit in the first place? It’s necessary for forgery prevention. More on that later.) Randomizing the high bit of $a$ forces us to randomize the high bit of $c$ as well, but that’s ok since we can choose $m$ to be as large as we want.

How do we randomize these bits while leaving the rest intact? We want to multiply by a random “noise” element that leaves the “signal” alone but randomly perturbs the rest. We start with a fully random element, $z \in \Z_N^*$, and write $\log_g(z) = (a_{z}, b_{z}, c_{z}, d_{z})$. Our noise element is then $w = \pm z^{2^{m-1}}$. To understand what this does, first look at $z^{2^{m-1}}$ alone:

```math
\begin{aligned}
\log_g(z^{2^{m-1}})
&= 2^{m-1} \, (a_{z}, b_{z}, c_{z}, d_{z}) \\
&= (0,\; 2^{m-1} b_{z},\; 2^{m-1} c_{z},\; 2^{m-1} d_{z}) \\
\end{aligned}
```

Since $m ≥ 3$, the factor $2^{m-1}$ is divisible by $4$, so the $C_4$ exponent vanishes. The $2^{m-1} c_{z} \bmod 2^m$ component is either $0$ or $2^{m-1}$ — in other words, all the lower bits are zero and the top bit is random. Since $2^{m-1}$ is coprime to both $p$ and $q$, the $C_p$ and $C_q$ exponents are fully randomized; moreover, since $p$ and $q$ are coprime to each other, they are *independently* randomized.

When the random sign is negative it adds $\log_g(-1) = (2, 0, 2^{m-1}, 0)$ to the exponent vector, flipping the top bit of both the $C_4$ and $C_{2^m}$ exponents; the high bit of $C_{2^m}$ was already random, so this doesn’t matter, but now the high bit of the $C_4$ coordinate is *also* randomized. Crucially, the high bits of these two coordinates are independently random. Overall this leaves us with a $w$ element that looks like this:

```math
\begin{aligned}
\log_g(w) = (2 a_w,\; e,\; 2^{m-1} c_w,\; e)
\end{aligned}
```
where $a_w$ and $c_w$ are independent random bits and $e \in \Z_{pq}$ also random.
where $a_w$ and $c_w$ are independent random bits and $e \in \Z_{pq}$ is also random (a single value that fixes the $C_p$ and $C_q$ exponents by the Chinese Remainder Theorem, since $p$ and $q$ are coprime).

Now that we have $w$, we multiply our client value by it after exponentiating by $t$, sending $x$ to $w x^t$:

```math
\begin{aligned}
\log_g(w x^t) = (ta + 2 a_w,\; tb + e,\; tc + 2^{m-1} c_w,\; td + e)
\end{aligned}
```

Reading this coordinate by coordinate confirms that only the geometric sample and the Jacobi parity survive:

- **$C_p$, $C_q$:** the $e$ terms make these freshly random, erasing all bulk identity; and regardless of the other values, including $t$, some $w$ produces any pair of values here.
- **$C_4$:** the parity of $ta + 2 a_w$ is the parity of $a$, since $t$ is odd—the Jacobi bit is preserved—while the high bit is scrambled by $a_w$.
- **$C_{2^m}$:** for $\tz(c) < m-1$ the sample $\tz(tc + 2^{m-1} c_w) = \tz(c)$ is untouched; at the very top, $c_w$ randomizes the top bit, merging the two rarest rungs ($\tz = m-1$ and $\tz = m$) into a single saturated level. So the readable geometric sample is $\min(\tz(c),\, m-1)$.

The parity of $a$—the Jacobi bit—is the one thing the noise deliberately preserves: $x \mapsto x^t$ scales it by the odd $t$ and $w$ shifts it by the even $2 a_w$, so it survives both. As it turns out, we need this parity bit preserved in order to prevent clients from artificially inflating their geometric samples, which is the subject of the next section.

## Fighting inflation

Clients are supposed to raise $x$ to an odd power, $t$. What happens if they raise it to an even $t$ instead? For any values $t$ and ${} c$ we have $\tz(tc) = \tz(t) + \tz(c)$. When $t$ is odd we have $\tz(t) = 0$, so $\tz(tc) = \tz(c)$. But if $t$ is even then $\tz(t) > 0$ and $x^t$ inflates the geometric sample. If a client sends $x^{2^{m-1}}$, for example, then it is guaranteed to force the sample to the saturated level $m-1$, which is supposed to be vanishingly rare, occurring with only about $1/2^{m-1}$ probability. If we require that $a = 1 \bmod 2$ in $x$, however, then the parity of $ta$ equals the parity of $t$, which lets us read the parity of $t$ from the first component of $wx^t$. Parity bit to the rescue! So we’d like to require clients to choose $x$ with odd $a$ and then check that $ta = 1$ in the $wx^t$ value that is sent. Any requests not satisfying this should be ignored for client count estimation purposes, since that indicates a malicious or malfunctioning client.

Readers may wonder how clients can check whether the $x$ that they’ve chosen has $a = 1 \bmod 2$. After all, the whole point of working in $\Z_N$ is that people who don’t know the factorization of $N$ can’t extract the components of $x \in \Z_N$. Someone who knows the factorization of $N$ can, of course, check this easily:

```math
\begin{aligned}
a \text{ even} ~~\iff~~ \fmod(x,P)^{2p} = 1
\end{aligned}
```

But the client doesn’t know $P$ and can’t check this. It is, however, possible to efficiently compute the [Jacobi symbol](https://en.wikipedia.org/wiki/Jacobi_symbol), $\Jacobi_N(x)$, without knowing the factorization of $N$. For this particular ring shape, the value of the Jacobi symbol is $(-1)^{a + c}$. In general, the Jacobi symbol is the total parity of all the log-coordinates whose moduli are even. This doesn’t tell us the value of $a$ by itself, but if $\Jacobi_N(x) = -1$ then we know that $a + c = 1 \bmod 2$. This isn’t exactly what we wanted, but it does give us what we need. Instead of requiring $a = 1$, we can require that
```math
\begin{aligned}
\Jacobi_N(wx^t) &= -1 && \iff \\
(-1)^{t(a + c)} &= -1 && \iff \\
t(a + c) &= 1 \bmod 2 && \iff \\
t = a + c &= 1 \bmod 2
\end{aligned}
```
In words, this requirement forces both $t$ and $a + c$ to be odd. We should also verify that conditioning on $\Jacobi_N(x) = -1$ does not change the distribution of $\tz(c)$. Since $a \in \{0,1\}$ is uniform and independent of $c$, for every value of $\tz(c)$ the probability mass is split equally between the two Jacobi values:

| $x \in$ | $\Z_N^*$ | $J_N^+$ | $J_N^-$ |
|----:|:---:|:---:|:---:|
| $\tz(c) = 0$ | $\tfrac{1}{2}$ | $\tfrac{1}{4}$  | $\tfrac{1}{4}$  |
| $\tz(c) = 1$ | $\tfrac{1}{4}$ | $\tfrac{1}{8}$  | $\tfrac{1}{8}$  |
| $\tz(c) = 2$ | $\tfrac{1}{8}$ | $\tfrac{1}{16}$ | $\tfrac{1}{16}$ |
| $\vdots$     | $\vdots$       | $\vdots$        | $\vdots$        |
| total        | $1$            | $\tfrac{1}{2}$  | $\tfrac{1}{2}$  |

This table uses $J_N^± = \set{x \in \Z_N^* \st \Jacobi_N(x) = ± 1}$.
Every row is split evenly between $J_N^+$ and $J_N^-$, so conditioning on $\Jacobi_N(x) = -1$ leaves the distribution of $\tz(c)$ intact. Equivalently, $\tz(c)$ and $\Jacobi_N(x)$ are independent random variables.

In summary, instead of requiring $x$ with $a = 1 \bmod 2$, which would work but cannot be checked by clients, we require $x$ with $\Jacobi_N(x) = -1$, which they can check. This implies that $\Jacobi_N(wx^t) = -1$ if and only if $t$ is odd. The server can check that this is the case and disregard any requests that don’t satisfy this parity requirement.

## Bucket brigade

We now have the ability for clients to blindly sample geometric values, but for full HyperLogLog sampling we also need a uniformly distributed bucket value. Can we modify our RSA ring to include a bucket value too? We certainly can:

```math
\begin{aligned}
N = P Q = (4 B p + 1)(2^m q + 1)
\end{aligned}
```

This is the full RSA ring shape for encrypted HyperLogLog sampling. As before, $P$, $Q$, $p$, and $q$ are distinct, odd primes, $m$ is the geometric sample cap (samples run from $0$ to $m-1$), and $B$ is the number of HyperLogLog buckets, which must be odd and coprime to everything else. The multiplicative structure of a HyperLogLog RSA ring is:

```math
\begin{aligned}
\Z_N^* \cong (C_4 \times C_B \times C_p) \times (C_{2^m} \times C_q)
\end{aligned}
```

The $C_B$ component encodes the encrypted bucket value. The multiplicative order of the $C_{2^m}$ component encodes the geometric sample value. The low bit of the $C_4$ component is used to ensure that clients don’t inflate geometric samples by raising values to even powers. The $C_p$ and $C_q$ components are what make values hard to decode—we don’t actually care about the values in these components. Since this is starting to be a lot of components and we don’t actually care about $C_p$ and $C_q$, we’ll combine them into a single cyclic component with order $pq$:

```math
\begin{aligned}
\Z_N^* \cong C_4 \times C_B \times C_{2^m} \times C_{pq}
\end{aligned}
```

Note that $C_a \times C_b \cong C_{ab}$ only holds if $a$ and $b$ are coprime, which is the case here since $p$ and $q$ are distinct primes.

The client randomly chooses and saves a persistent $x \in \Z_N^*$ value. For each new request, it randomly chooses $w \in \pm(\Z_N^*)^{B2^{m-1}}$ and $t = 1 \bmod{2B}$ and sends $y = wx^t$. It’s not too hard to intuitively see that the only meaningful pieces of information conveyed about $x$ are:

1. The parity bit
2. The bucket value
3. The geometric sample

We can analyze it in terms of generator exponents. Choose a semigenerator, $g$, in this new ring and write $\log_g(x) = (a, b, c, d)$. As in the geometric-only ring, the noise has $\log_g(w) = (2 a_w, 0, 2^{m-1} c_w, e)$ with $a_w$ and $c_w$ random bits—the extra factor of $B$ in the exponent only zeroes the $C_B$ component, leaving the bucket untouched. Then:

```math
\begin{aligned}
\log_g(wx^t)
&= (ta + 2a_w,\; b,\; tc + 2^{m-1} c_w,\; td + e) \\
\end{aligned}
```

Spelling these four components out:

- Parity bit: $ta + 2a_w = a \bmod 2$ since $t = 1 \bmod 2$
- Bucket index: $tb = b$ since $t = 1 \bmod B$
- Geometric sample: $\min(\tz(tc + 2^{m-1} c_w), m-1) = \min(\tz(tc), m-1) = \min(\tz(c), m-1)$
- The rest: $td + e \bmod{pq}$ is freshly random in each request

The core HyperLogLog Over RSA construction is now in place. The server constructs a ring with the shape $N = PQ = (4Bp+1)(2^m q+1)$, keeping $P$, $Q$, $p$, $q$ secret and publishing $N$ together with a semigenerator $g$. A client generates a persistent secret $x \in \Z_N^*$ with $\Jacobi_N(x) = -1$ and, for each request, sends a freshly randomized token $y = wx^t$ where $t \equiv 1 \bmod{2B}$ and $w$ is chosen from $\pm(\Z_N^*)^{B2^{m-1}}$. The server discards any request where $\Jacobi_N(y) \neq -1$, since that indicates a misbehaving client. For the remaining requests, it uses the secret factors, $P$ and $Q$, to extract the bucket index, $b$, from the $C_B$ component and the geometric sample, $k = \min(\tz(c),\, m-1)$, from the $C_{2^m}$ component, yielding a decrypted HLL sample pair, $(b, k)$. Clients cannot steer or bias samples without knowing the factorization of $N$. The next section completes the protocol by making it work nicely with resource class sharding.

## Master keys

The construction so far assumes that each client chooses a single random $x \in J_N^-$ that encodes a single HLL sample. But one more ingredient is needed for the full protocol. [Resource class sharding](03-counting-users.md#Resource-class-sharding) requires that each client has separate, statistically independent secrets for each one of an unbounded set of resource classes. If a client reused the same $x$ in every resource class, it would have the same HLL value in every class. Clients must use a different, unlinkable $x$ value for each resource class.

### Possible approaches

The naive solution is to store a separately generated $x$ per class. This works but involves a distasteful amount of storage: 1KB per $\Z_N^*$ value and perhaps 5k resource classes over a user’s lifetime, so nearly a megabyte of persistent storage (per server that the client talks to, which, for most clients, is only one). Another paranoid reason not to store a secret per class: such a cache provides someone who gets access to the client’s storage a complete record of every package they’ve ever installed. Of course, the attacker can also just look at what packages they have installed, but a cache of resource-class secrets would persist even when all traces of an installed package have been removed. If we provide a way to clear or purge the cache, we are essentially refreshing the client’s identity for counting purposes, which is not catastrophic, but is undesirable.

What we really want is a way for each client to generate a single “master key” and derive individual resource-class-specific $x$ values from the master key and resource class. The most obvious way to do this is to hash the master key and the resource class together into $\Z_N$:

```math
\begin{aligned}
x = \hash_N(\text{key}, \text{class}) \in \Z_N^*
\end{aligned}
```

Here $\hash_N$ is some cryptographic hash function taking a tuple of arguments to a uniformly distributed value in $\Z_N$. The complication is that we need to guarantee that $\Jacobi_N(x) = -1$ for inflation prevention. We could include a counter in the hash:

```math
\begin{aligned}
x_i = \hash_N(\text{key}, \text{class}, i)
\end{aligned}
```

and use the first $x_i$ such that $\Jacobi_N(x_i) = -1$. Since about half of the values in $\Z_N$ have a Jacobi symbol of $-1$, it should only take a few attempts to find a usable $x_i$. This works, but is somewhat inelegant and involves computing a variable number of (large) hash values and Jacobi symbols for each counter value. Keep in mind that this work has to be done for *each request* the Pkg client wants to make to a server, so having it be at all expensive or non-deterministic would be best to avoid.

### What we actually do

Can we design a scheme where the client generates a single pre-validated master key and constructs $x$ for each resource class so that it always has a negative Jacobi symbol? And ideally, every possible value in the negative Jacobi set would be reachable for some potential resource class. As it happens, this is possible. The key insight is that our particular shape of RSA ring, while not cyclic, is very close to it. A semigenerator $g$ projects onto every cyclic factor, so its powers already sweep the entire $C_{pq}$ bulk along with all of $C_B$ and $C_{2^m}$. The only issue is that the $C_4$ coefficient and the low two bits of the $C_{2^m}$ component are scanned in lock step. The subgroup generated by $g$ is index two in $J_N^+$ and index four in $\Z_N^*$ overall. If we combine $g$ with $-1$ we generate all of $J_N^+$, and if we combine with any element in $J_N^-$, we get all of $\Z_N^*$. Recall that $\Z_N^*$ has this multiplicative structure:

```math
\begin{aligned}
\Z_N^* \cong
C_4 \times C_B \times C_{2^m} \times C_{pq}
\end{aligned}
```

All semigenerators in this ring have $\Jacobi_N(g) = 1$ and multiplicative order $\ord(g) = \lambda(N) = B 2^m pq$. So $\gen{g}$ sits inside the positive-Jacobi subgroup $J_N^+$ but reaches only *half* of it: no power of $g$ hits the order-2 element of the $C_4$ factor, and that element is exactly $-1$. Fixing any $x_0$ with $\Jacobi_N(x_0) = -1$, the four cosets of $\gen{g}$ split $\Z_N^*$ as

```math
\begin{aligned}
J_N^+ = \gen{g} \sqcup -\gen{g}, \qquad J_N^- = x_0\gen{g} \sqcup -x_0\gen{g}
\end{aligned}
```

Ranging $k$ over all exponents, $x_0 g^k$ therefore covers the coset $x_0\gen{g}$ — half of $J_N^-$. That is enough. The missing half is $-x_0\gen{g}$, and $-1$ lives in the white-noise subgroup $W$ that the client mixes into every token, so the two halves are identified once the token is randomized. In other words, *modulo the noise*, every $x \in J_N^-$ is reachable as $x_0 g^k$ — which is all the anonymity guarantee needs. The client can easily pick a valid $x_0$ since all they have to check is that $\Jacobi_N(x_0) = -1$. The client cannot, on the other hand, check if $g$ is a semigenerator since it doesn’t know the factorization of $N$, but the server can do this and publish a common $g$ value along with $N$. Clients cannot check that $g$ is actually a semigenerator, but there’s no real harm done if it isn’t (we address this rigorously in our security analysis).

Putting it together:

- The server, when generating the ring, also chooses and publishes a common “semigenerator” element, $g \in \Z_N^*$;
- The client, when downloading the ring parameters for the first time, also chooses and saves a random $x_0 \in \Z_N^*$ with $\Jacobi_N(x_0) = -1$. This $x_0$ is the client’s master key.

Since half of the values in $\Z_N$ have negative Jacobi symbol, a viable $x_0$ is quick to find, and it only has to be done once for a new ring. Regardless of which $x_0$ the client chooses, every $x \in J_N^-$ is reachable as $x_0 g^k$ modulo the noise subgroup, as just discussed. The client’s choice of $x_0$ changes how exponents map to $x$ values in a way that we’ll explore below. Write the logarithm vector of the master key as:

```math
\begin{aligned}
\log(x_0) = (a, b, c, d)
\end{aligned}
```

Since $\Jacobi_N(x_0) = (-1)^{a + c} = -1$ we know that $a + c = 1 \bmod 2$. In other words, exactly one of $a$ or $c$ is odd. Other than this relation, the four values are completely random. For each resource class, the client derives a client- and class-specific hash value, $h$:

```math
\begin{aligned}
h &= \hash(x_0, \text{class})
\end{aligned}
```

This uses $x_0$ as salt string, rather than as a ring element; we will *also* use $x_0$ as a ring element later. Using $x_0$ as salt makes the hash value, $h$, client-specific in a highly unpredictable manner. The server cannot learn a client’s $h$ value, but even if they did, it doesn’t reveal anything about $h$ values in other classes. The hash function, $\hash$, used here, is different from the previously proposed $\hash_N$ — this hash’s output is used as an exponent, rather than a ring element. This seems like a minor difference, but we’ll see that this usage requires far fewer output bits and makes our derivation significantly faster. Next, we use $x_0$ and $h$ to derive the client’s personal resource class secret:

```math
\begin{aligned}
x_h &= x_0 g^h
\end{aligned}
```

This is the second usage of $x_0$, this time as a ring element. Note that the construction of $x_h$ guarantees that it has negative Jacobi symbol:

```math
\begin{aligned}
\Jacobi_N(x_h) = \Jacobi_N(x_0) \Jacobi_N(g)^h = -1
\end{aligned}
```

This $x_h$ plays the role of $x$ in previous sections, where we presumed it to be a random value with negative Jacobi symbol. As long as the hash values cover a sufficiently large range, $x_h$ is an arbitrary value in the coset $x_0\gen{g}$, which is half of $J_N^-$:

```math
\ord(g)
= \lambda(N)
= 2^m B p q
= \tfrac{1}{2}\norm{J_N^+}
= \tfrac{1}{4}\norm{\Z_N^*}
```

If $h$ can take on every value in $\Z_{\lambda(N)}$ then $x_h$ takes on every value in $x_0\gen{g}$ and white noise in $w$ supplies the other half of $J_N^-$. Naively, even this presents a problem: $\lambda(N)$ is large and necessarily unknown to the client. Fortunately, $h$ does not actually need to cover this entire range because the client doesn’t send $x_h$ as is: it actually sends $y = wx_h^t$ where $w \in W$ is random “white noise” and $t$ is a random exponent with $t = 1 \bmod 2B$. If $\log_g(x_0) = (a, b, c, d)$ and $\log_g(w) = (2 a_w, 0, 2^{m-1} c_w, e)$ with $a_w$ and $c_w$ random bits, then we have:

```math
\begin{aligned}
\log(y)
&= \log(w x_h^t)
= \log(w (x_0 g^h)^t) \\
&= t \, (a + h, b + h, c + h, d + h) + (2 a_w, 0, 2^{m-1} c_w, e) \\
&= (t(a + h) + 2 a_w, \; b + h, \; t(c + h) + 2^{m-1} c_w, \; t(d + h) + e) \\
\end{aligned}
```

The components of $y$ that convey information are:

- The parity bit: $t(a + h) + 2 a_w = a + h \bmod 2$
- The bucket index: $b + h \in \Z_B$
- The geometric sample: $\min(\tz(t(c + h) + 2^{m-1} c_w),\, m-1) = \min(\tz(c + h), m-1) \in \set{0, \dots, m-1}$

Here we can see that the real requirement on $h$ is that $(b + h, c + h)$ covers all of $\Z_B \times \Z_{2^m}$ fairly uniformly as $h$ takes on different values. To ensure this it’s sufficient to ensure that $h$ is sampled from a modulus, $M$ such that $\fmod(M, B2^m)$ is tiny relative to $M \div B2^m$. We could use an exact multiple of $B2^m$, which makes the modulus zero. But that’s inconvenient, since real world hashes have power-of-two outputs. Fortunately, if $M$ is sufficiently large, the bias (ratio) is negligible. For any modern hash function like SHA-256, this is the case. In our reference implementation, we actually use the output of HMAC-SHA-256 (keyed by the SHA-256 hash of $x_0$), truncated into a 128-bit integer value, which is still more than large enough to guarantee effectively uniform coverage.

What information do we make sure is **not** conveyed?

- The leading digits of $t(c+h) + 2^{m-1} c_w$ should be fully random
- The last component, $t(d + h) + e$, should be fully random

To ensure the former, we want $t$ to be able to take every odd residue class in $\Z_{2^m}$ which is accomplished by choosing $i \in [0, 2^{m-1})$ and letting $t = 2Bi + 1$. To ensure the latter, we want $w$ to take every possible value in $W$ which is accomplished by choosing $z \in \Z_N^*$ at random, picking a random sign, and letting $w = \pm z^{B2^{m-1}}$ (we don’t actually care about $t$ for this factor).

Putting it all together, for each request a client makes, this scheme requires the client to:

- Compute $h = \hash(x_0, \text{class})$ — one HMAC-SHA-256 operation
- Compute $g^h \bmod N$ — one modular exponentiation
- Compute $x_h = x_0 g^h \bmod N$ — one modular multiplication
- Generate random $z \in \Z_N^*$
- Compute $w = \pm z^{B2^{m-1}}$ — one modular exponentiation (and a negation)
- Generate random $i \in \Z_{2^{m-1}}$
- Compute $t = 2Bi + 1$ — a few small arithmetic operations
- Compute $x_h^t \bmod N$ — one modular exponentiation
- Compute $y = w x_h^t$ — one modular multiplication

When generating $z \in \Z_N^*$ the client technically needs to choose $z \in \set{1, \dots, N-1}$ and then also check that $\gcd(z, N) = 1$. Otherwise there is an astronomically slim chance that

```math
\Jacobi_N(y) = \Jacobi_N(w) = \Jacobi_N(z) = 0.
```

This renders the request invalid, so the server would disregard it since it fails the $\Jacobi_N(y) = -1$ check. This is so unlikely, however, that we will just take that chance on each request rather than spending time checking invertibility. As a practical matter, this will never happen—in the literal sense that for realistic values of $N$, it is so unlikely that even over centuries of trillions of clients using this protocol, chances are that none of them will ever choose $z$ that isn’t in $\Z_N^*$. If it does happen, the consequences are virtually nil: a single request will be ignored when computing a client count estimate that is approximate anyway. In other words, invertibility is absolutely not worth checking for here.

Here’s the actual code for this in our test implementation:
```julia
h = hash_resource_class(x₀, class)
x = modmul(x₀, powermod(g, h, N), N)
z = rand(rng, 1:N-1)
w = powermod(z, oftype(N, B) << (m-1), N)
if rand(rng, Bool)
    w = N - w
end
i = rand(rng, zero(N):(oftype(N, 1) << (m-1)) - one(N))
t = 2 * oftype(N, B) * i + one(N)
y = modmul(w, powermod(x, t, N), N)
```
As you can see, this is a very literal rendition of the steps described. On my M2 MacBook Air, for $1024$-bit $N$, this takes around $80$ microseconds; for $2048$-bit $N$, it takes $435$ microseconds, less than half a millisecond. In other words, these computations are easily fast enough to perform during each Pkg request — the network request itself will take orders of magnitude more time.
