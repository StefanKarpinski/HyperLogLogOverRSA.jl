# Encrypted Geometric Sampling

What alternative is there to signing HyperLogLog values? In the last section we encrypted HLL values. What if clients could sample encrypted HLL values for themselves? Of course, the way we discussed encrypting them previously, each distinct HLL value would have an equal chance of being chosen, which is entirely the wrong distribution—the whole point is that it's geometric. But can we encrypt HLL values so that sampling encrypted values has the right distribution? This requires more common HLL values to appear in distinct encrypted forms many times: the most common values would have to appear $2^{m-1}$ times as often as the least common ones, where $m$ is the maximum geometric sample value. This is certainly possible—we can just encrypt a 128-bit unsigned random value and derive the HLL value from it on the server side. But that's a massive client fingerprint. We need another crucial ingredient: the client needs to be able to randomize what they send so that any client with the same HLL value could have sent that value. This seemingly impossible trick turns out to be possible in the mathematical setting of RSA rings.

## RSA Ring Refresher

An RSA ring is a modular integer ring, $\Z_N$, where $N$ is a large, balanced semiprime—*i.e.* a product of two similarly sized primes: $N = P Q$. When $N$ is large enough, this is a hard factoring problem—hard enough that it is currently impossible when $\log_2(N) ≥ 1024$. For good measure, people often use $2048$-bit and $4096$-bit semiprimes these days. Constructing such an $N$, on the other hand, is straightforward and takes no more than a fraction of a second: find two distinct primes, $P$ and $Q$, of appropriate size, and multiply them to get $N = PQ$. The party that generates $N$ can safely publish it and keep $P$ and $Q$ secret. This gives them knowledge of the structure of $\Z_N$ that no one else has, yet everyone can perform computations in $\Z_N$. This was one of the first public-key encryption schemes and is still widely used. If all you want to do is public key cryptography, the additional structure of $\Z_N$ is actually a liability and has to be carefully tiptoed around. But for us the additional structure is crucial to making what we want to do possible.

The classic use of an RSA ring is for public-key cryptography. The party that generates $N$ chooses an exponent, $E$, and computes $D = E^{-1} \bmod{\lambda(N)}$. Here $\lambda(N)$ is the [Carmichael function](https://en.wikipedia.org/wiki/Carmichael_function), defined as the smallest exponent such that for all $x \in \Z_N$, $x^{\lambda(N)} = 1 \bmod N$. This function is crucial for understanding exponents in $\Z_N^*$ because exponents behave modularly as well, but with a modulus of $\lambda(N)$: *i.e.* $x^a = x^b \bmod N$ whenever $a = b \bmod{\lambda(N)}$. For semiprime $N = PQ$, we have $\lambda(N) = \lcm(P-1, Q-1)$, which is easy for someone who knows $P$ and $Q$ to compute, but infeasible for someone who only knows $N$. Since $DE = 1 \bmod{\lambda(N)}$ we have $x^{DE} = x^1 = x \bmod N$ for all $x \in \Z_N^*$. Thus, the exponents $E$ and $D$ form a public/private key pair:

- You can publish $E$ and anyone can encrypt a message $x \in \Z_N^*$ by exponentiating by $E$: $c = x^E \bmod N$.
- Only someone who knows $D$ can decrypt the message by exponentiating by $D$: $c^D = (x^E)^D = x^{DE} = x \bmod N$.

This system can also be used for signatures. To sign a message, $x$, the party who knows $D$ computes $s = x^D \bmod N$ and sends $(x, s)$ as a message/signature pair. Anyone with the public key can check the signature like this: 

```math
\begin{aligned}
s^E = (x^D)^E = x^{DE} = x \bmod N.
\end{aligned}
```

 If $s^E$ doesn't match the message, the signature is invalid. Anyone can check a signature, but only someone who knows $E$ can generate a valid one, so a valid signed message proves that the message is authentically from the party who knows $E$.

We will not actually be using RSA rings for public-key cryptography or signatures. Instead, we'll encode HyperLogLog values in $\Z_N^*$. Only the server who constructed $N$ and knows the factors $P$ and $Q$ can decode the HLL values. The structure of $\Z_N^*$ also allows the client to randomize what they send so that the HLL value and _only_ the HLL value is preserved by this randomization. This makes two clients with the same HLL value indistinguishable from the server's perspective—the server learns the client's HLL value and nothing else. Finally, this can all be done in such a way that the client only needs to generate and store a single "master key" value and can quickly and reproducibly generate resource-class-specific HLL values that are statistically independent from each other.

## Geometric RSA Rings

The first thing we'll tackle is encoding geometric samples in an RSA ring with the right frequency distribution. To do this, we define a special class of RSA rings with some additional structure. Consider the following special shape of "geometric RSA ring":

```math
\begin{aligned}
N = P Q = (2 p + 1)(2^m q + 1)
\end{aligned}
```

Here $P$, $Q$, $p$ and $q$ are all primes and $m$ is some positive integer. This ring's structure is $\Z_N \cong \Z_P \times \Z_Q$, which cannot be broken down further as rings since $P$ and $Q$ are prime. But if we drop the additive structure and consider only the multiplicative structure, then we *can* break it down further into a product of cyclic groups:

```math
\begin{aligned}
\Z_N^* \cong (C_2 \times C_p) \times (C_{2^m} \times C_q)
\end{aligned}
```

Here $C_n$ denotes a cyclic group of order $n$, meaning that there is some generator element, $g \in C_n$, such that $g^{n-1} ≠ 1$ and $g^n = 1$ and every element of $C_n$ is of the form $g^k$ for some $k$. When $n$ is prime, every element besides $g^0 = 1$ is also a generator; when $n$ is composite, it's more complicated. This structure connects back to the Carmichael function and the modular exponential structure of $\Z_N^*$. We choose a "semigenerator" element, $g \in \Z_N^*$, such that $\fmod(g,P)$ and $\fmod(g,Q)$ are generators in $\Z_P^*$ and $\Z_Q^*$, respectively. Each component of $g$ is a generator of that cyclic component, so if we take the elementwise logarithm of the cyclic components of $g$ with respect to itself then by definition we have:

```math
\begin{aligned}
\log_g(g) = (1, 1, 1, 1)
\in \Z_2 \times \Z_p \times \Z_{2^m} \times \Z_q
\end{aligned}
```

This is analagous to how $\log_x(x) = 1$ for any positive real $x$. An arbitrary unit, $x \in \Z_N^*$, can be expressed in terms of its logarithm vector with respect to $g$:

```math
\begin{aligned}
\log_g(x) = (a, b, c, d)
\in \Z_2 \times \Z_p \times \Z_{2^m} \times \Z_q
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

It's worth emphasizing that the exponent in each component is modular with respect to its own modulus. This is especially significant when exponentiating by $t$, which manifests as scaling all exponents by $t$, since scaling by the same value can have very different effects on components based on their respective moduli. The Chinese Remainder Theorem guarantees that a single value $t$ can be chosen to have independent effects on coprime components, which is very powerful. Components whose moduli are not coprime, on the other hand, will always be scaled in lock step up to the greatest common denominator of their moduli. For example, if we want to scale $b \in \Z_p$ by $t_p$ and $d \in \Z_q$ by $t_q$ then we can find $t$ such that $t = t_p \bmod p$ and $t = t_q \bmod q$ and have precisely the desired effect. If we want to scale $a \in \Z_2$ by $t_2$ and $c \in \Z_{2^m}$ by $t_{2^m}$, on the other hand, we are forced to have $t = t_2 = t_{2^m} \bmod 2$, so they cannot be chosen independently—$t_2$ and $t_{2^m}$ must have the same parity.

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

In other words, sampling an element from $C_{2^m}$ and taking its order is fundamentally the same operation as sampling an $m$-bit integer and counting its trailing zeros—specifically the trailing zeros of its logarithm with respect to a generator (any generator—it doesn't matter which one).

Suppose we take a random element of our ring, $x \in \Z_N$. If $x$ isn't in $\Z_N^*$ then we can't take its logarithm. It's straightforward to test if $x$ is invertible, however: just check that $\gcd(x,N) = 1$. For large, semiprime $N$, there is an astronomically small chance that $x$ isn't invertible. There are $PQ$ elements in $\Z_N$ and $(P-1)(Q-1) = PQ - P - Q + 1$ elements in $\Z_N^*$, so the probability of getting a non-invertible element at random is:

```math
\begin{aligned}
\frac{P + Q - 1}{PQ}
= \frac{1}{P} + \frac{1}{Q} - \frac{1}{PQ}
\end{aligned}
```

If we assume that $P \approx Q$ are large, which should be the case for RSA rings, then this is about $2/P$. When $P$ and $Q$ are large enough to be cryptographically secure, this is vanishingly small. Moreover, if someone happens to sample a non-zero $x$ that isn't coprime with $N$, they have lucked into being able to factor $N$, which is a bigger problem. Of course, a non-malicious client would test for invertibility of $x$ and try again, kindly forgetting that they have accidentally cracked our RSA ring. But when $N$ has thousands of bits, a trillion clients could do this constantly for a trillion years and none of them would ever get lucky. In what follows, we simply assume that clients generating random elements of $\Z_N$ always get invertible ones.

Given $x \in \Z_N^*$, write its logarithm as before:

```math
\begin{aligned}
\log_g(x) = (a, b, c, d)
\end{aligned}
```

We'd like to find $k = \tz(c)$, where $\tz$ is the "trailing zeros" function, sometimes also called the 2-adic valuation and denoted as $\nu_2$, but we'll use the more memorable function name, $\tz$, in this writeup. We can't actually efficiently compute logarithms in $\Z_N^*$, however—even knowing the factorization, the size of $P$ and $Q$ make computing discrete logarithms in $\Z_P^*$ and $\Z_Q^*$ infeasible. Fortunately, we don't have to explicitly compute logarithms. First, we can project onto $\Z_Q^*$ by computing $\fmod(x, Q) \in \Z_Q^*$ which has:

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

 We could raise this to $q^{-1} \bmod{2^m}$ in order to undo the scaling by $q$, but we don't actually need to: since $q$ is odd, we have $k = \tz(qc) = \tz(c)$ and $\ord(\fmod(x,Q)^q) = 2^k$. We can efficiently find $k$ by starting with $\fmod(x,Q)^q$ and squaring $\Z_Q^*$ until we reach one and counting how many squarings are needed.

This gives us a computationally efficient way for someone who knows the factorization of $N$ to turn a random, uniformly sampled element, $x \in \Z_N^*$, into a geometrically distributed sample value, $k = \tz(c)$. Crucially for our purposes, someone who doesn't know the factorization of $N$ cannot compute this geometric sample value.
