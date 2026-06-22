# Master Keys

At this point we have a clever scheme that allows clients to randomly choose $x \in \Z_N^*$ and have it encode a correctly distributed HyperLogLog value such that:

1. The client cannot decode or bias which HyperLogLog value they have sampled.
2. The client's identity is perfectly obscured by sending $w x^t$ instead of $x$, where $w$ and $t$ are random values such that ${} w \in W = (\Z_N^*)^{B2^m} {}$ and $t = 1 \bmod{2B}$.

This is already quite good. Recall, however that we want to use different, independent values of $x$ for different resource classes. Clients could simply generate and store different independently random $x$ values for each resource class. But as we discussed when considering signed HLLs, this entails a nontrivial amount of storage that looks quite bad: assuming 256 bits per resource class key and 1024 bits per $\Z_N^*$ value, that's 160 bytes per resource class; if the typical user installs 5k packages over time, that's 0.8MB of storage. Once again, this is feasible on modern hardware, but it would be great if each client could generate and store just a single reasonably sized secret.

In essence, we want a way for each client to generate a random "master key" and derive individual resource-class-specific $x$ values from the master key and resource class. The most obvious way to do this is to hash the master key and the resource class together into $\Z_N$:

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

and use the first $x_i$ such that $\Jacobi_N(x_i) = -1$. Since about half of the values in $\Z_N$ have Jacobi symbol of $-1$, it shouldn't take too many attempts to find a usable $x_i$. This works, but is somewhat inelegant and involves computing a variable number of (large) hash values and Jacobi symbols for each. Keep in mind that this work has to be done for each request the Pkg client wants to make to a server. Can we come up with a scheme where the client generates a single pre-validated master key and constructs $x$ for each resource class so that it always has a negative Jacobi symbol? And ideally, every possible value in the negative Jacobi set would be reachable for some potential resource class. As it happens, this is possible.

The key insight is that our particular shape of RSA ring, while not cyclic, is very nearly cyclic. Recall that it has this multiplicative structure:

```math
\begin{aligned}
\Z_N^* \cong
C_2 \times C_B \times C_{2^m} \times C_{pq}
\end{aligned}
```

The only shared factor among cyclic component orders is a single factor of two shared by $C_2$ and $C_{2^m}$. This means that if $g$ is a semigenerator, then every element in $\Z_N^*$ is either of the form $g^k$ or $x_0 g^k$ where $x_0$ is any fixed element with $\Jacobi_N(x_0) = -1$. Recall that a semigenerator, $g \in \Z_N^*$, has $\fmod(x,P)$ and $\fmod(x,Q)$ both generators. All semigenerators in this ring have $\Jacobi_N(g) = 1$, so every element of the form $g^k$ has positive Jacobi symbol while every element of the form $x_0 g^k$ has negative Jacobi symbol. This gives us a natural way to generate every negative Jacobi value: fix $x_0$ and $g$ and let $k$ range over exponent values: for every $x \in J_N^-$ there is some $k$ such that $x = x_0 g^k$. The client can easily pick a valid $x_0$ since all they have to check is that $\Jacobi_N(x_0) = -1$. On the other hand, the client cannot check if $g$ is a semigenerator since it doesn't know the factorization of $N$. The server can do this, however, and publish $g$ along with $N$.

This, then, is the generation part of our master key scheme:

- The server, when generating the ring, also chooses and publishes a common "semigenerator" element, $g \in \Z_N^*$;
- The client, when downloading the ring parameters for the first time, also chooses and saves a random $x_0 \in \Z_N^*$ with $\Jacobi_N(x_0) = -1$.

Since half of the values in $\Z_N$ have negative Jacobi symbol, a viable $x_0$ is quick to find, and it only has to be done when once for a new ring. Regardless of which $x_0$ the client chooses, every $x \in \Z_N^*$ with $\Jacobi_N(x) = -1$ has $x = x_0 g^k$ for some $k$. The client's choice of $x_0$ changes how exponents map to $x$ values in a way that we'll explore below — this is the client's master key. We'll write its logarithm vector as:

```math
\begin{aligned}
\log(x_0) = (a, b, c, d)
\end{aligned}
```

Since $\Jacobi_N(x_0) = (-1)^{a + c} = -1$ we know that $a + c = 1 \bmod 2$. In other words, exactly one of $a$ or $c$ is odd. Other than this relation, the values are completely random.

For each resource class, the client derives a class-specific value by multiplying their master key with a hashed power of the semigenerator:

```math
\begin{aligned}
x_h = x_0 g^h = x_0 g^{\hash(x_0,\,\text{class})}
\end{aligned}
```

This hash function, $\hash$, is different than the previously proposed $\hash_N$ — this hash value is used as an exponent, rather than a ring element. Later, we'll see that this requires far fewer output bits. In addition to using $x_0$ to construct values with negative Jacobi symbol, the client also uses $x_0$ as salt when hashing the resource class. This makes the hash value, $h$, client-specific in a highly unpredictable manner. Anything a server might learn about a client's $h$ value in one resource class—which should be nothing, but let's be cautious—doesn't reveal anything about the same client's $h$ values for other classes since the server cannot compute $x_0$ from $h$ since $\hash$ is a secure one-way function.

The construction of $x_h$ guarantees that it has negative Jacobi symbol:

```math
\begin{aligned}
\Jacobi_N(x_h) = \Jacobi_N(x_0) \Jacobi_N(g)^h = -1
\end{aligned}
```

This $x_h$ plays the role of $x$ in previous sections, where we presumed it to be a random value with negative Jacobi symbol. As long as the hash values cover a sufficiently large range, $x_h$ is an arbitrary value in $J_N^-$:

```math
\ord(g)
= \norm{J_N^+}
= \tfrac{1}{2}\norm{\Z_N^*}
%= \tfrac{1}{2}\varphi(N)
= 2^m B p q
= \lambda(N)
```

If $h$ can take on every value in $\Z_{\lambda(N)}$ then $x_h$ takes on every possible value in $J_N^-$. Naively, this presents a problem: $\lambda(N)$ is large and unknown to the client. Fortunately, $h$ does not actually need to cover this entire range because the client doesn't send $x_h$ as is: it actually sends $y = wx_h^t$ where $w \in W$ is random "white noise" and $t$ is a random exponent with $t = 1 \bmod 2B$. If $\log_g(x_0) = (a, b, c, d)$ and $\log_g(w) = (0, 0, 0, e)$ then we have:

```math
\begin{aligned}
\log(y)
&= \log(w x_h^t)
= \log(w (x_0 g^h)^t) \\
&= t \, (a + h, b + h, c + h, d + h) + (0, 0, 0, e) \\
&= (a + h, b + h, t(c + h), t(d + h) + e) \\
\end{aligned}
```

The components of $y$ that convey information are:

- The parity bit: $a + h \in \Z_2$
- The bucket index: $b + h \in \Z_B$
- The geometric sample: $\tz(t(c + h)) \in \set{0, \dots, m}$

What information do we want to make sure is _not_ conveyed?

- The leading digits of $t(c+h)$ should be fully random
- The last component, $t(d + h) + e$, should be fully random

To ensure the former, we want $t$ to be able to take every odd residue class in $\Z_{2^m}$ which is accomplished by choosing $i \in [0, 2^{m-1})$ and letting $t = 2Bi + 1$. To ensure the latter, we want $w$ to take every possible value in $W$ which is accomplished by choosing $z \in \Z_N^*$ at random and letting $w = z^{B2^m}$ (we don't actually care about $t$ for this). What's really required is that $(b + h, c + h)$ covers all of $\Z_B \times \Z_{2^m}$ as $h$ takes on different values. All we need is for $\hash$ to be a hash function with at least $\log_2(B) + m$ bits of output. In practice, we have $\log_2(B) \approx 12$ and $m = 63$ so we need more than $75$ bits of hash output. In our test implementation, we use the first $128$ bits of the SHA256 hash of $x_0$ with the resource class.

Putting it all together, for each request a client makes, this scheme requires the client to:

- Compute $h = \hash(x_0, \text{class})$ — one SHA256 hash operation
- Compute $g^h \bmod N$ — one modular exponentiation
- Compute $x_h = x_0 g^h \bmod N$ — one modular multiplication
- Generate random $z \in \Z_N^*$
- Compute $w = z^{B2^m}$ — one modular exponentiation
- Generate random $i \in \Z_{2^{m-1}}$
- Compute $t = 2Bi + 1$ a few small arithmetic operations
- Compute $x_h^t \bmod N$ — one modular exponentiation
- Compute $y = w x_h^t$ — one modular multiplication

When generatting $z \in \Z_N^*$ the client technically needs to choose $z \in \set{1, \dots, N-1}$ and then also check that $\gcd(z, N) = 1$. Otherwise there is an astronomically slim chance that

```math
\Jacobi_N(y) = \Jacobi_N(w) = \Jacobi_N(z) = 0.
```

This renders the request invalid, so the server would disregard it since it fails the $\Jacobi_N(y) = -1$ check. This is so unlikely, however, that we will just take that chance on each request rather than spending time checking invertibility. As a practical matter, this will never happen—in the literal sense that for realistic values of $N$, it is so unlikely that even over centuries of trillions of clients using this protocol, chances are that none of them will ever choose $z$ that isn't in $\Z_N^*$. If it does happen, the consequences are virtually nil: a single request will be ignored when computing a client count estimate that is approximate anyway. In other words, invertibility is absolutely not worth checking for here.

Here's the actual code for this in our test implementation:
```julia
h = hash_resource_class(x₀, class)
x = modmul(x₀, powermod(g, h, N), N)
z = rand(rng, 1:N-1)
w = powermod(z, widen(B) << m, N)
i = rand(rng, 0:(1 << (m-1)) - 1)
t = widemul(2B, i) + 1
y = modmul(w, powermod(x, t, N), N)
```
As you can see, this is a very literal rendition of the steps described. On my M2 MacBook Air, for $1024$-bit $N$, this takes around $80$ microseconds; for $2048$-bit $N$, it takes $435$ microseconds, less than half a millisecond. In other words, these computations are easily fast enough to perform during each Pkg request — the network request itself will take orders of magnitude more time.
