# HyperLogLog Over RSA

This site is two things at once:

- a **reference implementation** — the [`HyperLogLogOverRSA.jl`](https://github.com/StefanKarpinski/HyperLogLogOverRSA.jl) Julia package, a working test/demo implementation of the protocol; and
- a **writeup** — a from-scratch explanation of the *motivation*, *design*, and *security arguments* behind it.

**HyperLogLog Over RSA** is a protocol for estimating the number of unique clients making requests to a service — for instance, the number of distinct installs hitting an open-source project's package servers — *without* anyone, even someone with full access to the server logs, being able to identify or track individual clients. It combines HyperLogLog cardinality estimation with the algebraic structure of RSA rings, so that clients sample encrypted HyperLogLog values they cannot themselves decode or bias.

## Reading the writeup

For the design and the privacy model:

- [Executive Summary](executive-summary.md) — the short version.
- The writeup proper, in order, starts at [Anonymously Counting Users](01-counting-users.md) and builds to the [Protocol Summary](13-protocol-summary.md).
- [For Cryptographers](for-cryptographers.md) — a condensed, notation-first tour for readers who want to check the math.

## Using the package

For running the protocol:

- [Reference Implementation](reference.md) — installation, a usage example, and the full API.

```julia
using HyperLogLogOverRSA
ring   = Ring(2^12-1, 63, 1024)               # generate a HyperLogLog RSA ring
cert   = RingCert(ring)                        # certificate that the ring is fingerprint-free
client = Client(cert)                          # verify cert, pick a random secret
y      = hll_generate(client, "/package/123")  # encrypted HLL token for a request
hll_decode(ring, y)                            # server-side: (bucket, geometric) HLL value
```

!!! note
    The implementation is a test/demo of the protocol — intended to make the design concrete and checkable, not a hardened production library.
