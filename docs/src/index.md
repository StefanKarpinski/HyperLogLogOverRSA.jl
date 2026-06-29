# HyperLogLog Over RSA

**HyperLogLog Over RSA** is a new protocol for estimating the number of unique clients making requests to a service — for instance, the number of distinct installs hitting an open-source project’s package servers — *without* anyone, even someone with full access to the server logs, being able to identify or track individual clients. It combines HyperLogLog cardinality estimation with the algebraic structure of RSA rings, so that clients sample encrypted HyperLogLog values they cannot themselves decode or bias.

This site includes:

- A **writeup** — a from-scratch explanation of the *motivation*, *design*, and *security arguments* behind it;
- A **reference implementation** — the [`HyperLogLogOverRSA.jl`](https://github.com/StefanKarpinski/HyperLogLogOverRSA.jl) Julia package, a working test/demo implementation of the protocol.

## Reading the writeup

For the design and the privacy model:

- [Executive Summary](01-executive-summary.md) — the short version;
- [For Cryptographers](02-for-cryptographers.md) — a condensed, notation-first tour for readers who want to check the math;
- The full writeup starts [here](03-counting-users.md) and proceeds in short chapters.

## Using the package

For running the protocol:

- [Reference Implementation](reference.md) — installation, a usage example, and the full API.

```julia
using HyperLogLogOverRSA
ring   = Ring(2^12-1, 63, 1024)                # generate a HyperLogLog RSA ring
cert   = RingCert(ring)                        # certificate that the ring is fingerprint-free
client = Client(cert)                          # verify cert, pick a random secret
y      = hll_generate(client, "/package/123")  # encrypted HLL token for a request
hll_decode(ring, y)                            # server-side: (bucket, geometric) HLL value
```

!!! note
    The implementation is a test/demo of the protocol — intended to make the design concrete and checkable, not a hardened production library.
