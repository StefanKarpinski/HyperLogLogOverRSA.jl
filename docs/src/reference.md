# Reference Implementation

The [`HyperLogLogOverRSA.jl`](https://github.com/StefanKarpinski/HyperLogLogOverRSA.jl) package is a test/demo implementation of the full protocol described in the writeup. It exists to make the design concrete and checkable — not to be a hardened production library.

## Installation

```julia-repl
pkg> add https://github.com/StefanKarpinski/HyperLogLogOverRSA.jl
```

## Usage

```julia-repl
julia> using HyperLogLogOverRSA

julia> ring = Ring(2^12-1, 63, 1024)   # generate a HyperLogLog RSA ring
Ring(B=4095, m=63, N=…)

julia> cert = RingCert(ring)           # certificate that the ring is fingerprint-free
RingCert(B=4095, m=63, N=…)

julia> client = Client(cert)           # check the certificate, pick a random secret x₀
Client(B=4095, m=63, N=…, x₀=…)

julia> y = hll_generate(client)        # encrypted HLL token, default resource class
…

julia> hll_decode(ring, y)             # ring holder decodes the (bucket, geometric) value
(1342, 0)

julia> y = hll_generate(client)        # a fresh, unlinkable token for the same client/class
…

julia> hll_decode(ring, y)             # …decodes to the same HLL value
(1342, 0)

julia> y = hll_generate(client, "/package/123")   # an independent value per resource class
…

julia> hll_decode(ring, y)
(1902, 3)
```

The server should log the `y` values and only decode them later, offline — the
factorization of `N` never needs to live on the public-facing servers. See the
[Protocol Summary](13-protocol-summary.md) for the end-to-end flow.

## API

```@docs
Ring
RingCert
Client
hll_generate
hll_decode
bucket_map
```
