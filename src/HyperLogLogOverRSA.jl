module HyperLogLogOverRSA

using Random

# Default source of randomness: cryptographic entropy. Every randomized function
# takes an `rng::AbstractRNG` keyword defaulting to this, so callers (e.g. tests)
# can pass a seeded rng for reproducibility. Because each method uses its `rng`
# argument directly in `rand`, Julia specializes on the concrete rng type — no
# dynamic dispatch on the hot path.
const DEFAULT_RNG = RandomDevice()

include("PrimePairs.jl")
include("Jacobi.jl")
include("Hashing.jl")
include("Ring.jl")
include("RingCert.jl")
include("Client.jl")
include("Estimate.jl")

export
    Ring,
    RingCert,
    Client,
    bucket_map,
    hll_generate,
    hll_decode,
    hll_estimate

end # module HyperLogLogOverRSA
