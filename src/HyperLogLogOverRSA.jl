module HyperLogLogOverRSA

using Random
const rng = RandomDevice()

include("PrimePairs.jl")
include("Jacobi.jl")
include("Hashing.jl")
include("Ring.jl")
include("RingCert.jl")
include("Client.jl")

export
    Ring,
    RingCert,
    Client,
    bucket_map,
    hll_generate,
    hll_decode

end # module HyperLogLogOverRSA
