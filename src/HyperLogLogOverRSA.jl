module HyperLogLogOverRSA

using Random
const rng = RandomDevice()

include("Jacobi.jl")
include("Ring.jl")
include("Client.jl")

export
    Ring,
    Client,
    rand_semigenerator,
    bucket_map,
    hll_generate,
    hll_decode

end # module HyperLogLogOverRSA
