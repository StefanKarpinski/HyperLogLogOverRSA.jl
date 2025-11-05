module HyperLogLogOverRSA

using SHA
using Base: SHA1
using Random
const rng = RandomDevice()

include("Jacobi.jl")
include("Ring.jl")
include("Certificate.jl")
include("Client.jl")

export
    Ring,
    RingCert,
    Client,
    rand_semigenerator,
    bucket_map,
    hll_generate,
    hll_decode

end # module HyperLogLogOverRSA
