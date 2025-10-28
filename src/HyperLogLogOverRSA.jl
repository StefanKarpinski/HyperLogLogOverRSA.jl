module HyperLogLogOverRSA

include("Jacobi.jl")
include("Ring.jl")

using SHA

function hash_resource_class(str::AbstractString)
    bytes = sha224(str)
    h = zero(UInt64)
    for i = 1:8
        h <<= 8
        h |= bytes[i]
    end
    return h
end

export
    Ring,
    rand_elt,
    rand_power,
    rand_semigenerator,
    rand_jacobi_twist,
    hash_resource_class

end # module HyperLogLogOverRSA
