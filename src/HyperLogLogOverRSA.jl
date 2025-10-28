module HyperLogLogOverRSA

include("Jacobi.jl")
include("Ring.jl")

using SHA

function hash_resource_class(str::AbstractString)
    h = BigInt()
    bytes = sha224(str)
    ccall((:__gmpz_import, Base.GMP.libgmp), Cvoid,
        (Base.GMP.MPZ.mpz_t, Csize_t, Cint, Csize_t, Cint, Csize_t, Ptr{UInt8}),
        h, length(bytes), 1, 1, 1, 0, bytes)
        # order=+1, size=1, endian=+1 (big), nails=0
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
