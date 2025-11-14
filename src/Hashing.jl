using SHA
using Base: SHA1
using Base.GMP.MPZ: mul_2exp!, add_ui!

function shift_add_byte!(x::Integer, b::UInt8)
    x <<= 8
    x |= b
end

function shift_add_byte!(x::BigInt, b::UInt8)
    mul_2exp!(x, 8)
    add_ui!(x, b)
end

function ring_hash(
    N :: Integer,
    keys :: Union{Integer,AbstractString,Symbol}...;
    untwist :: Integer = zero(N),
)
    prefix = sprint() do io
        print(io, N)
        for key in keys
            T = key isa Integer ? "int" :
            key isa AbstractString ? "str" : "sym"
            print(io, '\0', T, string(key))
        end
    end
    L = Base.top_set_bit(N) + 1
    x = zero(N)
    for i = 1:cld(L, 512)
        for b in sha512("$prefix\0$i\0")
            x = shift_add_byte!(x, b)
        end
    end
    x = mod(x, N)
    if !iszero(untwist) && jacobi(x, N) == -1
        x = mod(widemul(untwist, x), N)
    end
    return x
end

function hash_resource_class(salt::Any, class::Any)
    bytes = sha256("$salt\0$class\0")
    h = zero(UInt128)
    for i = 1:8
        h <<= 8
        h |= bytes[i]
    end
    return h
end
