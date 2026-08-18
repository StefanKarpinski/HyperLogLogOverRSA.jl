using SHA
using Base: SHA1
using Base.GMP.MPZ: mul_2exp!, add_ui!

function shift_add_byte!(x::BigInt, b::UInt8)
    mul_2exp!(x, 8)
    add_ui!(x, b)
end

# Hash a tuple of keys deterministically into an element of ℤ_N. The digest
# stream is accumulated at full precision (BigInt) so the result depends only on
# the *value* of N, not its in-memory integer type; a fixed-width accumulator
# would truncate the fold and make the same numeric N hash differently as
# Int64/Int128 vs BigInt (e.g. across a TOML reparse).
function hash_into_ring(N::Integer, keys::Union{Integer,AbstractString,Symbol}...)
    prefix = sprint() do io
        print(io, N)
        for key in keys
            T = key isa Integer ? "int" :
            key isa AbstractString ? "str" : "sym"
            print(io, '\0', T, string(key))
        end
    end
    L = Base.top_set_bit(N) + 1
    x = zero(BigInt)
    for i = 1:cld(L, 512)
        for b in sha512("$prefix\0$i\0")
            x = shift_add_byte!(x, b)
        end
    end
    return oftype(N, mod(x, N))
end

# Hash into J_N^+ (positive Jacobi symbol), as the certificate challenge scheme
# requires. A hashed element already has positive Jacobi symbol half the time;
# when it is negative we multiply by 2, which lies in J_N^- whenever N ≡ 5 mod 8
# (the ring shape guarantees this), mapping it into J_N^+. A zero Jacobi symbol
# means the hash shares a factor with N — so N is factorable — and we reject N
# rather than resample: for a valid semiprime this never happens (density
# ~2/√N), and when it does the modulus is broken.
function hash_into_J₊(N::Integer, keys::Union{Integer,AbstractString,Symbol}...)
    x = hash_into_ring(N, keys...)
    j = jacobi(x, N)
    j == 0 && throw(ArgumentError("hash is a non-unit ⇒ N is factorable (N=$N)"))
    return j == -1 ? oftype(N, mod(widemul(oftype(N, 2), x), N)) : x
end

# Per-class exponent: HMAC-SHA-256 keyed by a hash of the client's secret x₀,
# with the first 16 bytes read big-endian into a 128-bit integer. Keying by x₀
# means the client can't bias its own draw and each class is an independent
# draw. The ring holder never recomputes this — decoding uses the factorization
# directly — so any keyed hash of the class works; HMAC is a standard keyed MAC.
function hash_resource_class(x₀::Integer, class::Any)
    key = sha256(digits(UInt8, x₀; base=256))
    bytes = hmac_sha2_256(key, codeunits(string(class)))
    h = zero(UInt128)
    for i = 1:sizeof(h)
        h <<= 8
        h |= bytes[i]
    end
    return h
end
