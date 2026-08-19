using SHA
using Base: SHA1
using Base.GMP.MPZ: mul_2exp!, add_ui!

function shift_add_byte!(x::BigInt, b::UInt8)
    mul_2exp!(x, 8)
    add_ui!(x, b)
end

# Extra bits generated beyond the bit-length of N before reducing mod N, so the
# reduction is nearly uniform. Reducing a uniform draw from [0, 2^H) mod N has
# relative per-residue bias < N/2^H < 2^(top_set_bit(N) - H), so H - top_set_bit(N)
# is the "bits of uniformity"; 128 keeps the bias (≤ 2^-128) negligible against
# even the strongest certificate strength we target. Generating only enough bits
# to cover N (the earlier behavior) leaves ~1 bit of headroom — a ~2× bias —
# whenever top_set_bit(N) sits just under a 512-bit block boundary (e.g. a
# 2047-bit modulus, a common size for a product of two 1024-bit primes).
const HASH_MARGIN = 128

# Number of 512-bit SHA-512 blocks folded to hash into ℤ_N with that margin.
hash_blocks(N::Integer) = cld(Base.top_set_bit(N) + HASH_MARGIN, 512)

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
    x = zero(BigInt)
    for i = 1:hash_blocks(N)
        for b in sha512("$prefix\0$i\0")
            x = shift_add_byte!(x, b)
        end
    end
    return oftype(N, mod(x, N))
end

# Hash into ℤ_N \ J_N^-: a value whose Jacobi symbol is not -1 — so in J_N^+, or
# (with probability ~2/√N) a non-unit. A raw hash already has Jacobi symbol ≠ -1
# half the time; when it is -1 we multiply by 2, which flips it into J_N^+ because
# 𝒥_N(2) = -1 for N ≡ 5 mod 8 (the ring shape guarantees this). This is a *total*
# function — it never fails — and that is safe: the certificate's covering lemma
# (one of x, y, xy is a quadratic residue) holds for the whole set {𝒥_N ≠ -1} on a
# squarefree two-prime N, not just J_N^+, and push_sqrt_mod_N handles a non-unit
# fine (modsqrt(0, P) = 0). A non-unit would expose a factor of N, so it cannot
# arise for a valid semiprime; callers use the result without checking.
function hash_into_J₊(N::Integer, keys::Union{Integer,AbstractString,Symbol}...)
    x = hash_into_ring(N, keys...)
    jacobi(x, N) == -1 ? oftype(N, mod(widemul(oftype(N, 2), x), N)) : x
end

# Derive the per-class generator f. It is a deterministic
# function of N (recomputable by every client, so a malicious server cannot vary
# it per client as a tracking tag): the B-th power of a hashed J_N^+ element.
# Being a B-th power, f's C_B coordinate is zero, so x₀ f^h keeps the client's
# own bucket b₀ fixed across every resource class while the rank still shards;
# and 𝒥_N(f) = 𝒥_N(hash)^B = +1, so x₀ f^h stays in J_N^- for every h. The ring
# generator additionally checks that f is a unit (jacobi(f, N) ≠ 0 — a non-unit :f
# hash, negligibly rare, would break token generation) and that its C_{2^m}
# coordinate is odd (jacobi(f, Q) == -1, needed for the rank to shard),
# regenerating N otherwise. Given a certified ring, a client uses f as-is.
derive_f(N::Integer, B::Integer) = powermod(hash_into_J₊(N, :f), oftype(N, B), N)

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
