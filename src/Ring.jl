using Random
using Primes

## Representing & generating RSA Ring instances

struct Ring{T<:Integer}
    # general shape
    B :: Int # bucket factor (odd)
    m :: Int # max geometric sample size

    # specific values
    p :: T # 1st inner prime
    q :: T # 2nd inner prime
end

function Ring{T}(
    B :: Integer, # bucket factor — must be odd
    m :: Integer, # max geometric sample size
    L :: Integer; # bit length of modulus
    rng :: AbstractRNG = Random.GLOBAL_RNG,
) where {T<:Integer}
    m ≥ 1 || throw(ArgumentError("m must be ≥ 1"))
    isodd(B) || throw(ArgumentError("B must be odd"))
    N_min = (one(T) << (L-1)) + 1
    N_max = (one(T) << L) - 1
    C, D = 2T(B), one(T) << m
    swap = rand(rng, Bool) # generate left or right first?
    if swap
        C, D = D, C
    end
    P_min = (one(T) << (L÷2-1)) + 1
    P_max = (one(T) << (L÷2)) - 1
    P, p = gen_prime_pair(P_min, P_max, C; rng)
    Q_min = cld(N_min, P)
    Q_max = fld(N_max, P)
    Q, q = gen_prime_pair(Q_min, Q_max, D; rng)
    if swap
        p, q = q, p
    end
    Ring{T}(B, m, p, q)
end

ring_type(L::Integer) =
    L ≤ 8   ? UInt8   :
    L ≤ 16  ? UInt16  :
    L ≤ 32  ? UInt32  :
    L ≤ 64  ? UInt64  :
    L ≤ 128 ? UInt128 : BigInt

function Ring(
    B :: Integer, # bucket factor — must be odd
    m :: Integer, # max geometric sample size
    L :: Integer; # bit length of modulus
    rng :: AbstractRNG = Random.GLOBAL_RNG,
)
    Ring{ring_type(L)}(B, m, L; rng)
end

function factors(ring::Ring{T}) where {T<:Integer}
    P = 2T(ring.B) * ring.p + true
    Q = (one(T) << ring.m) * ring.q + true
    P, Q
end

modulus(ring::Ring) = prod(factors(ring))

# don't print prime factors to avoid accidentally leaking them
function Base.show(io::IO, ring::Ring)
    show(io, typeof(ring))
    print(io, "(")
    N = modulus(ring)
    show(io, "B=$(ring.B), m=$(ring.m), N=$N")
    print(io, ")")
end

## Generating prime factors for Ring

"""
    gen_prime_pair(
        P_min :: Integer,
        P_max :: Integer,
        scale :: Integer;
        rng :: AbstractRNG = Random.GLOBAL_RNG,
    ) -> P, p

Return a pair of primes, `(P, p)`, such that:

- `P in P_min:P_max`
- `P == scale*p + 1`

Requires that `scale` is even.
"""
function gen_prime_pair(
    P_min :: Integer,
    P_max :: Integer,
    scale :: Integer;
    rng :: AbstractRNG = Random.GLOBAL_RNG,
)
    iseven(scale) || throw(ArgumentError("scale factor must be even"))

    # range of p values
    P_min = nextprime(P_min)
    P_max = prevprime(P_max)
    # P == s*p + 1 <=> p == (P - 1)/s
    p_min = cld(P_min - true, scale)
    p_max = fld(P_max - true, scale)
    # adjust to actual primes
    p_min = nextprime(p_min)
    p_max = prevprime(p_max)
    p_range = p_min:p_max

    isempty(p_range) && throw(ArgumentError("infeasible range"))

    while true
        # sample an odd value
        p = rand(rng, p_range) | true
        # check if it's good
        isprime(p) || continue
        P = scale*p + true
        isprime(P) || continue
        return P, p
    end
end
