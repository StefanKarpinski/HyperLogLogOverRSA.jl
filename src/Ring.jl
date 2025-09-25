using Random
using Primes

## Representing & generating RSA Ring instances

struct Ring{T<:Integer}
    # general shape
    m :: Int # max geometric sample size
    W :: Int # work factor (prime)
    B :: Int # bucket factor (prime)

    # specific values
    p :: T # 1st inner prime
    q :: T # 2nd inner prime
end

function Ring{T}(
    m :: Integer, # max geometric sample size
    W :: Integer, # work factor — use next prime
    B :: Integer, # bucket factor — use next prime
    L :: Integer; # bit length of modulus
    rng :: AbstractRNG = Random.GLOBAL_RNG,
) where {T<:Integer}
    m ≥ 1 || throw(ArgumentError("m must be ≥ 1"))
    W = nextprime(W)
    B = nextprime(B)
    𝟚 = T(2)
    N_min = 𝟚^(L-1)+1
    N_max = 𝟚^L-1
    swap = rand(rng, Bool) # generate left or right first?
    C, D = swap ? (B, W) : (W, B)
    P_min = 𝟚^(L÷2-1)+1
    P_max = 𝟚^(L÷2)-1
    P, p = gen_prime_pair(P_min, P_max, T(C) << m; rng)
    Q_min = cld(N_min, P)
    Q_max = fld(N_max, P)
    Q, q = gen_prime_pair(Q_min, Q_max, T(D) << m; rng)
    swap && ((p, q) = (q, p))
    Ring{T}(m, W, B, p, q)
end

function factors(ring::Ring{T}) where {T<:Integer}
    M = T(2)^ring.m
    P = T(M*ring.W)*ring.p + true
    Q = T(M*ring.B)*ring.q + true
    P, Q
end

modulus(ring::Ring) = prod(factors(ring))

# don't print prime factors to avoid accidentally leaking them
function Base.show(io::IO, ring::Ring)
    show(io, typeof(ring))
    print(io, "(")
    N = modulus(ring)
    show(io, "m=$(ring.m), W=$(ring.W), B=$(ring.B), N=$N")
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
