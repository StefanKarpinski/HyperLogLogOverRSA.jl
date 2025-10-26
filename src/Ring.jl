using Random
using Primes

## Representing & generating RSA Ring instances

"""
Ring structure:

    N = P Q = (2 B p + 1)(2^m q + 1)

    ℤ_N ≅ ℤ_P × ℤ_Q
    ℤ_P* ≅ C_2 × C_B × C_p
    ℤ_Q* ≅ C_(2^m) × C_q

"""
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
    # argument checks
    isodd(B) || throw(ArgumentError("B must be odd"))
    m > 0 || throw(ArgumentError("m must be positive"))
    L > 0 || throw(ArgumentError("L must be positive"))

    # range of N values
    N_max = one(T) << L - 1
    N_min = one(T) << (L-1) + 1

    # coefficients
    C, D = 2T(B), one(T) << m
    swap = rand(rng, Bool) # generate left or right first?
    if swap
        C, D = D, C
    end

    # generate (P, p) pair
    P_max = one(T) << cld(L,2) - 1
    P_min = one(T) << fld(L-1,2) + 1
    local P, p
    while true
        P, p = gen_prime_pair(P_min, P_max, C; rng)
        B ∉ (P, p) && break
    end

    # generate (Q, q) pair
    Q_min = cld(N_min, P)
    Q_max = fld(N_max, P)
    local Q, q
    while true
        Q, q = gen_prime_pair(Q_min, Q_max, D; rng)
        Q ∉ (B, P, p) && q ∉ (B, P, p) && break
    end

    # swap primes if coefficients were swapped
    if swap
        p, q = q, p
    end

    # construct Ring object
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

Base.getproperty(ring::Ring, name::Symbol) =
    name === :N ? ring.P*ring.Q :
    name === :P ? 2*ring.p*ring.B + true :
    name === :Q ? ring.q << ring.m + true :
    name === :λ ? ring.B*ring.p*(ring.q << ring.m) :
        getfield(ring, name)

modulus(ring::Ring) = ring.N
factors(ring::Ring) = (ring.P, ring.Q)
lambda(ring::Ring) = ring.λ

# don't print prime factors to avoid accidentally leaking them
Base.show(io::IO, ring::Ring) =
    print(io, "Ring(B=$(ring.B), m=$(ring.m), N=$(ring.N))")

function find_g(ring::Ring)
    P, Q = factors(ring)
    # find generator for ℤ_P*
    local g_P
    range_P = 1:P-1
    λ_P = 2*ring.B*ring.p
    while true
        g_P = rand(range_P)
        powermod(g_P, ring.B*ring.p, P) ≠ 1 &&
        powermod(g_P, 2*ring.B, P) ≠ 1 &&
        all(powermod(g_P, λ_P ÷ r, P) ≠ 1
            for (r, _) in factor(ring.B)) && break
    end
    # find generator for ℤ_Q*
    local g_Q
    range_Q = 1:Q-1
    while true
        g_Q = rand(range_Q)
        powermod(g_Q, ring.q << (ring.m-1), Q) ≠ 1 &&
        powermod(g_Q, one(Q) << ring.m, Q) ≠ 1 && break
    end
    # combine into g ∈ ℤ_N*
    _, u, v = gcdx(P, Q)
    g = mod(g_P*v*Q + g_Q*u*P, P*Q)
    return g
end

function find_x(ring::Ring)
    N = ring.N
    range = 0:N-1
    while true
        x = rand(range)
        jacobi(x, N) == -1 && return x
    end
end

bucket_map(ring::Ring{T}) where {T<:Integer} =
    Dict(powermod(one(T) << b, 2ring.p, ring.P) => b for b = 0:ring.B.-1)

function decode_bucket(
    ring :: Ring{T},
    y    :: Integer;
    bmap :: Dict{T,Int} = bucket_map(ring),
) where {T<:Integer}
    bmap[powermod(y, 2ring.p, ring.P)]
end

function decode_sample(
    ring :: Ring{T},
    y    :: T,
) where {T<:Integer}
    z = powermod(y, ring.q, ring.Q) # z = y^q mod Q
    k = ring.m
    while !isone(z)
        z = powermod(z, 2, ring.Q) # z <- z^2 mod Q
        k -= 1
    end
    return k
end

## Generating prime pairs

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
