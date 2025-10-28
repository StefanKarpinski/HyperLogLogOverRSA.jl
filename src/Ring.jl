using Random
using Primes

## Representing & generating RSA Ring instances

"""
Ring structure:

    N = P Q = (2 B p + 1)(2^m q + 1)

    ℤ_N ≅ ℤ_P × ℤ_Q
    ℤ_P^* ≅ C_2 × C_B × C_p
    ℤ_Q^* ≅ C_(2^m) × C_q

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

    # ranges of prime factors
    P_scale, Q_scale = T(2B), one(T) << m
    P_min = Q_min = one(T) << fld(L-1,2) + 1
    P_max = Q_max = fld(N_max, P_min)
    # for small log₂(N) we need to check for feasibility, which
    # also speeds up the search when there are few solutions
    # for large log₂(N) it's unnecessary (there are many solutions)
    # and inefficient since narrowing the prime range is very slow
    while true
        done = true
        if L < 128
            P_min, P_max = narrow_prime_range(P_scale, P_min, P_max)
            P_min ≤ P_max || throw(ArgumentError("infeasible ring spec"))
        end
        Q_min′ = cld(N_min, P_max)
        Q_max′ = fld(N_max, P_min)
        if Q_min′ > Q_min; Q_min = Q_min′; done = false; end
        if Q_max′ < Q_max; Q_max = Q_max′; done = false; end
        if L < 128
            Q_min, Q_max = narrow_prime_range(Q_scale, Q_min, Q_max)
            Q_min ≤ Q_max || throw(ArgumentError("infeasible ring spec"))
        end
        P_min′ = cld(N_min, Q_max)
        P_max′ = fld(N_max, Q_min)
        if P_min′ > P_min; P_min = P_min′; done = false; end
        if P_max′ < P_max; P_max = P_max′; done = false; end
        done && break
    end
    # feasible solutions exist:
    @assert N_min ≤ P_min*Q_max ≤ N_max
    @assert N_min ≤ P_max*Q_min ≤ N_max

    B_factors = sort!(collect(keys(factor(B))))
    if L < 128
        # check that one of these is usable (unique primes)
        p_min = div(P_min - 1, P_scale)
        p_max = div(P_max - 1, P_scale)
        q_min = div(Q_min - 1, P_scale)
        q_max = div(Q_max - 1, P_scale)
        allunique([B_factors; P_min; Q_max; p_min; q_max]) ||
        allunique([B_factors; P_max; Q_min; p_max; q_min]) ||
            throw(ArgumentError("infeasible ring spec"))
        # giving up here is overly conservative, but we want to be sure
        # that some usable solution exists before we start sampling
        # otherwise we could end up looping forever
    end

    swap = rand(rng, Bool) # generate left or right first?
    if swap
        P_scale, Q_scale = Q_scale, P_scale
        P_min, Q_min = Q_min, P_min
        P_max, Q_max = Q_max, P_max
    end

    local P, Q, p, q
    while true
        # generate (P, p) pair
        while true
            P, p = gen_prime_pair(P_scale, P_min, P_max; rng)
            allunique([B_factors; P; p]) && break
        end
        # range of second prime factor
        Q_min′ = max(Q_min, cld(N_min, P))
        Q_max′ = fld(N_max, P) # previous bound doesn't matter for max
        Q_min′, Q_max′ = narrow_prime_range(Q_scale, Q_min′, Q_max′)
        # generate (Q, q) pair
        Q, q = gen_prime_pair(Q_scale, Q_min′, Q_max′; rng)
        allunique([B_factors; P; p; Q; q]) && break
    end

    # swap primes if coefficients were swapped
    if swap
        p, q = q, p
    end

    # construct the Ring object
    Ring{T}(B, m, p, q)
end

ring_type(L::Integer) = L < 64 ? Int64 : L < 128 ? Int128 : BigInt

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
    name === :P ? 2*ring.p*ring.B + 1 :
    name === :Q ? ring.q << ring.m + 1 :
    name === :λ ? ring.B*ring.p*(ring.q << ring.m) :
        getfield(ring, name)

modulus(ring::Ring) = ring.N
factors(ring::Ring) = (ring.P, ring.Q)
lambda(ring::Ring) = ring.λ

# don't print prime factors to avoid accidentally leaking them
Base.show(io::IO, ring::Ring) =
    print(io, "Ring(B=$(ring.B), m=$(ring.m), N=$(ring.N))")

function rand_semigenerator(ring::Ring)
    P, Q = factors(ring)
    # find generator for ℤ_P^*
    range_P = 1:P-1
    Bp = ring.B*ring.p
    𝟚B = 2*ring.B
    λ_P = 2*ring.B*ring.p
    λ_P_rs = [λ_P ÷ r for r in keys(factor(ring.B))]
    local g_P
    while true
        g_P = rand(range_P)
        powermod(g_P, Bp, P) ≠ 1 &&
        powermod(g_P, 𝟚B, P) ≠ 1 &&
        all(powermod(g_P, λ_P_r, P) ≠ 1 for λ_P_r in λ_P_rs) && break
    end
    # find generator for ℤ_Q^*
    range_Q = 1:Q-1
    q𝟚ᵐ⁻¹ = ring.q << (ring.m-1)
    𝟚ᵐ = one(Q) << ring.m
    local g_Q
    while true
        g_Q = rand(range_Q)
        powermod(g_Q, q𝟚ᵐ⁻¹, Q) ≠ 1 &&
        powermod(g_Q, 𝟚ᵐ, Q) ≠ 1 && break
    end
    # combine into g ∈ ℤ_N^*
    _, u, v = gcdx(P, Q)
    g = mod(g_P*v*Q + g_Q*u*P, P*Q)
    return g
end

function rand_jacobi_twist(ring::Ring)
    N = ring.N
    range = 0:N-1
    while true
        τ = rand(range)
        jacobi(τ, N) == -1 && return τ
    end
end

rand_elt(ring::Ring) = rand(1:ring.N-1)

rand_power(B::Integer, N::Integer) = 2B*rand(0:(N-1)÷B-1) + 1
rand_power(ring::Ring) = rand_power(ring.B, ring.N)

function bucket_map(ring::Ring{T}) where {T<:Integer}
    # find first g_B that generates the C_B part of ℤ_P^*
    P = ring.P; Pm1 = P - 1
    B = ring.B
    g_B = 2
    while g_B < P
        all(powermod(g_B, Pm1 ÷ p, P) ≠ 1 for p in keys(factor(B))) && break
        g_B += 1
    end
    𝟚p = 2ring.p
    Dict(powermod(g_B, 𝟚p*b, P) => b for b = 0:B-1)
end

function hll_bucket(
    ring :: Ring{T},
    x    :: Integer;
    bmap :: Dict{T,Int} = bucket_map(ring),
) where {T<:Integer}
    bmap[powermod(x, 2ring.p, ring.P)]
end

function hll_geometric(ring::Ring, x::Integer)
    y = powermod(x, ring.q, ring.Q) # y = x^q mod Q
    iszero(y) && throw(ArgumentError("invalid x ∉ ℤ_N^*"))
    k = ring.m
    while !isone(y)
        y = powermod(y, 2, ring.Q) # y <- y^2 mod Q
        k -= 1
    end
    return k
end

function hll_value(
    ring :: Ring{T},
    x    :: Integer;
    bmap :: Dict{T,Int} = bucket_map(ring),
) where {T<:Integer}
    b = hll_bucket(ring, x; bmap)
    k = hll_geometric(ring, x)
    return b, k
end

## Generating prime pairs

# P == scale*p + 1 <=> p == (P - 1)/scale

function next_paired_prime(scale::Integer, P::Integer, P_max::Integer)
    while P ≤ P_max
        P = scale*nextprime(cld(P - 1, scale)) + 1
        isprime(P) && break
        P = nextprime(P)
    end
    if P ≤ P_max
        @assert isprime(P)
        @assert isprime((P-1) ÷ scale)
    end
    return P
end

function prev_paired_prime(scale::Integer, P_min::Integer, P::Integer)
    while P_min ≤ P
        P = scale*prevprime(fld(P - 1, scale)) + 1
        isprime(P) && break
        P = prevprime(P)
    end
    if P_min ≤ P
        @assert isprime(P)
        @assert isprime((P-1) ÷ scale)
    end
    return P
end

function narrow_prime_range(scale::Integer, P_min::Integer, P_max::Integer)
    iseven(scale) || throw(ArgumentError("scale factor must be even"))
    P_max = prev_paired_prime(scale, P_min, P_max)
    P_min = next_paired_prime(scale, P_min, P_max)
    return P_min, P_max
end

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
    scale :: Integer,
    P_min :: Integer,
    P_max :: Integer;
    rng :: AbstractRNG = Random.GLOBAL_RNG,
)
    iseven(scale) || throw(ArgumentError("scale factor must be even"))

    p_min = cld(P_min - 1, scale)
    p_max = fld(P_max - 1, scale)
    p_min += iseven(p_min)
    p_max -= iseven(p_max)
    p_range = p_min:p_max

    while true
        # sample an odd value
        p = rand(rng, p_range) | 1
        # check if it's good
        isprime(p) || continue
        P = scale*p + 1
        isprime(P) || continue
        return P, p
    end
end
