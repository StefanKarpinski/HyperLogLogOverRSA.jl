using Primes

## Representing & generating RSA Ring instances

"""
    Ring(B, m, L; rng=…) -> Ring

A HyperLogLog RSA ring: a modulus `N` together with the secret structure that
lets its holder decode HyperLogLog values while no one else can.

The ring has the form

    N = P Q = (4 B p + 1)(2^m q + 1)

so that, multiplicatively,

    ℤ_N^* ≅ C_4 × C_B × C_(2^m) × C_(p q)

where `B` (which must be odd) is the number of HyperLogLog buckets and `m` is one
more than the maximum geometric sample value; `L` is the target bit-length of `N`.

The `C_4` factor carries no HyperLogLog value: the noise subgroup randomizes it,
so it contributes nothing a log holder can use. It is there for its arithmetic
alone — it makes `P ≡ 5 mod 8`, hence `N ≡ 5 mod 8`, hence `jacobi(-1, N) = 1`. Otherwise
`-1` is an element of `J_N^-` with known log-coordinates, which allows forging tokens.

The primes `P`, `Q`, `p`, `q` are secret — and are deliberately omitted when a
`Ring` is shown — so only a holder of the `Ring` can decode HLL values via
[`hll_decode`](@ref). Publish a [`RingCert`](@ref) so that clients can use the
ring without learning its factorization.

`rng` (any `AbstractRNG`, default a `RandomDevice`) is the source of randomness
for choosing the primes; pass a seeded RNG for a reproducible ring.

By default the ring is *certifiable*: it meets every criterion the protocol
needs beyond the basic shape. Currently that is one thing — its canonical
per-class generator (see [`Client`](@ref)) shards the geometric rank, which the
generator enforces by regenerating `N` until it holds. Pass
`certifiable = false` to skip these criteria and return any modulus of the right
shape — useful for structural work, where they play no role.
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
    B :: Integer, # bucket factor — must be odd
    m :: Integer, # max geometric sample size
    L :: Integer; # bit length of modulus
    rng :: AbstractRNG = DEFAULT_RNG,
    certifiable :: Bool = true,
) where {T<:Integer}
    !certifiable && return generate_ring(T, B, m, L, rng)
    for _ in 1:1000
        ring = generate_ring(T, B, m, L, rng)
        # certifiable ⟺ the canonical f is a shardable unit: a unit (jacobi(f, N)
        # ≠ 0 — otherwise the negligibly-rare non-unit :f hash would break tokens)
        # whose C_{2^m} coordinate is odd (rank shards) ⟺ jacobi(f, Q) == -1
        f = derive_f(ring.N, ring.B)
        jacobi(f, ring.N) != 0 && jacobi(f, ring.Q) == -1 && return ring
    end
    throw(ArgumentError("no modulus with a shardable f for spec (B=$B, m=$m, L=$L)"))
end

function generate_ring(::Type{T}, B::Integer, m::Integer, L::Integer, rng::AbstractRNG) where {T<:Integer}
    # argument checks
    isodd(B) || throw(ArgumentError("B must be odd"))
    m ≥ 3 || throw(ArgumentError("m must be ≥ 3"))
    L > 0 || throw(ArgumentError("L must be positive"))

    # range of N values
    N_max = one(T) << L - 1
    N_min = one(T) << (L-1) + 1

    # ranges of prime factors
    P_scale, Q_scale = T(4B), one(T) << m
    P_min = Q_min = one(T) << fld(L-1,2) + 1
    P_max = Q_max = fld(N_max, P_min)
    # for small log₂(N) we need to check for feasibility, which
    # also speeds up the search when there are few solutions
    # for large log₂(N) it's unnecessary (there are many solutions)
    # and inefficient since narrowing the prime range is very slow
    narrow_if_L_less = 128
    while true
        done = true
        if L < narrow_if_L_less
            P_min, P_max = narrow_prime_range(P_scale, P_min, P_max)
            P_min ≤ P_max || throw(ArgumentError("infeasible ring spec"))
        end
        Q_min′ = cld(N_min, P_max)
        Q_max′ = fld(N_max, P_min)
        if Q_min′ > Q_min; Q_min = Q_min′; done = false; end
        if Q_max′ < Q_max; Q_max = Q_max′; done = false; end
        if L < narrow_if_L_less
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
    if L < narrow_if_L_less
        # check that one of these is usable (unique primes)
        p_min = div(P_min - 1, P_scale)
        p_max = div(P_max - 1, P_scale)
        q_min = div(Q_min - 1, Q_scale)
        q_max = div(Q_max - 1, Q_scale)
        allunique([B_factors; P_min; Q_max; p_min; q_max]) ||
        allunique([B_factors; P_max; Q_min; p_max; q_min]) ||
            throw(ArgumentError("infeasible ring spec"))
        # giving up here is overly conservative, but we want to be sure
        # that some usable solution exists before we start sampling
        # otherwise we could end up looping forever
    end

    swap = false
    local P, Q, p, q
    while true
        if rand(rng, Bool)
            swap = !swap
            P_scale, Q_scale = Q_scale, P_scale
            P_min, Q_min = Q_min, P_min
            P_max, Q_max = Q_max, P_max
        end
        # generate (P, p) pair
        while true
            P, p = gen_prime_pair(P_scale, P_min, P_max; rng)
            allunique([B_factors; P; p]) && break
        end
        # range of second prime factor
        Q_min′ = max(Q_min, cld(N_min, P))
        Q_max′ = fld(N_max, P) # previous bound doesn't matter for max
        if L < narrow_if_L_less
            Q_min′, Q_max′ = narrow_prime_range(Q_scale, Q_min′, Q_max′)
        end
        # generate (Q, q) pair
        Q, q = gen_prime_pair(Q_scale, Q_min′, Q_max′; rng)
        allunique([B_factors; P; p; Q; q]) || continue
        break
    end
    # unswap if needed
    if swap
        p, q = q, p
    end

    # construct the Ring object
    Ring{T}(B, m, p, q)
end

ring_type(L::Integer) = L < 64 ? Int64 : L < 128 ? Int128 : BigInt

function Ring(
    B :: Integer, # bucket factor — must be odd
    m :: Integer, # max geometric sample size
    L :: Integer; # bit length of modulus
    rng :: AbstractRNG = DEFAULT_RNG,
    certifiable :: Bool = true,
)
    Ring{ring_type(L)}(B, m, L; rng, certifiable)
end

Base.getproperty(ring::Ring, name::Symbol) =
    name === :N ? ring.P*ring.Q :
    name === :P ? 4*ring.p*ring.B + 1 :
    name === :Q ? ring.q << ring.m + 1 :
    name === :λ ? ring.B*ring.p*(ring.q << ring.m) :
        getfield(ring, name)

modulus(ring::Ring) = ring.N
factors(ring::Ring) = (ring.P, ring.Q)
lambda(ring::Ring) = ring.λ

# don't print prime factors to avoid accidentally leaking them
Base.show(io::IO, ring::Ring) =
    print(io, "Ring(B=$(ring.B), m=$(ring.m), N=$(ring.N))")

function rand_semigenerator(ring::Ring; rng::AbstractRNG = DEFAULT_RNG)
    P, Q = factors(ring)
    # find generator for ℤ_P^*
    range_P = 1:P-1
    P_1 = 4*ring.B*ring.p # P - 1
    𝟜B = 4*ring.B
    P_1_rs = [P_1 ÷ r for r in keys(factor(ring.B))]
    local g_P
    while true
        g_P = rand(rng, range_P)
        powermod(g_P, P_1 ÷ 2, P) ≠ 1 && # generates C_4
        powermod(g_P, 𝟜B, P) ≠ 1 &&      # nontrivial in C_p
        all(powermod(g_P, P_1_r, P) ≠ 1 for P_1_r in P_1_rs) && break # generates C_B
    end
    @assert jacobi(g_P, P) == -1
    # find generator for ℤ_Q^*
    range_Q = 1:Q-1
    q𝟚ᵐ⁻¹ = ring.q << (ring.m-1)
    𝟚ᵐ = one(Q) << ring.m
    local g_Q
    while true
        g_Q = rand(rng, range_Q)
        powermod(g_Q, q𝟚ᵐ⁻¹, Q) ≠ 1 &&
        powermod(g_Q, 𝟚ᵐ, Q) ≠ 1 && break
    end
    @assert jacobi(g_Q, Q) == -1
    # combine into g ∈ ℤ_N^*
    N = P*Q
    _, u, v = gcdx(P, Q)
    uP = mod(widemul(u, P), N)
    vQ = mod(widemul(v, Q), N)
    g = oftype(N, mod(g_P*vQ + g_Q*uP, N))
    @assert mod(g, P) == g_P
    @assert mod(g, Q) == g_Q
    @assert jacobi(g, N) == 1
    return g
end

"""
    bucket_map(ring::Ring) -> Dict

Precompute the map from `C_B` representatives in `ℤ_P^*` to bucket indices
`0:B-1`, so that repeated [`hll_decode`](@ref) calls need not recompute it.
"""
function bucket_map(ring::Ring{T}) where {T<:Integer}
    # find first g_B that generates the C_B part of ℤ_P^*
    P = ring.P
    B = ring.B
    P_1 = P-1
    g_B = 2
    while g_B < P
        all(powermod(g_B, P_1 ÷ p, P) ≠ 1 for p in keys(factor(B))) && break
        g_B += 1
    end
    𝟜p = 4ring.p
    Dict(powermod(g_B, 𝟜p*b, P) => b for b = 0:B-1)
end

function hll_bucket(
    ring :: Ring{T},
    x    :: Integer;
    bmap :: Dict{T,Int} = bucket_map(ring),
) where {T<:Integer}
    bmap[powermod(x, 4ring.p, ring.P)]
end

function hll_geometric(ring::Ring, x::Integer)
    y = powermod(x, ring.q, ring.Q) # y = x^q mod Q
    iszero(y) && throw(ArgumentError("invalid x ∉ ℤ_N^*"))
    k = ring.m
    while !isone(y)
        y = powermod(y, 2, ring.Q) # y <- y^2 mod Q
        k -= 1
    end
    # The noise subgroup randomizes the C_2^m component's top bit, so the two top
    # rungs of the ladder are indistinguishable and the sample caps at m-1.
    return min(k, ring.m - 1)
end

"""
    hll_decode(ring::Ring, y; bmap=bucket_map(ring)) -> (bucket, geometric)

Decode an encrypted HyperLogLog token `y` to its `(bucket, geometric)` value,
using the secret factorization carried by `ring`. Pass a precomputed `bmap`
(from [`bucket_map`](@ref)) to amortize bucket decoding across many tokens.
"""
function hll_decode(
    ring :: Ring{T},
    x    :: Integer;
    bmap :: Dict{T,Int} = bucket_map(ring),
) where {T<:Integer}
    b = hll_bucket(ring, x; bmap)
    k = hll_geometric(ring, x)
    return b, k
end
