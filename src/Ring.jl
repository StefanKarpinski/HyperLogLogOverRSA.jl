using Primes

## Representing & generating RSA Ring instances

"""
    Ring(B, m, L; rng=…) -> Ring

A HyperLogLog RSA ring: a modulus `N` together with the secret structure that
lets its holder decode HyperLogLog values while no one else can.

The ring has the form

    N = P Q = (2 B p + 1)(2^m q + 1)

so that, multiplicatively,

    ℤ_N^* ≅ C_(2B) × C_(2^m) × C_(p q)

where `B` (which must be ≡ 2 mod 4) is the number of HyperLogLog buckets and `m`
is the maximum geometric sample value; `L` is the target bit-length of `N`.

`B ≡ 2 mod 4` makes `P ≡ 5 mod 8`, hence `N ≡ 5 mod 8`, hence `jacobi(-1, N) = 1`
— which keeps `-1` out of `J_N^-` and so denies a forger a public element with
known logarithms. See the security analysis for why that matters.

The primes `P`, `Q`, `p`, `q` are secret — and are deliberately omitted when a
`Ring` is shown — so only a holder of the `Ring` can decode HLL values via
[`hll_decode`](@ref). Publish a [`RingCert`](@ref) so that clients can use the
ring without learning its factorization.

`rng` (any `AbstractRNG`, default a `RandomDevice`) is the source of randomness
for choosing the primes; pass a seeded RNG for a reproducible ring.
"""
struct Ring{T<:Integer}
    # general shape
    B :: Int # bucket count (≡ 2 mod 4)
    m :: Int # max geometric sample size

    # specific values
    p :: T # 1st inner prime
    q :: T # 2nd inner prime
end

function Ring{T}(
    B :: Integer, # bucket count — must be ≡ 2 mod 4
    m :: Integer, # max geometric sample size
    L :: Integer; # bit length of modulus
    rng :: AbstractRNG = DEFAULT_RNG,
) where {T<:Integer}
    # argument checks
    mod(B, 4) == 2 || throw(ArgumentError("B must be ≡ 2 mod 4"))
    m ≥ 3 || throw(ArgumentError("m must be ≥ 3"))
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
    B :: Integer, # bucket count — must be ≡ 2 mod 4
    m :: Integer, # max geometric sample size
    L :: Integer; # bit length of modulus
    rng :: AbstractRNG = DEFAULT_RNG,
)
    Ring{ring_type(L)}(B, m, L; rng)
end

Base.getproperty(ring::Ring, name::Symbol) =
    name === :N ? ring.P*ring.Q :
    name === :P ? 2*ring.p*ring.B + 1 :
    name === :Q ? ring.q << ring.m + 1 :
    name === :λ ? (ring.B >> 1)*ring.p*(ring.q << ring.m) :
        getfield(ring, name)

modulus(ring::Ring) = ring.N
factors(ring::Ring) = (ring.P, ring.Q)
lambda(ring::Ring) = ring.λ

# don't print prime factors to avoid accidentally leaking them
Base.show(io::IO, ring::Ring) =
    print(io, "Ring(B=$(ring.B), m=$(ring.m), N=$(ring.N))")

function rand_semigenerator(ring::Ring; rng::AbstractRNG = DEFAULT_RNG)
    P, Q = factors(ring)
    # find generator for ℤ_P^*: its C_2B part must generate (checked prime by
    # prime of 2B, which for B ≡ 2 mod 4 covers the 2-part and the odd part
    # alike) and its C_p part must be nontrivial
    range_P = 1:P-1
    𝟚B = 2*ring.B
    P_1 = 𝟚B*ring.p # P - 1
    P_1_rs = [P_1 ÷ r for r in keys(factor(𝟚B))]
    local g_P
    while true
        g_P = rand(rng, range_P)
        powermod(g_P, 𝟚B, P) ≠ 1 &&
        all(powermod(g_P, P_1_r, P) ≠ 1 for P_1_r in P_1_rs) && break
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

Precompute the map from `C_(2B)` representatives in `ℤ_P^*` to their logarithms
`0:2B-1`, so that repeated [`hll_decode`](@ref) calls need not recompute it.
"""
function bucket_map(ring::Ring{T}) where {T<:Integer}
    # find first g that generates the C_2B part of ℤ_P^*
    P = ring.P
    𝟚B = 2*ring.B
    P_1 = P-1
    rs = collect(keys(factor(𝟚B)))
    g = 2
    while g < P
        all(powermod(g, P_1 ÷ r, P) ≠ 1 for r in rs) && break
        g += 1
    end
    γ = powermod(oftype(P, g), ring.p, P) # kills C_p, generates C_2B
    Dict(powermod(γ, a, P) => a for a = 0:𝟚B-1)
end

# The C_2B logarithm of `x`, up to a fixed odd scalar fixed by `bmap`'s choice of
# generator — which only permutes bucket labels.
function hll_log(
    ring :: Ring{T},
    x    :: Integer;
    bmap :: Dict{T,Int} = bucket_map(ring),
) where {T<:Integer}
    bmap[powermod(x, ring.p, ring.P)]
end

# Returns the geometric sample `k` together with the order-4 element of ⟨x^q⟩,
# or zero when `x` has no such element (`k ≥ m-1`). Squaring down to 1 counts the
# order of the C_2^m part, and the second-to-last value it passes through is
# ι^u, where u is the leading unit of that part's logarithm and ι is one of the
# two order-4 elements mod Q — so the walk hands us `u` for free.
function hll_geometric(ring::Ring, x::Integer)
    Q = ring.Q
    y = powermod(x, ring.q, Q) # y = x^q mod Q
    iszero(y) && throw(ArgumentError("invalid x ∉ ℤ_N^*"))
    k = ring.m
    o2 = o4 = zero(y)
    while !isone(y)
        o4, o2 = o2, y
        y = powermod(y, 2, Q) # y <- y^2 mod Q
        k -= 1
    end
    return k, o4
end

"""
    hll_decode(ring::Ring, y; bmap=bucket_map(ring)) -> (bucket, geometric)

Decode an encrypted HyperLogLog token `y` to its `(bucket, geometric)` value,
using the secret factorization carried by `ring`. Pass a precomputed `bmap`
(from [`bucket_map`](@ref)) to amortize bucket decoding across many tokens.

The bucket index splits as `β + (B/2)·s`: `β` is the token's logarithm in the
odd part of `C_(2B)`, and `s` is the one bit of the 2-part that survives token
re-randomization. Neither `w` nor an odd `t` can touch 2-torsion, so `s` is
information the ring holder can read either way; it is part of the HyperLogLog
value rather than a fingerprint precisely because it is declared here. It is
pinned to `0` for `k ≥ m-1`, where the two orbits merge and only `B/2` buckets
are reachable — a `2^-(m-1)` share of tokens.
"""
function hll_decode(
    ring :: Ring{T},
    x    :: Integer;
    bmap :: Dict{T,Int} = bucket_map(ring),
) where {T<:Integer}
    a = hll_log(ring, x; bmap)
    k, o4 = hll_geometric(ring, x)
    # α = a mod 4 and the leading unit u of the C_2^m logarithm are each defined
    # only up to an odd scalar, but t is odd so t^2 = 1 mod 4: the product αu is
    # invariant under re-randomization even though neither factor is.
    s = 0
    if !iszero(o4)
        u = 2*o4 < ring.Q ? 1 : 3 # canonical order-4 element gets u = 1
        s = (a % 4)*u % 4 ≥ 2 ? 1 : 0
    end
    B′ = ring.B >> 1
    return a % B′ + B′*s, k
end
