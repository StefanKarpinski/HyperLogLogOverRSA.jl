const B_max = 2^12
const m_max = 128         # per-request arithmetic runs in the ring's integer type
const L_max = 2^20
const α_min = exp2(112)

"""
    Client(cert::RingCert; rng=…) -> Client

Verify a published [`RingCert`](@ref) and, if every check passes, construct a
client holding the public ring parameters together with a freshly chosen random
secret `x₀` (a Jacobi "twist" element). Throws an `ArgumentError` if the
certificate fails any check. `rng` (default a `RandomDevice`) is the source of
randomness for `x₀`; pass a seeded RNG to get a reproducible client.

Call [`hll_generate`](@ref) to produce the encrypted HyperLogLog token to send
with a request.

In this reference implementation the `Client` object *is* the client's persisted
identity: each instance holds its own random `x₀`, so two instances are two
distinct clients. Persisting a client across sessions would just mean serializing
the object, which is out of scope here.
"""
struct Client{T<:Integer}
    B :: Int # bucket factor (odd)
    m :: Int # max geometric sample size
    N :: T   # ring modulus
    f :: T   # per-class generator, derived deterministically from N (semisharding)
    x₀ :: T  # client-specific random Jacobi twist element
end

function Client(
    B :: Int,
    m :: Int,
    N :: T,
    f :: T;
    rng :: AbstractRNG = DEFAULT_RNG,
) where {T<:Integer}
    x₀ = rand_jacobi_twist(N; rng)
    Client(B, m, N, f, x₀)
end

function Client(cert::RingCert; rng::AbstractRNG = DEFAULT_RNG)
    B, m, N = cert.B, cert.m, cert.N

    # check shape parameters
    B ≤ B_max ||
        throw(ArgumentError("cert: B too large: $B"))
    isodd(B) ||
        throw(ArgumentError("cert: B even: $B"))
    m ≤ m_max ||
        throw(ArgumentError("cert: m too large: $m"))
    m ≥ 3 ||
        throw(ArgumentError("cert: m < 3: $m"))

    # check modulus properties
    Base.top_set_bit(N) ≤ L_max ||
        throw(ArgumentError("cert: N too large: $N"))
    mod8(N) == 5 ||
        throw(ArgumentError("cert: N ≠ 5 mod 8: $N"))
    gcd(B, N) == 1 ||
        throw(ArgumentError("cert: gcd(B, N) ≠ 1: $N"))
    gcd(B, N-1) == 1 ||
        throw(ArgumentError("cert: gcd(B, N-1) ≠ 1: $N"))

    # check that cert contains enough square roots
    (8/5)^length(cert.sqrts) ≥ α_min ||
        throw(ArgumentError("cert: too few sqrts: $(length(cert.sqrts))"))

    # check provided square roots
    for (i, r) in enumerate(cert.sqrts)
        r² = powermod(r, 2, N)
        x = hash_into_J₊(N, :sqrt_x, i)
        x == r² && continue
        y = hash_into_J₊(N, :sqrt_y, i)
        y == r² && continue
        z = modmul(x, y, N)
        z == r² && continue
        throw(ArgumentError("cert: invalid sqrt (N=$N)"))
    end

    # cert is valid, N is safe; recompute the canonical semisharding generator
    f = derive_f(N, B)
    Client(B, m, N, f; rng)
end

Base.show(io::IO, c::Client) =
    print(io, "Client(B=$(c.B), m=$(c.m), N=$(c.N), x₀=$(c.x₀))")

"""
    hll_generate(client::Client, class="/registries"; rng=…) -> Integer

Produce a fresh, randomized encrypted HyperLogLog token `y = w xᵗ` for the given
resource `class`, to send along with a request. Every call re-randomizes the
token, so two tokens from the same client are unlinkable; yet they all decode
(by the ring holder) to that client's single, stable HLL value for the class.
The client cannot itself decode or bias the value it samples.

This *semishards*: using the per-class generator `f` (whose `C_B` coordinate is
zero) the client's **bucket is its own `b₀` in every class**, while the
**geometric rank shards** per class. Only the rank — the part that carries the
rare, identifying values — is randomized across classes; the bucket is left
fixed on purpose (see the resource-class sharding discussion).

`rng` (default a `RandomDevice`) supplies the per-token randomness. It does not
affect the decoded value — only the token's unlinkable encoding — so seeding it
makes token output reproducible without changing what the token decodes to.
"""
function hll_generate(client::Client, class::Any="/registries"; rng::AbstractRNG = DEFAULT_RNG)
    B, m, N, f, x₀ = client.B, client.m, client.N, client.f, client.x₀
    h = hash_resource_class(x₀, class)         # h = H(x₀, class)
    x = modmul(x₀, powermod(f, h, N), N)       # x = x₀ f^h (semisharding)
    z = rand(rng, 1:N-1)                       # z ∈ [1, N)
    w = powermod(z, oftype(N, B) << (m-1), N)  # w = z^(B 2^(m-1))
    w = randsign(rng, w, N)                    # w = ±z^(B 2^(m-1))
    i = rand(rng, zero(N):(oftype(N, 1) << (m-1)) - one(N))  # i ∈ [0, 2^(m-1))
    t = 2 * oftype(N, B) * i + one(N)          # t ≡ 1 mod 2B
    y = modmul(w, powermod(x, t, N), N)        # y = w x^t
end
