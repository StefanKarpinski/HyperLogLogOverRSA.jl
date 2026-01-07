const B_max = 2^16
const m_max = 128
const L_max = 2^20
const α_min = exp2(50)

struct Client{T<:Integer}
    B :: Int # bucket factor (odd)
    m :: Int # max geometric sample size
    N :: T   # ring modulus
    g :: T   # server-provided, common ring semigenerator
    x₀ :: T  # client-specific random Jacobi twist element
end

function Client(
    B :: Int,
    m :: Int,
    N :: T,
    g :: T,
) where {T<:Integer}
    x₀ = rand_jacobi_twist(N)
    Client(B, m, N, g, x₀)
end

function Client(cert::RingCert)
    B, m, N, g = cert.B, cert.m, cert.N, cert.g

    # check shape parameters
    B ≤ B_max ||
        throw(ArgumentError("cert: B too large: $B"))
    isodd(B) ||
        throw(ArgumentError("cert: B even: $B"))
    m ≤ m_max ||
        throw(ArgumentError("cert: m too large: $m"))
    1 < m ||
        throw(ArgumentError("cert: m ≤ 1: $m"))

    # check modulus properties
    Base.top_set_bit(N) ≤ L_max ||
        throw(ArgumentError("cert: N too large: $N"))
    mod4(N) == 3 ||
        throw(ArgumentError("cert: N ≠ 3 mod 4: $N"))
    gcd(B, N) == 1 ||
        throw(ArgumentError("cert: gcd(B, N) ≠ 1: $N"))
    gcd(B, N-1) == 1 ||
        throw(ArgumentError("cert: gcd(B, N-1) ≠ 1: $N"))

    # check semigenerator Jacobi symbol
    jacobi(g, N) == 1 ||
        throw(ArgumentError("cert: invalid semigenerator: $g"))

    # check that cert contains enough square roots
    (8/5)^length(cert.sqrts) ≥ α_min ||
        throw(ArgumentError("cert: too few sqrts: $(length(cert.sqrts))"))

    # check provided square roots
    τ = fixed_twist(N)
    for (i, r) in enumerate(cert.sqrts)
        r² = powermod(r, 2, N)
        x = ring_hash(N, :sqrt_x, i; untwist=τ)
        x == r² && continue
        y = ring_hash(N, :sqrt_y, i; untwist=τ)
        y == r² && continue
        z = modmul(x, y, N)
        z == r² && continue
        throw(ArgumentError("cert: invalid sqrt (N=$N)"))
    end

    # cert is valid, N is safe
    Client(B, m, N, g)
end

Base.show(io::IO, c::Client) =
    print(io, "Client(B=$(c.B), m=$(c.m), N=$(c.N), x₀=$(c.x₀))")

function hll_generate(client::Client, class::Any="/registries")
    N, B, m, g, x₀ = client.N, client.B, client.m, client.g, client.x₀
    h = hash_resource_class(x₀, class)        # h = H(x₀, class)
    x = modmul(x₀, powermod(g, h, N), N)      # x = x₀ g^h
    z = rand(rng, 1:N-1)                      # z ∈ [1, N)
    w = powermod(z, widen(B) << m, N)         # w = z^(B 2^m)
    i = rand(rng, 0:(Int64(1) << (m-1)) - 1)  # i ∈ [0, 2^(m-1))
    t = widemul(2B, i) + 1                    # t = 1 mod 2B
    y = modmul(w, powermod(x, t, N), N)       # y = w x^t
end
