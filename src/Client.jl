struct Client{T<:Integer}
    B :: Int # bucket factor (odd)
    m :: Int # max geometric sample size
    N :: T   # ring modulus
    g :: T   # server-provided, common ring semigenerator
    x₀ :: T  # client-specific random Jacobi twist element
    salt :: SHA1 # pre-computed hash of x₀
end

function Client(
    B :: Int,
    m :: Int,
    N :: T,
    g :: T,
) where {T<:Integer}
    x₀ = rand_jacobi_twist(N)
    salt = SHA1(sha1(string(x₀, base=62)))
    Client(B, m, N, g, x₀, salt)
end

function Client(cert::RingCert)
    N = cert.N

    # check shape parameters
    isodd(cert.B) ||
        throw(ArgumentError("cert: B even: $(cert.B))"))
    cert.m > 1 ||
        throw(ArgumentError("cert: m ≤ 1: $(cert.m)"))

    # check basic modulus properties
    (N & 3) == 3 ||
        throw(ArgumentError("cert: N ≠ 3 mod 4 (N=$N)"))
    isprime(N >> 1) ||
        throw(ArgumentError("cert: (N-1)/2 not prime (N=$N)"))
    powermod(2, N >> 1, N) ∉ (1, N-1) ||
        throw(ArgumentError("cert: 2 fails to witness compositeness (N=$N)"))

    # check that cert contains enough roots
    length(cert.roots) ≥ ROOT_SAMPLES ||
        throw(ArgumentError("cert: too few Nth roots (N=$N)"))
    length(cert.sqrts) ≥ SQRT_MINIMUM ||
        throw(ArgumentError("cert: too few sqrts (N=$N)"))

    # check N not divisible by odd p ≤ TRIALDIV_MAX
    for p = 3:2:TRIALDIV_MAX
        N % p ≠ 0 && continue
        isprime(N ÷ p) && break # N semiprime w. very small factor
        throw(ArgumentError("cert: not semiprime N = $p*$(N ÷ p)"))
    end

    # check provided Nth roots
    for (i, r) in enumerate(cert.roots)
        powermod(r, N, N) == ring_hash(N, :Nth_root, i) ||
            throw(ArgumentError("cert: invalid Nth root (N=$N)"))
    end

    # check provided square roots
    i = 0
    τ = hash_twist(N)
    for (j, r) in enumerate(cert.sqrts)
        j > SQRT_MINIMUM && break # checked enough sqrts
        r² = powermod(r, 2, N)
        while ring_hash(N, :sqrt, i+=1; untwist=τ) ≠ r²
            i ≤ SQRT_SAMPLES ||
                throw(ArgumentError("cert: sqrt test failed (N=$N)"))
        end
    end

    # cert is valid, N is safe
    Client(cert.B, cert.m, cert.N, cert.g)
end

Base.show(io::IO, c::Client) =
    print(io, "Client(B=$(c.B), m=$(c.m), N=$(c.N), x₀=$(c.x₀))")

function hll_generate(client::Client, class::AbstractString)
    N, B, g, x₀ = client.N, client.B, client.g, client.x₀
    h = hash_resource_class(client.salt, class)
    x = mod(x₀ * powermod(g, h, N), N)    # x = x₀ g^h
    t = 2B*rand(rng, 0:(N-1)÷(2B)-1) + 1  # t ∈ ℤ_N st t = 1 mod 2B
    y = powermod(x, t, N)                 # y = x^t
end
