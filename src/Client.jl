struct Client{T<:Integer}
    B :: Int # bucket factor (odd)
    m :: Int # max geometric sample size
    N :: T   # ring modulus
    g :: T   # common published ring semigenerator
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

    (N & 3) == 3 ||
        throw(ArgumentError("modulus: N ≠ 3 mod 4 (N=$N)"))
    isprime(N >> 1) ||
        throw(ArgumentError("modulus: (N-1)/2 not prime (N=$N)"))
    powermod(2, N >> 1, N) ∉ (1, N-1) ||
        throw(ArgumentError("modulus: 2 fails to witness compositeness (N=$N)"))

    # test N not divisible by any odd p ≤ TRIAL_DIV_MAX
    for p = 3:2:TRIAL_DIV_MAX
        N % p ≠ 0 && continue
        isprime(N ÷ p) && break # N is just a very small semiprime
        throw(ArgumentError("semiprimality: divisible by $p (N=$N)"))
    end

    # check the provided Nth roots
    length(cert.roots) ≥ NTH_ROOT_SAMPLES ||
        throw(ArgumentError("modulus: too few Nth roots (N=$N)"))
    for i = 1:NTH_ROOT_SAMPLES
        r = cert.roots[i]
        x = ring_hash(N, :Nth_root, i)
        powermod(r, N, N) == x ||
            throw(ArgumentError("modulus: invalid Nth root (N=$N)"))
    end

    # find a deterministic twist element
    τ = zero(N)
    for i = 1:N
        τ = ring_hash(N, :twist, i)
        jacobi(τ, N) == -1 && break
    end

    # check the provided square roots
    sqrts = Dict(cert.sqrts)
    i = j = 0
    while i < SQRT_SAMPLES
        i += 1
        x = ring_hash(N, :sqrt, i)
        if jacobi(x, N) == -1
            x = mod(widemul(τ, x), N)
        end
        r = get(sqrts, i, nothing)
        r === nothing && continue
        powermod(r, 2, N) == x ||
            throw(ArgumentError("semiprimality: invalid sqrt (N=$N)"))
        j += 1
    end
    j ≥ SQRT_MINIMUM ||
        throw(ArgumentError("semiprimality: too few sqrts (N=$N)"))

    Client(cert.B, cert.m, cert.N, cert.g)
end

Base.show(io::IO, c::Client) =
    print(io, "Client(B=$(c.B), N=$(c.N)), x₀=$(c.x₀))")

function hll_generate(client::Client, class::AbstractString)
    N, B, g, x₀ = client.N, client.B, client.g, client.x₀
    h = hash_resource_class(client.salt, class)
    x = mod(x₀ * powermod(g, h, N), N)    # x = x₀ g^h
    t = 2B*rand(rng, 0:(N-1)÷(2B)-1) + 1  # t ∈ ℤ_N st t = 1 mod 2B
    y = powermod(x, t, N)                 # y = x^t
end
