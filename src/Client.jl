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

    # check that N has ≤ two distinct prime factors (p ≥ 1 - ε)
    sqrts = Dict(cert.sqrts)
    i = j = k = 0
    while j < SQRT_SAMPLES
        i += 1
        x = ring_hash(N, i)
        jacobi(x, N) == 1 || continue
        j += 1
        r = get(sqrts, i, nothing)
        r === nothing && continue
        powermod(r, 2, N) == x ||
            throw(ArgumentError("semiprimality: invalid sqrt (N=$N)"))
        k += 1
    end
    k ≥ SQRT_MINIMUM ||
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

function hash_resource_class(salt::SHA1, class::AbstractString)
    bytes = sha1("$salt\0$class\0")
    h = zero(UInt64)
    for i = 1:8
        h <<= 8
        h |= bytes[i]
    end
    return h
end
