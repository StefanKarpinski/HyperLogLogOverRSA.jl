using Distributions

const ε = exp2(-30)
const SQRT_SAMPLES = ceil(Int, -21.96*log2(ε))
const SQRT_MINIMUM = quantile(Binomial(SQRT_SAMPLES, 0.5), ε)

struct RingCert{T<:Integer}
    # general shape
    B :: Int # bucket factor (odd)
    m :: Int # max geometric sample size

    # ring-specific info
    N :: T # modulus
    g :: T # semigenerator

    # square roots of hash-generated elements
    sqrts :: Vector{Pair{Int,T}}
end

function RingCert(ring::Ring{T}) where {T<:Integer}
    P, Q = factors(ring)
    N = P*Q

    # test required properties
    (N & 3) == 3 ||
        throw(ArgumentError("modulus: N ≠ 3 mod 4 (N=$N)"))
    isprime(N >> 1) ||
        throw(ArgumentError("modulus: (N-1)/2 not prime (N=$N)"))
    powermod(2, N >> 1, N) ∉ (1, N-1) ||
        throw(ArgumentError("modulus: 2 fails to witness compositeness (N=$N)"))

    # generate a semigenerator element
    g = rand_semigenerator(ring)

    # Bézout & CRT coefficients
    _, u, v = gcdx(P, Q)
    uP = mod(u*P, N)
    vQ = mod(v*Q, N)

    # compute sqrts that exist
    sqrts = Vector{Pair{Int,T}}()
    i = j = 0
    while j < SQRT_SAMPLES
        i += 1
        x = ring_hash(N, i)
        jacobi(x, N) == 1 || continue
        j += 1
        r_P = modsqrt(x, P)
        r_P === nothing && continue
        r_Q = modsqrt(x, Q)
        r_Q === nothing && continue
        r = mod(r_P*vQ + r_Q*uP, N)
        @assert mod(r^2, N) == x
        push!(sqrts, i => r)
    end
    length(sqrts) ≥ SQRT_MINIMUM ||
        throw(ArgumentError("ring: fails semiprimality test (N=$N)"))

    return RingCert(ring.B, ring.m, N, g, sqrts)
end

Base.show(io::IO, cert::RingCert) =
    print(io, "RingCert(B=$(cert.B), m=$(cert.m), N=$(cert.N))")

"""
    modsqrt(x::Integer, p::Integer) -> Union{Integer,Nothing}

Return r in 0..p-1 such that r^2 ≡ x (mod p), where p is prime.
If no square root exists, return `nothing`.

Fast paths:
- p == 2
- p % 4 == 3: r = x^((p+1)/4) mod p
- p % 8 == 5: Atkin’s shortcut

Otherwise uses Tonelli–Shanks (always works for odd primes).

The returned root is canonical: r ≤ p - r.
"""
function modsqrt(x::Integer, p::Integer)
    p ≥ 2 && isprime(p) ||
        throw(ArgumentError("Modulus must be prime"))

    x = mod(x, p)
    p == 2 && return x # only 0 or 1
    x == 0 && return zero(p) # sqrt(0) = 0

    # Legendre symbol via Euler’s criterion
    l = powermod(x, (p - 1) >> 1, p)
    l == p - 1 && return nothing # ≡ -1 mod p
    @assert isone(l) # guaranteed for prime p

    # p ≡ 3 (mod 4): one-shot
    if (p & 3) == 3
        r = powermod(x, (p + 1) >> 2, p)
        return min(r, p - r)
    end

    # p ≡ 5 (mod 8): Atkin's method
    if (p & 7) == 5
        r = powermod(x, (p + 3) >> 3, p)
        # if wrong square, multiply by 2^((p-1)/4)
        if r^2 % p ≠ x
            r = mod(r*powermod(2, (p - 1) >> 2, p), p)
        end
        return min(r, p - r)
    end

    # Tonelli–Shanks for the general case
    # write p-1 = q * 2^s with q odd
    q = p - 1
    s = trailing_zeros(q)
    q >>= s

    # find a quadratic non-residue z (Legendre = -1)
    z = 2
    while powermod(z, (p - 1) >> 1, p) == 1
        z = (z + 1) % p
    end

    c = powermod(z, q, p)
    r = powermod(x, (q + 1) >> 1, p)
    t = powermod(x, q, p)
    m = s

    while t ≠ 1
        # find minimal i in [1, m-1] with t^(2^i) == 1
        i = 1
        t′ = t^2 % p
        while t′ ≠ 1
            i += 1
            t′ = t′^2 % p
            i ≥ m && break
        end

        # b = c^(2^(m-i-1))
        b = c
        for _ in 1:m-i-1
            b = mod(b^2, p)
        end

        r = mod(r * b, p)
        b′ = mod(b^2, p)
        t = mod(t * b′, p)
        c = b′
        m = i
    end

    return min(r, p - r)
end
