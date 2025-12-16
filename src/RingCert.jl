const L = 2^10
const ε = exp2(-40)
const ROOT_SAMPLES = ceil(Int, -log2(ε)/(log2(L)-log2(3)+1))
const SQRT_SAMPLES = ceil(Int, -log2(ε)/(log2(8)-log2(5)))

struct RingCert{T<:Integer}
    # general shape
    B :: Int # bucket factor (odd)
    m :: Int # max geometric sample size

    # ring-specific info
    N :: T # modulus
    g :: T # semigenerator

    # roots of hash-generated elements
    roots :: Vector{T} # Nth roots
    sqrts :: Vector{T} # square roots
end

function RingCert(ring::Ring{T}) where {T<:Integer}
    B = ring.B
    P, Q = factors(ring)
    N = P*Q

    # test modulus
    mod4(N) == 3 ||
        throw(ArgumentError("modulus: N ≠ 3 mod 4 (N=$N)"))
    gcd(B, N-1) == 1 ||
        throw(ArgumentError("modulus: gcd(B, N-1) ≠ 1 (N=$N)"))

    # generate a semigenerator element
    g = rand_semigenerator(ring)

    # compute Nth roots of hashed elements
    roots = T[]
    D = invmod(N, lambda(ring)) # ∀ x: (x^N)^D == x mod N
    for i = 1:ROOT_SAMPLES
        x = ring_hash(N, :root, i)
        r = powermod(x, D, N) # Nth root of x
        @assert powermod(r, N, N) == x
        push!(roots, r)
    end

    # Bézout & CRT coefficients
    _, u, v = gcdx(P, Q)
    uP = mod(widemul(u, P), N)
    vQ = mod(widemul(v, Q), N)

    # closure for square roots modulo N
    function push_sqrt_mod_N(x::Integer)
        r_P = modsqrt(x, P); r_P === nothing && return false
        r_Q = modsqrt(x, Q); r_Q === nothing && return false
        r = mod(r_P*vQ + r_Q*uP, N) # sqrt of x
        @assert powermod(r, 2, N) == x
        push!(sqrts, r)
        return true
    end

    # compute sqrts of hashed elements
    sqrts = T[]
    τ = hash_twist(N)
    for i = 1:SQRT_SAMPLES
        x = ring_hash(N, :sqrt_x, i; untwist=τ)
        push_sqrt_mod_N(x) && continue
        y = ring_hash(N, :sqrt_y, i; untwist=τ)
        push_sqrt_mod_N(y) && continue
        z = modmul(x, y, N)
        push_sqrt_mod_N(z) && continue
        throw(ArgumentError("ring: fails semiprimality test (N=$N)"))
    end

    return RingCert(ring.B, ring.m, N, g, roots, sqrts)
end

Base.show(io::IO, cert::RingCert) =
    print(io, "RingCert(B=$(cert.B), m=$(cert.m), N=$(cert.N))")

function hash_twist(N::Integer)
    i = 0
    while true
        τ = ring_hash(N, :twist, i += 1)
        jacobi(τ, N) == -1 && return τ
    end
end

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
        if modmul(r, r, p) ≠ x
            r = modmul(r, powermod(2, (p - 1) >> 2, p), p)
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
        z = mod(z + 1, p)
    end

    c = powermod(z, q, p)
    r = powermod(x, (q + 1) >> 1, p)
    t = powermod(x, q, p)
    m = s

    while t ≠ 1
        # find minimal i in [1, m-1] with t^(2^i) == 1
        i = 1
        t′ = modmul(t, t,  p)
        while t′ ≠ 1
            i += 1
            t′ = modmul(t′, t′,  p)
            i ≥ m && break
        end

        # b = c^(2^(m-i-1))
        b = c
        for _ in 1:m-i-1
            b = modmul(b, b,  p)
        end

        r = modmul(r, b, p)
        b′ = modmul(b, b,  p)
        t = modmul(t, b′, p)
        c = b′
        m = i
    end

    return min(r, p - r)
end

modmul(a::Integer, b::Integer, m::Integer) =
    oftype(m, mod(widemul(a, b), m))
