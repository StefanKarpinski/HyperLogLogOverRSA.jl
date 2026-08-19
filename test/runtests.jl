include("setup.jl")

@testset "Prime pairs" begin
    @testset "correct usage" begin
        for (lo, hi, s) in Any[
                (10, 100, 6)
                (0xf, 0xffff, 22)
                (big(2)^(256-1), big(2)^256-1, 2^32*nextprime(2^12))
            ]
            P, p = gen_prime_pair(s, lo, hi)
            @test lo ≤ P ≤ hi
            @test isprime(P)
            @test isprime(p)
            @test P == s*p + 1
        end
    end
    @testset "error handling" begin
        @test_throws ArgumentError gen_prime_pair(7, 10, 100)
        @test_throws ArgumentError gen_prime_pair(6, 100, 10)
        @test_throws ArgumentError gen_prime_pair(12, 8, 10)
    end
end

@testset "Jacobi symbol" begin
    J = Dict{Tuple{Int,Int},Int}()
    primes = [3, 5, 17, 19]
    for P in primes
        squares = sort!(unique(mod(x^2, P) for x in 1:P-1))
        for x in 0:P-1
            j = x == 0 ? 0 : x in squares ? 1 : -1
            @test j == jacobi(x,P)
            J[P,x] = j
        end
    end
    for P in primes, Q in primes
        P < Q || continue
        N = P*Q
        for x in 1:N-1
            j = J[P,mod(x,P)]*J[Q,mod(x,Q)]
            @test j == jacobi(x,N)
        end
    end
    for P in primes, Q in primes, R in primes
        P < Q < R || continue
        N = P*Q*R
        for x in 1:N-1
            j = J[P,mod(x,P)]*J[Q,mod(x,Q)]*J[R,mod(x,R)]
            @test j == jacobi(x,N)
        end
    end
end

@testset "modsqrt" begin
    for p in primes(1000), x = 0:p-1
        R = [r for r = 0:p-1 if mod(r^2, p) == x]
        r = modsqrt(x, p)
        if isempty(R)
            @test r === nothing
        else
            @test r == minimum(R)
            @test mod(r^2, p) == x
        end
    end
end

@testset "hash_into_ring type-independence" begin
    # The folded element must depend only on the value of N, not its in-memory
    # integer type — otherwise a TOML reparse (which picks the type by magnitude)
    # would rederive different certificate challenges and reject a valid cert.
    for Nval in (2663261, Int128(2)^126 + 15, Int128(2)^100 - 25)
        types = Nval ≤ typemax(Int64) ? (Int64, Int128, BigInt) : (Int128, BigInt)
        for keys in ((:sqrt_x, 1), (:sqrt_y, 7))
            ref = hash_into_ring(big(Nval), keys...)
            for T in types
                x = hash_into_ring(T(Nval), keys...)
                @test x == ref              # same value regardless of type
                @test typeof(x) == T        # returned in N's type, in [0, N)
                @test 0 ≤ x < Nval
            end
        end
    end
end

@testset "hash_into_J₊" begin
    # Maps hashes into J_N^+ (as certificate challenges require) and is likewise
    # type-independent. Uses a prime ≡ 5 mod 8: then 2 ∈ J_N^- (so the ×2 twist
    # flips a negative Jacobi symbol into J_N^+) and there are no non-units.
    p = 9223372036854775549       # an Int64 prime ≡ 5 mod 8
    @test isprime(p) && p % 8 == 5 && jacobi(2, p) == -1
    for keys in ((:sqrt_x, 1), (:sqrt_y, 7), (:sqrt_x, 42))
        ref = hash_into_J₊(big(p), keys...)
        @test jacobi(ref, big(p)) == 1        # lands in J_N^+
        for T in (Int64, Int128, BigInt)
            x = hash_into_J₊(T(p), keys...)
            @test x == ref
            @test typeof(x) == T
        end
    end
    # It is a total function into ℤ_N \ J_N^-: for a modulus N ≡ 5 mod 8 (where
    # 𝒥_N(2) = -1, so the ×2 twist works) the result always has Jacobi symbol
    # ≠ -1, even when the raw hash is a non-unit — no throw, no `nothing`.
    for N in 5:8:301, i in 1:8
        @test jacobi(hash_into_J₊(N, :t, i), N) != -1
    end
    @test jacobi(hash_into_J₊(21, :t, 1), 21) == 0  # raw hash is a non-unit, returned as-is
end

@testset "cert rejection" begin
    # A client rejects a certificate whose square roots don't verify — the
    # client-side check is where an unsound (e.g. >2-prime) modulus is caught.
    ring = Ring(2^5+1, 8, 63); cert = RingCert(ring)
    @test Client(cert) isa Client                  # the honest cert is accepted
    bad = RingCert(cert.B, cert.m, cert.N, [cert.sqrts[1] + 1; cert.sqrts[2:end]])
    @test_throws ArgumentError Client(bad)         # a corrupted square root is rejected
end

@testset "hash_into_ring uniformity margin" begin
    # Enough bits are folded beyond the bit-length of N that the modular bias of
    # the reduction is ≤ 2^-HASH_MARGIN. The worst cases are bit-lengths just
    # under a 512-bit block boundary (b ≡ 511 mod 512, e.g. 2047), where the
    # earlier "just cover N" scheme left only ~1 bit of headroom.
    for b in (8, 62, 511, 512, 1023, 2046, 2047, 2048, 4095, 4096)
        N = big(2)^(b - 1) + 1                 # a b-bit odd modulus
        @test Base.top_set_bit(N) == b
        H = 512 * hash_blocks(N)               # bits generated before reduction
        @test H - b ≥ HASH_MARGIN              # ≥128 bits of uniformity
    end
end

# false &&
@testset "Ring sructure" begin
    @testset "basics" begin
        for bits in [55, 63, 64]
            ring = Ring{UInt64}(2^5+1, 8, bits)
            check_ring(ring)
            @test leading_zeros(modulus(ring)) == 64-bits
        end
        ring = Ring(2^5+1, 8, 63)
        @test ring isa Ring{Int64}
        check_ring(ring)
        ring = Ring(2^5+1, 8, 64)
        @test ring isa Ring{Int128}
        check_ring(ring)
        ring = Ring(2^5+1, 8, 127)
        @test ring isa Ring{Int128}
        check_ring(ring)
        ring = Ring(2^5+1, 8, 128)
        @test ring isa Ring{BigInt}
        check_ring(ring)
    end
    # generate some small rings for comprehensive structural testing. These
    # exercise the group structure and decoding, which are independent of the
    # semisharding f, so they are built with `certifiable = false`: at L=20 some
    # (B,m) specs have no shardable modulus (which the default `Ring` rightly
    # rejects), but any modulus of the right shape works for structure.
    rings = Ring{Int}[]
    for B = (3, 5, 9, 11), m = 3:5
        ring = Ring(B, m, 20; certifiable = false)
        check_ring(ring)
        push!(rings, ring)
    end
    @testset "Jacobi classification" for ring in rings
        # jacobi(x) ==  0 <=> not invertible
        # jacobi(x) == +1 <=> x = ±g^k for some k
        # jacobi(x) == -1 <=> x = ±x₀*g^k for some k
        N, λ = ring.N, ring.λ
        g = rand_semigenerator(ring)
        x₀ = rand_jacobi_twist(N)
        @test jacobi(g, N) == +1
        @test jacobi(x₀, N) == -1
        J₀ = [x for x in 0:N-1 if gcd(x, N) ≠ 1]
        # ⟨g⟩ has index 2 in J_N^+, not index 1: no semigenerator reaches the
        # order-2 element of C_4, and -1 is exactly what lies outside ⟨g⟩
        G = [powermod(g, k, N) for k in 0:λ-1]
        @test allunique(G)
        @test !(N-1 in G)
        J₊ = sort!([G; mod.(N .- G, N)])
        J₋ = sort!(mod.(x₀ .* J₊, N))
        @test all(jacobi(x, N) ==  0 for x in J₀)
        @test all(jacobi(x, N) == +1 for x in J₊)
        @test all(jacobi(x, N) == -1 for x in J₋)
        J = [J₀; J₊; J₋]
        @test allunique(J)
        @test length(J) == N
        @test all(0 ≤ x < N for x in J)
    end
    @testset "HyperLogLog frequencies" for ring in rings
        B, m, N, = ring.B, ring.m, ring.N
        pq = ring.p * ring.q
        counts = fill(0, B, m)
        bmap = bucket_map(ring)
        for x in 0:N-1
            jacobi(x, N) == -1 || continue
            b, k = hll_decode(ring, x; bmap)
            counts[b+1,k+1] += 1
        end
        # an exact B × m rectangle: every cell one orbit, with the two top rungs
        # of the ladder merged into the capped level m-1
        @test counts == [
            pq << max(2, m-k)
            for b = 0:B-1, k = 0:m-1
        ]
        @test sum(counts) == count(x -> jacobi(x, N) == -1, 0:N-1)
    end
end

# false &&
@testset "HLL gen & decode" begin
    specs = [
        Int64  => (2^5+1, 8, 63)
        Int128 => (2^9+1, 16, 127)
        BigInt => (2^12-1, 32, 512)
    ]
    for (T, (B, m, L)) in specs
        ring = Ring(B, m, L); cert = RingCert(ring)
        check_ring(ring)
        @test ring isa Ring{T}
        @test cert isa RingCert{T}
        client = Client(cert)
        @test client isa Client{T}
        bmap = bucket_map(ring)
        for uuid = 1:100
            Y = [hll_generate(client, "/package/$uuid") for _ = 1:100]
            H = [hll_decode(ring, y; bmap) for y in Y]
            @test allunique(Y)  # tokens are unlinkable
            @test allequal(H)   # but decode to one stable (b,k) per class
        end
    end
end

# false &&
@testset "semisharding" begin
    # The client's bucket is its own b₀ in every class, while the rank shards.
    ring = Ring(2^12-1, 8, 63; rng = Xoshiro(1)); cert = RingCert(ring)
    bmap = bucket_map(ring)
    client = Client(cert; rng = Xoshiro(2))
    vals = [hll_decode(ring, hll_generate(client, "/c/$i"); bmap) for i = 1:200]
    @test allequal(first.(vals))  # bucket b₀ fixed across all classes
    @test !allequal(last.(vals))  # rank shards across classes
    # distinct clients spread across buckets — this is what lets the server count
    buckets = [hll_decode(ring, hll_generate(Client(cert)), bmap=bmap)[1] for _ = 1:64]
    @test length(unique(buckets)) > 1
end

# false &&
@testset "HLL estimate" begin
    # End-to-end cardinality recovery: n distinct clients each send a few requests
    # in ONE resource class, and the estimate recovers n unique clients. Under
    # semisharding a client shares its bucket across classes, so it is distinct
    # clients — not distinct classes — that populate a per-class sketch. Seeded
    # for determinism; HLL relative error is ≈ 1.04/√B ≈ 1.6%.
    rng = Xoshiro(0)
    ring = Ring(2^12-1, 16, 63; rng); cert = RingCert(ring)
    check_ring(ring)
    n = 5000
    clients = [Client(cert; rng) for _ = 1:n]
    Y = [hll_generate(c, "/registries"; rng) for c in clients for _ = 1:rand(rng, 1:3)]
    @test allunique(Y)  # every emitted token is freshly randomized
    n̂ = hll_estimate(ring, Y)
    @test abs(n̂ - n)/n ≤ 0.05  # ≈ 3·RSE
    @test n̂ ≤ length(Y)  # never more unique clients than requests
end

@testset "HLL estimate request-count cap" begin
    # The batch inflation/decode oracle: a malicious client sends x₀·g^h with
    # h ≡ 0 (mod 2^m) — pinning the geometric coordinate at c₀, so every token
    # has the same rank k — while h = 2^m·j sweeps all B buckets (gcd(2^m,B)=1).
    # The resulting flat sketch estimates to ≈ (1/ln2)·B·2^k ≥ 1.44·B from only B
    # requests. Capping the report at the request count clips that back to B, so
    # an attacker cannot inflate the count above the number of requests it sends.
    rng = Xoshiro(1)
    ring = Ring(2^7-1, 8, 128; rng)
    B, m, N = ring.B, ring.m, ring.N
    g  = rand_semigenerator(ring; rng)
    x₀ = rand_jacobi_twist(N; rng)
    batch = [modmul(x₀, powermod(g, oftype(N, 2)^m * j, N), N) for j = 0:B-1]
    bmap = bucket_map(ring)
    sketch = [hll_decode(ring, y; bmap) for y in batch]
    @test length(unique(first.(sketch))) == B   # every bucket covered exactly once
    @test allequal(last.(sketch))               # all one rank ⇒ flat sketch
    n̂ = hll_estimate(ring, batch)
    @test n̂ == B                                # capped at the request count …
    @test α_∞ * B * exp2(last(sketch[1]) + 1) > B  # … below the uncapped estimate
end
