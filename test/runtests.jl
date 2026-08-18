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
    # The certificate challenge must depend only on the value of N, not its
    # in-memory integer type — otherwise a TOML reparse (which picks the type by
    # magnitude) would rederive different challenges and reject a valid cert.
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
        # the twist path must be type-independent too
        τ = big(2)
        @test hash_into_ring(big(Nval), :sqrt_x, 3; untwist=τ) ==
              hash_into_ring(Int128(Nval), :sqrt_x, 3; untwist=Int128(2))
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
    # generate some small rings for comprehensive testing
    rings = Ring{Int}[]
    for B = (3, 5, 9, 11), m = 3:5
        ring = Ring(B, m, 20)
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
    rings = [
        Int64  => Ring(2^5+1, 8, 63)
        Int128 => Ring(2^9+1, 16, 127)
        BigInt => Ring(2^12-1, 32, 512)
    ]
    for (T, ring) in rings
        check_ring(ring)
        @test ring isa Ring{T}
        cert = RingCert(ring)
        @test cert isa RingCert{T}
        client = Client(cert)
        @test client isa Client{T}
        bmap = bucket_map(ring)
        for uuid = 1:100
            Y = [hll_generate(client, "/package/$uuid") for _ = 1:100]
            H = [hll_decode(ring, y; bmap) for y in Y]
            @test allunique(Y)
            @test allequal(H)
        end
    end
end

# false &&
@testset "HLL estimate" begin
    # End-to-end cardinality recovery through the whole protocol, with a seeded
    # rng for determinism (most seeds pass; pick a new one if the PRNG changes).
    # The rand(1:3) repeats collapse — same client+class decodes to one value —
    # so the estimate counts the n distinct classes, not requests. HLL's relative
    # error is ≈ 1.04/√B ≈ 1.6%; this seed lands at ~1.2%.
    rng = Xoshiro(0)
    ring = Ring(2^12-1, 16, 63; rng)
    check_ring(ring)
    client = Client(RingCert(ring; rng); rng)
    n = 5000
    Y = [hll_generate(client, "/package/$id"; rng) for id = 1:n for _ = 1:rand(rng, 1:3)]
    @test allunique(Y)                  # every emitted token is freshly randomized
    n̂ = hll_estimate(ring, Y)
    @test abs(n̂ - n)/n ≤ 0.05           # ≈ 3·RSE
end
