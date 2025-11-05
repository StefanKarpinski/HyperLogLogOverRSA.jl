using Test
using Primes
using HyperLogLogOverRSA
using HyperLogLogOverRSA:
    gen_prime_pair, jacobi, modulus, factors, lambda, modsqrt,
    rand_semigenerator, rand_jacobi_twist, bucket_map

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

# false &&
@testset "Ring sructure" begin
    @testset "basics" begin
        for bits in [55, 63, 64]
            ring = Ring{UInt64}(2^5+1, 8, bits)
            @test leading_zeros(modulus(ring)) == 64-bits
            @test all(isprime, factors(ring))
            @test isprime(ring.p)
            @test isprime(ring.q)
        end
        ring = Ring(2^5+1, 8, 63)
        @test ring isa Ring{Int64}
        ring = Ring(2^5+1, 8, 64)
        @test ring isa Ring{Int128}
        ring = Ring(2^5+1, 8, 127)
        @test ring isa Ring{Int128}
        ring = Ring(2^5+1, 8, 128)
        @test ring isa Ring{BigInt}
    end
    # generate some small rings for comprehensive testing
    rings = Ring{Int}[]
    for log_B = 2:5, m = 2:5
        B = 2^log_B + 1
        ring = Ring(B, m, 20, certifiable=false)
        ring in rings || push!(rings, ring)
    end
    @testset "Jacobi classification" for ring in rings
        # jacobi(x) ==  0 <=> not invertible
        # jacobi(x) == +1 <=> x = g^k for some k
        # jacobi(x) == -1 <=> x = x₀*g^k for some k
        N, λ = ring.N, ring.λ
        g = rand_semigenerator(ring)
        x₀ = rand_jacobi_twist(N)
        @test jacobi(g, N) == +1
        @test jacobi(x₀, N) == -1
        J₀ = [x for x in 0:N-1 if gcd(x, N) ≠ 1]
        J₊ = sort!([powermod(g, k, N) for k in 0:λ-1])
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
        counts = fill(0, B, m+1)
        bmap = bucket_map(ring)
        for x in 0:N-1
            jacobi(x, N) == -1 || continue
            b, k = hll_decode(ring, x; bmap)
            counts[b+1,k+1] += 1
        end
        @test counts == [
            pq << max(0, m-k-1)
            for b = 0:B-1, k = 0:m
        ]
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
