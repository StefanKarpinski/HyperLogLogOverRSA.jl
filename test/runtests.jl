using Test
using Primes
using HyperLogLogOverRSA
using HyperLogLogOverRSA:
    gen_prime_pair, jacobi, modulus, factors, lambda, find_g, find_x

const 𝟚 = BigInt(2)

@testset "Prime pair generation" begin
    @testset "correct usage" begin
        for (lo, hi, s) in Any[
                (10, 100, 6)
                (0xf, 0xffff, 22)
                (𝟚^(256-1), 𝟚^256-1, 2^32*nextprime(2^12))
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

@testset "Ring" begin
    @testset "basics" begin
        for bits in [55, 63, 64]
            ring = Ring{UInt64}(2^5+1, 8, bits)
            @test leading_zeros(modulus(ring)) == 64-bits
            @test all(isprime, factors(ring))
            @test isprime(ring.p)
            @test isprime(ring.q)
        end
        ring = Ring(2^12+1, 16, 63)
        @test ring isa Ring{Int64}
        ring = Ring(2^12+1, 16, 64)
        @test ring isa Ring{Int128}
        ring = Ring(2^12+1, 16, 127)
        @test ring isa Ring{Int128}
        ring = Ring(2^12+1, 16, 128)
        @test ring isa Ring{BigInt}
    end
    rings = Ring{Int}[]
    for log_B = 2:5, m = 2:5
        B = 2^log_B + 1
        ring = Ring(B, m, 20)
        ring in rings || push!(rings, ring)
    end
    @testset "structure" for ring in rings
        N, λ = ring.N, ring.λ
        g = find_g(ring)
        x = find_x(ring)
        @test jacobi(g, N) == +1
        @test jacobi(x, N) == -1
        # classification of elements by Jacobi symbol:
        #   jacobi(y) ==  0 <=> not invertible
        #   jacobi(y) == +1 <=> g^k for some k
        #   jacobi(y) == -1 <=> xg^k for some k
        J₀ = [x for x in 0:N-1 if gcd(x, N) ≠ 1]
        J₊ = sort!([powermod(g, k, N) for k in 0:λ-1])
        J₋ = sort!(mod.(x .* J₊, N))
        # test full disjoint coverage
        J = [J₀; J₊; J₋]
        @test allunique(J)
        @test length(J) == N
        @test all(0 ≤ y < N for y in J)
        # HyperLogLog statistics

    end
end
