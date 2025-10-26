using Test
using Primes
using AnonymousUserEstimation
using AnonymousUserEstimation:
    gen_prime_pair, jacobi, modulus, factors, lambda, find_g, find_x

const 𝟚 = BigInt(2)

@testset "Prime pair generation" begin
    @testset "correct usage" begin
        for (lo, hi, s) in Any[
                (10, 100, 6)
                (0xf, 0xffff, 22)
                (𝟚^(256-1), 𝟚^256-1, 2^32*nextprime(2^12))
            ]
            P, p = gen_prime_pair(lo, hi, s)
            @test lo ≤ P ≤ hi
            @test isprime(P)
            @test isprime(p)
            @test P == s*p + 1
        end
    end
    @testset "error handling" begin
        @test_throws ArgumentError gen_prime_pair(10, 100, 7)
        @test_throws ArgumentError gen_prime_pair(100, 10, 6)
        @test_throws ArgumentError gen_prime_pair(8, 10, 12)
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
        ring = Ring(2^12+1, 16, 64)
        @test ring isa Ring{UInt64}
        ring = Ring(2^12+1, 16, 65)
        @test ring isa Ring{UInt128}
        ring = Ring(2^12+1, 16, 128)
        @test ring isa Ring{UInt128}
        ring = Ring(2^12+1, 16, 129)
        @test ring isa Ring{BigInt}
    end
    @testset "structure" for _ = 1:10
        ring = Ring(2^4+1, 4, 20)
        N, λ = ring.N, ring.λ
        J⁰ = [x for x in 0:N-1 if gcd(x, N) ≠ 1]
        g = find_g(ring)
        @test jacobi(g, N) == +1
        J⁺ = sort!([powermod(g, k, N) for k in 0:λ-1])
        @test allunique(J⁺)
        x = find_x(ring)
        @test jacobi(x, N) == -1
        J⁻ = sort!(mod.(x .* J⁺, N))
        @test allunique(J⁻)
        # test disjointness
        @test isempty(J⁰ ∩ J⁺)
        @test isempty(J⁰ ∩ J⁻)
        @test isempty(J⁺ ∩ J⁻)
        # test full coverage
        @test isempty(symdiff(J⁰ ∪ J⁻ ∪ J⁺, 0:N-1))
    end
end
