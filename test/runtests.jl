using Test
using Primes
using AnonymousUserEstimation
using AnonymousUserEstimation: gen_prime_pair, modulus, factors

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

@testset "Ring" begin
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
