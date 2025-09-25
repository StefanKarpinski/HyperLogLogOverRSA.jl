using AnonymousUserEstimation
using Primes
using Test

const 𝟚 = BigInt(2)

@testset "Prime pair generation" begin
    gen_prime_pair = AnonymousUserEstimation.gen_prime_pair
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

# @testset "Ring" begin
#     ring = Ring{UInt64}()
# end
