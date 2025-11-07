using Test
using Primes
using HyperLogLogOverRSA
using HyperLogLogOverRSA:
    gen_prime_pair, jacobi, modulus, factors, lambda, modsqrt,
    rand_semigenerator, rand_jacobi_twist, bucket_map

function check_ring(ring::Ring)
    @test isprime(ring.P)
    @test isprime(ring.Q)
    @test isprime(ring.p)
    @test isprime(ring.q)
    @test ring.N == ring.P*ring.Q
    @test ring.P == 2*ring.B*ring.p + 1
    @test ring.Q == big(2)^ring.m*ring.q + 1
end
