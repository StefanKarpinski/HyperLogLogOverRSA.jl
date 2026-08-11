using Test
using Primes
using Random
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
    # B ≡ 2 mod 4 makes P ≡ 5 mod 8 and so N ≡ 5 mod 8, which puts -1 in J_N^+:
    # a forger gets no public element of J_N^- with known logarithms
    @test mod(ring.B, 4) == 2
    @test mod(ring.P, 8) == 5
    @test mod(ring.Q, 8) == 1
    @test mod(ring.N, 8) == 5
    @test jacobi(ring.N - 1, ring.N) == +1
end
