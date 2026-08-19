using Test
using Primes
using Random
using HyperLogLogOverRSA
using HyperLogLogOverRSA:
    gen_prime_pair, jacobi, modulus, factors, lambda, modsqrt,
    rand_semigenerator, rand_jacobi_twist, bucket_map,
    hash_into_ring, hash_into_J₊, hash_blocks, HASH_MARGIN, modmul, α_∞

function check_ring(ring::Ring)
    @test isprime(ring.P)
    @test isprime(ring.Q)
    @test isprime(ring.p)
    @test isprime(ring.q)
    @test ring.N == ring.P*ring.Q
    @test ring.P == 4*ring.B*ring.p + 1
    @test ring.Q == big(2)^ring.m*ring.q + 1
    # P ≡ 5 mod 8 and Q ≡ 1 mod 8 give N ≡ 5 mod 8, hence jacobi(-1, N) = +1: the
    # one element with public logarithms stays out of J_N^-
    @test mod(ring.P, 8) == 5
    @test mod(ring.Q, 8) == 1
    @test mod(ring.N, 8) == 5
    @test jacobi(ring.N - 1, ring.N) == +1
end
