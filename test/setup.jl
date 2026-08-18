using Test
using Primes
using Random
using HyperLogLogOverRSA
using HyperLogLogOverRSA:
    gen_prime_pair, jacobi, modulus, factors, lambda, modsqrt,
    rand_semigenerator, rand_jacobi_twist, bucket_map,
    hash_into_ring, hash_into_J₊, hash_blocks, HASH_MARGIN, modmul, α_∞,
    derive_f, f_shards

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

# Generate a ring of the given size together with its certificate, retrying ring
# generation when the canonical semisharding generator f happens not to shard —
# about half of moduli, which the server rejects and regenerates.
function valid_ring_cert(B, m, L; rng = Random.default_rng())
    while true
        ring = Ring(B, m, L; rng)
        try
            return ring, RingCert(ring)
        catch e
            (e isa ArgumentError && occursin("does not shard", e.msg)) || rethrow(e)
        end
    end
end
