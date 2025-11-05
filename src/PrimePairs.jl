## Generating prime pairs (P, p) such that:
#
#   P == scale*p + 1 <=> p == (P - 1)/scale
#

function next_paired_prime(scale::Integer, P::Integer, P_max::Integer)
    while P ≤ P_max
        P = scale*nextprime(cld(P - 1, scale)) + 1
        isprime(P) && break
        P = nextprime(P)
    end
    if P ≤ P_max
        @assert isprime(P)
        @assert isprime((P-1) ÷ scale)
    end
    return P
end

function prev_paired_prime(scale::Integer, P_min::Integer, P::Integer)
    while P_min ≤ P
        P = scale*prevprime(fld(P - 1, scale)) + 1
        isprime(P) && break
        P = prevprime(P)
    end
    if P_min ≤ P
        @assert isprime(P)
        @assert isprime((P-1) ÷ scale)
    end
    return P
end

function narrow_prime_range(scale::Integer, P_min::Integer, P_max::Integer)
    iseven(scale) || throw(ArgumentError("scale factor must be even"))
    P_max = prev_paired_prime(scale, P_min, P_max)
    P_min = next_paired_prime(scale, P_min, P_max)
    return P_min, P_max
end

"""
    gen_prime_pair(
        P_min :: Integer,
        P_max :: Integer,
        scale :: Integer,
    ) -> P, p

Return a pair of primes, `(P, p)`, such that:

- `P in P_min:P_max`
- `P == scale*p + 1`

Requires that `scale` is even.
"""
function gen_prime_pair(
    scale :: Integer,
    P_min :: Integer,
    P_max :: Integer,
)
    iseven(scale) || throw(ArgumentError("scale factor must be even"))

    p_min = cld(P_min - 1, scale)
    p_max = fld(P_max - 1, scale)
    p_min += iseven(p_min)
    p_max -= iseven(p_max)
    p_range = p_min:p_max

    while true
        # sample an odd value
        p = rand(rng, p_range) | 1
        # check if it's good
        isprime(p) || continue
        P = scale*p + 1
        isprime(P) || continue
        return P, p
    end
end
