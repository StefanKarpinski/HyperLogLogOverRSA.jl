## binary/reciprocity-based Jacobi

"""
    jacobi(x::Integer, N::Integer) -> Int

Compute the Jacobi symbol (x|N) for odd positive `N`.
Returns `-1`, `0`, or `1`.
"""
function jacobi(x::Integer, N::Integer)::Int
    # domain checks
    N > 0 || throw(ArgumentError("N must be positive"))
    isodd(N) || throw(ArgumentError("N must be odd"))

    # reduce x modulo N; handle zero quickly
    x = mod(x, N)
    iszero(x) && return 0

    # main loop
    s = 1
    while !iszero(x)
        # strip powers of two from x: (2|N) = -1 if N ≡ 3 or 5 (mod 8)
        z = trailing_zeros(x)
        if !iszero(z)
            if isodd(z)
                Nm8 = mod8(N)
                if Nm8 == 3 || Nm8 == 5
                    s = -s
                end
            end
            x >>= z
        end

        # quadratic reciprocity: flip sign if x ≡ N ≡ 3 (mod 4)
        if mod4(x) == 3 && mod4(N) == 3
            s = -s
        end

        # reduce and swap
        x, N = mod(N, x), x

        # early out: if x == 0, symbol is 0 unless N == 1
        iszero(x) && return N == 1 ? s : 0
    end

    # if we get here, gcd(original x, N) ≠ 1 unless N == 1
    return N == 1 ? s : 0
end

## Generate a random Jacobi "twist" element

"""
    rand_jacobi_twist(N::Integer)

Generate random τ ∈ ℤ_N with jacobi(τ) == -1.
"""
function rand_jacobi_twist(N::Integer)
    range = 0:N-1
    while true
        τ = rand(rng, range)
        jacobi(τ, N) == -1 && return τ
    end
end

## Small inline helpers with BigInt-specialized fast paths

@inline mod4(x::Integer) = Int(x & (4-1))
@inline mod8(x::Integer) = Int(x & (8-1))

@inline mod4(x::BigInt) = Int(
    Base.GMP.MPZ.tstbit(x, 0) +
    Base.GMP.MPZ.tstbit(x, 1) << 1
)
@inline mod8(x::BigInt) = Int(
    Base.GMP.MPZ.tstbit(x, 0) +
    Base.GMP.MPZ.tstbit(x, 1) << 1 +
    Base.GMP.MPZ.tstbit(x, 2) << 2
)
