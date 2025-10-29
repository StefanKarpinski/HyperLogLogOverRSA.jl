using SHA
using Base: SHA1

struct Client{T<:Integer}
    B :: Int # bucket factor (odd)
    N :: T   # ring modulus
    g :: T   # common published ring semigenerator
    x₀ :: T  # client-specific random Jacobi twist element
    salt :: SHA1 # pre-computed hash of x₀
end

function Client(
    B :: Int,
    N :: T,
    g :: T,
) where {T<:Integer}
    x₀ = rand_jacobi_twist(N)
    salt = SHA1(sha1(string(x₀, base=62)))
    Client(B, N, g, x₀, salt)
end

Client(ring::Ring) = Client(ring.B, ring.N, rand_semigenerator(ring))

Base.show(io::IO, c::Client) =
    print(io, "Client(B=$(c.B), N=$(c.N)), x₀=$(c.x₀))")

function hll_generate(client::Client, class::AbstractString)
    N, B, g, x₀ = client.N, client.B, client.g, client.x₀
    h = hash_resource_class(client.salt, class)
    x = mod(x₀ * powermod(g, h, N), N)    # x = x₀ g^h
    t = 2B*rand(rng, 0:(N-1)÷(2B)-1) + 1  # t ∈ ℤ_N st t = 1 mod 2B
    y = powermod(x, t, N)                 # y = x^t
end

function hash_resource_class(salt::SHA1, class::AbstractString)
    bytes = sha1("$salt\0$class\0")
    h = zero(UInt64)
    for i = 1:8
        h <<= 8
        h |= bytes[i]
    end
    return h
end
