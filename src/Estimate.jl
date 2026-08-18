## Improved HyperLogLog cardinality estimator
#
# Implements the count-estimation step of the protocol (Server step 6 in the
# writeup). Given the encrypted tokens logged for one resource class, the ring
# holder decodes each to its (bucket, geometric) value, aggregates the maximum
# geometric sample per bucket into a HyperLogLog sketch, and estimates the
# unique client count. Repeated tokens from the same client collapse, since they
# all decode to that client's single HLL value for the class.
#
# The estimate is capped at the number of requests aggregated: there cannot be
# more unique clients in a set of requests than there are requests. Besides being
# trivially correct, this cap is what bounds a malicious client's inflation to
# the number of requests it sends — a flat sketch covering every bucket at rank k
# estimates to ≈ (1/ln2)·B·2^k from only B requests, and the cap clips that back
# to B (see the "Malicious clients" section of the security analysis).
#
# The estimator is Otmar Ertl's improved estimator, "New Cardinality Estimation
# Methods for HyperLogLog Sketches" (2017), arXiv:1706.07290: a closed-form
# estimator, effectively unbiased across the whole cardinality range with no
# separate small/large-range corrections. Our decoded geometric sample
# k ∈ {0,…,m-1} is Ertl's register value r = k + 1 (so r = 0 marks an empty
# bucket and r = m a saturated one), making our B buckets a standard HyperLogLog
# sketch with saturation level q = m - 1.

# σ(x) = x + Σ_{k≥1} x^(2^k) 2^(k-1) — correction term for empty registers.
function σ(x::Float64)
    x ≥ 1 && return Inf
    y, p, w = x, x, 1.0
    while true
        p *= p          # p = x^(2^k)
        Δ = p*w         # w = 2^(k-1)
        y + Δ == y && return y
        y += Δ
        w *= 2
    end
end

# τ(x) = (1/3)[(1-x) - Σ_{k≥1} (1 - x^(2^-k))^2 2^-k] — correction for saturated
# registers.
function τ(x::Float64)
    (x ≤ 0 || x ≥ 1) && return 0.0
    y, r, w = 1 - x, x, 0.5
    while true
        r = sqrt(r)         # r = x^(2^-k)
        Δ = (1 - r)^2 * w   # w = 2^-k
        y - Δ == y && return y/3
        y -= Δ
        w *= 0.5
    end
end

const α_∞ = 1/(2log(2))

"""
    hll_estimate(ring::Ring, tokens) -> Float64

Estimate the number of unique clients behind a collection of encrypted HLL
`tokens` — the `y` values logged for a single resource class. Each token is
decoded with `ring` (sharing one [`bucket_map`](@ref) across the batch) and
aggregated into a HyperLogLog sketch; repeated tokens from the same client
collapse automatically, since they all decode to that client's single
`(bucket, geometric)` value for the class.

Uses the improved estimator of Ertl (2017): unbiased across the whole
cardinality range, with relative standard error ≈ `1.04/√B`. The result is
capped at the number of tokens aggregated — there cannot be more unique clients
than requests — which also bounds a malicious client's count inflation to the
number of requests it sends.
"""
function hll_estimate(ring::Ring, tokens::AbstractVector)
    B, m = ring.B, ring.m
    bmap = bucket_map(ring)
    reg = fill(-1, B)                    # per-bucket max geometric; -1 = empty
    for y in tokens
        b, k = hll_decode(ring, y; bmap)
        reg[b+1] = max(reg[b+1], k)
    end
    # histogram over Ertl register values r = k + 1 ∈ {0,…,m}, stored at
    # C[r+1]: C[1] is the empty count (r=0), C[m+1] the saturated count (r=m)
    C = zeros(Int, m+1)
    for k in reg
        C[k+2] += 1
    end
    d = B*σ(C[1]/B)                       # empty registers (r = 0)
    for r = 1:m-1
        d += C[r+1]*exp2(-r)             # normal registers (r = 1…m-1)
    end
    d += B*τ(1 - C[m+1]/B)*exp2(1-m)      # saturated registers (r = m)
    # cap at the request count: unique clients ≤ requests, and this is what
    # bounds inflation from a curated flat-sketch batch (see security analysis)
    return min(α_∞ * B^2 / d, float(length(tokens)))
end
