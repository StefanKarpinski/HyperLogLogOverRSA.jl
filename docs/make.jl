using Documenter
using HyperLogLogOverRSA

# LaTeX macros mirroring the writeup's preamble (preamble.sty), registered with
# KaTeX so the math in the writeup pages renders. Keys are the macro names; `#1`
# etc. are macro arguments. Unused entries are harmless.
const MATH_MACROS = Dict(
    raw"\qq"      => raw"\qquad",
    raw"\N"       => raw"\mathbb{N}",
    raw"\Z"       => raw"\mathbb{Z}",
    raw"\Q"       => raw"\mathbb{Q}",
    raw"\R"       => raw"\mathbb{R}",
    raw"\C"       => raw"\mathbb{C}",
    raw"\nfrac"   => raw"#1/#2",
    raw"\half"    => raw"\tfrac{1}{2}",
    raw"\minmax"  => raw"\operatorname{minmax}",
    raw"\eps"     => raw"\operatorname{eps}",
    raw"\fr"      => raw"\operatorname{frac}",
    raw"\fmod"    => raw"\operatorname{mod}",
    raw"\lcm"     => raw"\operatorname{lcm}",
    raw"\invmod"  => raw"\operatorname{invmod}",
    raw"\sig"     => raw"\operatorname{sig}",
    raw"\bitlen"  => raw"\operatorname{bitlen}",
    raw"\fld"     => raw"\operatorname{fld}",
    raw"\cld"     => raw"\operatorname{cld}",
    raw"\fldmod"  => raw"\operatorname{fldmod}",
    raw"\cldmod"  => raw"\operatorname{cldmod}",
    raw"\divides" => raw"\mid",
    raw"\inter"   => raw"\,\cap\,",
    raw"\union"   => raw"\,\cup\,",
    raw"\and"     => raw"\,\land\,",
    raw"\or"      => raw"\,\lor\,",
    raw"\xor"     => raw"\,\veebar\,",
    # NB: do NOT override \not — KaTeX builds \neq, \notin, and the unicode
    # ≠ / ∉ on top of the built-in \not, so redefining it renders them as
    # "¬ =" / "¬ ∈". The writeup uses ≠ and \notin and never a standalone \not
    # for logical negation (use \lnot / \neg for that).
    raw"\E"       => raw"\exists\:",
    raw"\A"       => raw"\forall\:",
    raw"\st"      => raw"\:\middle|\:",
    raw"\bra"     => raw"\left[#1\right]",
    raw"\set"     => raw"\left\{#1\right\}",
    raw"\seq"     => raw"\left\langle #1 \right\rangle",
    raw"\gen"     => raw"\left\langle #1 \right\rangle",
    raw"\norm"    => raw"\left\lvert #1 \right\rvert",
    raw"\Norm"    => raw"\left\lVert #1 \right\rVert",
    raw"\ceil"    => raw"\left\lceil #1 \right\rceil",
    raw"\floor"   => raw"\left\lfloor #1 \right\rfloor",
    raw"\round"   => raw"\operatorname{round}(#1)",
    raw"\float"   => raw"\operatorname{float}",
    raw"\tz"      => raw"\operatorname{tz}",
    raw"\lz"      => raw"\operatorname{lz}",
    raw"\ord"     => raw"\operatorname{ord}",
    raw"\Im"      => raw"\operatorname{Im}",
    raw"\Ker"     => raw"\operatorname{Ker}",
    raw"\img"     => raw"\operatorname{img}",
    raw"\hll"     => raw"\mathcal{H}",
    raw"\Jacobi"  => raw"\mathcal{J}",
    raw"\hash"    => raw"\mathsf{H}",
)

makedocs(
    sitename = "HyperLogLog Over RSA",
    authors  = "Stefan Karpinski",
    modules  = [HyperLogLogOverRSA],
    checkdocs = :exports,  # only require docstrings for exported names
    format = Documenter.HTML(
        mathengine = Documenter.KaTeX(Dict(:macros => MATH_MACROS)),
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical  = "https://karpinski.org/HyperLogLogOverRSA.jl",
        assets     = ["assets/custom.css"],
    ),
    pages = [
        "Home"               => "index.md",
        "Executive Summary"  => "01-executive-summary.md",
        "For Cryptographers" => "02-for-cryptographers.md",
        "The Writeup" => [
            "Anonymously Counting Users" => "03-counting-users.md",
            "HyperLogLog Over RSA"       => "04-hyperloglog-over-rsa.md",
            "Security Analysis (Proofs)" => "05-security-analysis.md",
            "Protocol Summary"           => "06-protocol-summary.md",
        ],
        "Reference Implementation" => "reference.md",
    ],
)

deploydocs(
    repo      = "github.com/StefanKarpinski/HyperLogLogOverRSA.jl",
    devbranch = "main",
)
