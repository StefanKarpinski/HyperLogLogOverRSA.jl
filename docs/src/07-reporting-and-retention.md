# Reporting, Publishing & Retention

Client count estimates are not of any use unless they get published. HyperLogLog allows estimating the number of unique clients in arbitrary slices of requests, but some care should be taken in what slices are actually reported.

It is tempting to think that, because each decoded record is only a small HyperLogLog value, the decoded logs are harmless and could simply be published as raw request-level records. They should not be. In any class small enough that clients are unlikely to share $(b, k)$ pairs, these decoded pairs uniquely identify clients, in much the same way that being a FreeBSD user effectively uniquely identifies someone (we *know* what you’re up to, Alex). Publishing decoded HyperLogLog pairs, would also give an attacker an oracle for decoding their HLL values, which, at the very least, would allow them to learn which pairs are rare.

What *is* safe to expose is **aggregate** output for large enough slices: per-slice cardinality, rather than the per-request $(b, k)$ values themselves. This is safe so long as a **minimum-cardinality floor** is enforced—report no slice whose estimated cardinality falls below a threshold $T$, folding everything beneath it into a single “other” bucket. How large $T$ must be has a clean answer: a published estimate keeps a client anonymous exactly when it cannot reveal whether that client is present, and we can measure how close it comes. Adding or removing one client shifts the expected estimate by $1$, against a standard error of about $1.04\,n/\sqrt{B}$, so a single client’s detectability is the ratio $\sqrt{B}/(1.04\,n)$. Once that ratio drops well below $1$, the estimator’s own sampling noise drowns out any individual contribution. Setting $n = \sqrt{B}$ ($\approx 64$ for $B = 2^{12}$) makes that ratio about $1$, so at a cardinality of $\sqrt{B}$ a lone client is still fully visible; bringing it down to a quarter or an eighth takes only a small multiple more, a floor of four to eight times $\sqrt{B}$, a few hundred clients, with a round $T = 1024$ at the conservative end. The bound is for a single published estimate, though; an adversary who differences overlapping slices or repeats queries averages the noise away, and must be limited separately.

These rules also close the inflation channel discussed in [Malicious clients](05-security-analysis.md#Malicious-clients): an attacker who can read fine-grained per-slice counts gains a decoding oracle for their own forged tokens, and denying counts on sub-floor slices removes it. So the cardinality floor pays for itself three times over—anonymity in small populations, safety of published output, and resistance to inflation attacks.

## Three levels of anonymity

It helps to separate three stages at which anonymity matters, because they are protected differently and are easy to conflate:

1. **Raw decoded logs** — individual $(\text{class}, b, k)$ records, freshly decoded, still bearing whatever coarse metadata (OS, version, region, coarse time) the request carried. This is the most sensitive form.
2. **Stored aggregates** — what is retained after decoding: ideally not the raw records at all, but pre-aggregated HyperLogLog sketches, per class and per slice.
3. **Published output** — the per-slice cardinalities actually released, behind the minimum-cardinality floor.

The cryptography guarantees anonymity unconditionally up to the moment of decoding: no one without the factorization — not the public-facing servers, not a log leak, not a subpoena of front-end logs — can link a client’s tokens or read its values at all. Everything *beyond* raw decoded logs, however, depends on the decoder operating correctly. This is the one place the protocol asks for trust, and it is worth being honest about: the design collapses the trusted surface from “the entire serving fleet and anyone who ever touches a request log” down to a single, identifiable, offline decode step. That is a large and worthwhile reduction — but it is trust-*minimized*, not trust-*free*, and the rules below are what that one trusted step must follow.

The three levels are genuinely distinct; level 2 does not collapse into level 3. The minimum-cardinality floor is a *publish-time* operation: it suppresses a slice below threshold $T$, but two below-floor slices can have a union *above* $T$, so anything that might later be recombined must be retained in a form more revealing than anything published. You cannot both keep a freely re-sliceable cube and have it already be floor-safe. What you *can* do is bound how long the most sensitive form exists.

## Aggregate eagerly

The single most effective operational rule is to **decode and aggregate in one step, and never persist raw decodable logs.** HyperLogLog aggregation is *within-class* — for each bucket, keep the maximum geometric sample — and the instant a request is folded into its class’s sketch, the fact that this bucket in class A and that bucket in class B came from the same request is gone. The cross-class join key that the [request-bundle attack](05-security-analysis.md#Request-bundles-and-cross-class-correlation) depends on is destroyed at aggregation time. Per-class, per-slice sketches are then safe to retain; raw $(\text{class}, b, k)$ records never need to be.

Concretely: the factorization lives only on an offline decoder, which reads request logs, decodes each token, folds it immediately into the appropriate sketches, and discards the decoded record. The window in which a raw, cross-class-linkable record exists shrinks to the decode-to-fold latency. This is what bounds the request-bundle residual that semisharding does not fully remove, and it is what keeps a stable bucket cohort from ever being assembled into a per-client profile: the per-client records simply do not persist to be joined.

There is a genuine limit here, worth stating rather than papering over. Retaining only floored, published aggregates gives up ad-hoc re-slicing; retaining a freely mergeable cube keeps re-slicing but holds something more sensitive than what is published, inside the decoder’s trust boundary. Where a deployment lands on that spectrum is a deliberate choice. For Julia’s purposes a fixed menu of reports — unique clients per package, splits by OS, version, and region for popular classes, at daily and coarser granularities — is retained as floored aggregates, so that what is stored is essentially what is published.

## Key regeneration and long-term unlinkability

Semisharding, eager aggregation, and the floor all bound what can be learned *within* a window of activity. None of them bounds linkability *across* long spans, because a client’s master key $x_0$ — and thus its bucket $b_0$ and its per-class values — is otherwise stable for as long as the client keeps its stored ring. The remedy is for the client to **regenerate $x_0$ periodically**, on a fixed period $T$ (Julia uses one year).

Regeneration, rather than folding a time epoch into the class hash, is the right mechanism because it *destroys* information rather than relabeling it: once last period’s $x_0$ is discarded, that period’s values cannot be reconstructed by anyone, including the client itself. It is also the only mitigation here that is client-enforced and independent of how the server behaves — every other rule above relies on the decoder acting correctly.

To avoid a synchronized rotation — every client rolling over on New Year’s Day, which would both spike the counts and reveal that everyone rotated at once — each client draws its *first* rotation time uniformly from $(0, T]$ and rotates every $T$ thereafter, advancing the expiry by whole periods so the schedule stays on the client’s own uniformly-random phase rather than resetting to whenever it happened to notice. Install times are not uniform (releases, semester starts), so phasing on them would not do; drawing the phase explicitly guarantees uniformity regardless.

Uniform-phase, fixed-period rotation turns identity into a renewal process with a *known* survival function: the probability that a client’s identity persists across a gap $g$ is exactly

```math
\begin{aligned}
S(g) = \max\!\left(0,\; 1 - g/T\right),
\end{aligned}
```

linear, and zero at $g = T$. A retention or cohort query at lag $k$ observes a returning client only when it did not rotate in the gap — with probability $S(k)$ — so the true count of returning clients is the observed overlap divided by $S(k)$, a correction with no free parameters since $T$ is a protocol constant. As $k \to T$ the divisor vanishes and the estimate’s variance grows gracefully; past $T$ there is simply no signal.

That last point is a theorem, not an accident: **no estimate can span more time than the client samples it rests on.** Because rotation makes a client’s successive identities independent, two populations that agree on all structure within $T$ but differ in how identities group into physical clients across $T$ produce *identically distributed* sketches — so no estimator can tell them apart. A long-window unique count, or a retention curve past a year, would have to distinguish exactly those, and cannot: the information is not merely hard to extract, it is absent. The maximum linkability window and the maximum measurable retention span are therefore the same quantity $T$, and pushing one out pushes the other out with it. This is a feature — the privacy horizon and the analytics horizon are a single dial, and no analysis can outrun the privacy guarantee.

What remains measurable over arbitrarily long spans is anything built from a *sequence* of within-$T$ measurements: per-period unique counts plotted over years, growth rates, the evolving mix of operating systems or versions, seasonal patterns. What is lost is *individual continuity* beyond a period — multi-year cohort retention, a deduplicated count of humans over five years — which is exactly the individual-tracking capability the horizon is meant to deny. The loss is precisely coextensive with the privacy gain.

## Ephemeral clients

One population deserves a note in any published methodology: ephemeral clients — CI runners, containers, anything that does not persist `~/.julia` — regenerate a fresh master key on essentially every run. They inflate the raw unique-client count (each run counts as a new “client”) and blur what “unique client” even means. Julia’s Pkg client can identify most of these from environment variables (`CI`, `GITHUB_ACTIONS`, and the like) and flag or separate them, so that human and automated counts do not contaminate each other. This matters a great deal for *interpreting* the counts, and — as noted in [the cohort analysis](05-security-analysis.md#The-bucket-as-a-cohort) — essentially not at all for how anonymous a persistent human client is, since ephemeral identities form no stable cohort.
