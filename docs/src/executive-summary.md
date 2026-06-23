# Executive Summary

I'm proposing a new (as in novel) way to estimate how many unique clients use Julia—and each package, sliced by OS, version, and region—*without* being able to track anyone. It stacks three layers:

1. Resource-class sharding
2. HyperLogLog sampling
3. RSA encryption

The first two layers provide the privacy: you can't follow a client across packages, and every client hides in a crowd of lookalikes in most classes. The RSA layer doesn't change the privacy model, but stops clients from faking their samples to skew the counts. So if you accept the privacy of sharding + HLL and believe two short proofs about the RSA layer, then you should be ok with the whole protocol.

I'd like to implement the client side of this protocol in Julia 1.14, and then use the data the server collects to publish client-count estimates alongside other Pkg server stats.

## Overview

We'd like to estimate how many people use Julia—and how many use each individual package—without tracking anyone. This isn't idle curiosity: open source funding and project survival often hinge on demonstrating impact ("installed on 50,000 systems last year"), and package authors deserve to know whether anyone is actually using their work. We'd also like aggregate breakdowns by operating system, Julia version, and region. The obstacle is that the obvious way to count unique installs—have each client send a persistent random ID—was rightly rejected as a privacy hazard: it would let anyone with access to Pkg server logs follow every request a client ever makes.

Somewhat surprisingly, there isn’t any really good technique for this in the literature—at least not that I or my research assistants, Claude and GPT, could find. To address this gap (and frankly to scratch an annoying itch), I’ve developed a new protocol for counting clients while keeping individuals anonymous. I’m tentatively calling it HyperLogLog Over RSA, and I’m posting here to give people a chance to look it over, comment, and raise any objections they might have to using it in Pkg.

The protocol’s privacy model is built in three layers, and the useful thing is that you can evaluate them independently. The headline result is this: **the privacy of the full system is exactly the privacy of the first two layers.** The third layer—the RSA cryptography—changes nothing about what the server learns about a client; it only keeps the clients from messing with the counts. So if you're comfortable with layers 1 and 2, and you trust the two proofs about layer 3, then you're comfortable with the whole protocol.

## Layer 1 — Resource class sharding

Instead of one identity per client, each client uses a *separate, unlinkable* identity for each "resource class": the top-level registry, each registry, each package, each artifact. Concretely, a client keeps a single master key and derives its per-class value by cryptographically hashing that key with the class name—so it never has to store anything per class. Because the classes are cryptographically independent, the server can count unique clients *within* a class but can never link a client's activity *across* classes. You cannot follow a user from one package to another.

This already removes the scariest capability of the rejected unique-ID scheme—cross-package tracking—and it would help even if we kept full, unique per-class IDs. Sharding on its own buys real privacy. What it doesn't fix is that, within a single class, each client would still be uniquely identifiable if we used unique IDs. But we don’t have to do that to get decent estimates.

## Layer 2 — HyperLogLog

HyperLogLog is a pretty famous cardinality estimation algorithm. In our context we can use it to provide on-average anonymity within resource classes. We replace a large, unique ID with a tiny, deliberately lossy sample: a uniformly random 12-bit bucket index plus a small (6-bit), geometrically distributed value. That's about 18 bits of information about each client. In a class with many more clients than there are buckets ($B \approx 4096$), that is far too little to pick anyone out: the average client shares its sample with a healthy crowd of lookalikes and is genuinely anonymous. HyperLogLog's near-magical property is that one can still estimate the number of *unique* clients from these samples, to a tunable accuracy, even though it does not seem like there are nearly enough bits for that.

The honest catch is in the tail. The geometric value is heavily skewed: most clients draw one of the small, common values and vanish into the crowd, but a few unlucky clients draw a large, rare value that almost nobody else shares—and that makes them essentially unique within the class. So within a single class, HyperLogLog is anonymous *on average*: the typical client is well hidden, but the most identifiable clients are not.

This is exactly the unfairness that sharding (Layer 1) cures. Because a client's sample is independent in every class, the clients who stand out *here* are not the ones who stand out *there*: every client draws a rare, identifiable value in a few of its many classes and a common, anonymous one in the overwhelming majority. So instead of a fixed unlucky minority being trackable everywhere, each client is conspicuous in a handful of its classes and invisible in the rest—and even where it is conspicuous, that only links its own requests *within that one class* (*i.e.* one package's update cadence), never across classes. Sharding turns "anonymous for the average client" into "anonymous on average (across classes) for everyone."

Two honest limits remain. A class—or a fine-grained slice of one—with only a handful of clients in it can't hide anyone, no matter the scheme: there's simply almost no one to blend in with. You can probably enumerate every FreeBSD user from the OS field alone (sorry, Alex). That's inherent to counting a small population rather than anything specific to HyperLogLog, and we handle it by declining to report counts for slices below a minimum size. And all of the above concerns the value a client sends, not network-level identifiers like its IP address—those carry their own tracking risk and must be handled separately, as they are today.

If you're still reading and you can accept the privacy model of HyperLogLog counting with resource class sharding, you've accepted the entire privacy model. Everything below is about integrity—ensuring that badly behaved clients cannot mess with the estimates too much—not privacy.

## Layer 3 — RSA, for integrity only

The two layers above completely define what the server can learn. The RSA layer's sole job is to stop *misbehaving clients* from skewing the counts. With plaintext HyperLogLog, one client could spoof a handful of requests carrying large geometric values and inflate an estimate arbitrarily. RSA fixes this by letting clients sample a HyperLogLog value *without being able to see or choose it*—so the only way to push a count up is to send more requests, at a near-optimal cost of about 1.44 requests per unit of inflation (the naive unique-ID scheme's cost is exactly 1).

Two proofs make the privacy promise precise:

1. **Anonymity is preserved exactly.** The encryption and the per-request re-randomization reveal *no more* about a client than a plaintext HyperLogLog value would. Layer 3 cannot make privacy any worse than layers 1 and 2 already are. (Proven in [Section 9](09-proof-of-anonymity.md).)
2. **The server can't smuggle in a fingerprint.** A maliciously constructed RSA modulus *could* secretly tag clients, so the server must publish a certificate proving its modulus is "fingerprint-free." Clients verify this certificate before using a ring. (Proven in [Section 12](12-malicious-servers.md).)

If those proofs hold—and reviewing them is exactly the scrutiny this writeup invites—then RSA is purely a guardrail on honest counting, invisible to the privacy analysis.

## Two modest costs

Compared to the rightly rejected unique-ID scheme, this approach gives up two things, neither of which matters much in practice:

1. Counts can't be aggregated *across* resource classes. You can estimate the users of package A and of package B separately, but not of "A or B". There's little demand for cross-package correlations anyway, so this is no great loss to us.
2. Estimates are approximate rather than exact, with small and tunable error (more buckets ⇒ tighter error, at a few extra bits per sample). This is also ok since we don't need to know the exact number of clients—being within a few percentage points is fine.

## Why this is worth doing

For Julia specifically, this would let us estimate total active installs, per-package user counts (real impact numbers for package authors and the funders who support them), and breakdowns by OS, version, and region. Because the design is built to be privacy-preserving, we can make it opt-out rather than opt-in, while still providing a way to turn it off for those who don't want it. Opt-out is fairly necessary for this kind of thing, since most people leave the default alone (though we don't know exactly how many). A persistent identifier like a client UUID counts as personal data under regulations like the GDPR and CCPA, so a counting scheme we can turn on by default has to be one that genuinely protects users; that is precisely what this design aims to be.

If you want to dig into the details: HyperLogLog is covered in [Section 2](02-hyperloglog.md), resource class sharding in [Section 3](03-resource-class-sharding.md), and the two proofs in [Section 9](09-proof-of-anonymity.md) (anonymity) and [Section 12](12-malicious-servers.md) (fingerprint-freedom).
