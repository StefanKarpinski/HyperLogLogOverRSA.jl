# HyperLogLog

HyperLogLog (HLL) is an algorithm for answering the question "How many unique values are there in this collection of values with repeats?" The naive way to do this, of course, would be to collect all the values and when you see a value, check if it's in the collection of seen values. But that requires saving all the encountered values, which requires $O(n)$ memory. HyperLogLog allows accurately estimating the number of unique values with $O(\log(\log(n)))$ memory—whence the name. This grows so slowly that it is, for most practical purposes, constant. In most use cases a fixed upper bound on $n$ is assumed, which makes the memory usage actually constant. Being able to count unique values with effectively constant memory, even approximately, is amazing.

It's worth emphasizing that this is not one of those algorithms that has great theoretical properties but is impractical to use in reality. HyperLogLog is extremely practical and widely used. The constant factor in the asymptotic memory usage is quite small, and relative standard error is around $1.6\%$ in a common configuration and can be made smaller with very modest increases in memory footprint.

Our approach treats each client install as a "unique value" to count and uses HyperLogLog to estimate the cardinality of the set of those unique values. Recall that in the "gold standard" client counting algorithm, each client generates a random unique value. Instead of sending that unique random value, each client can derive a HyperLogLog value from it, and send that instead. Individual HLL samples are very small—only 18 bits in a common configuration—so they don’t encode nearly enough information to uniquely identify all clients. Moreover, since clients generate these values independently, the [birthday paradox](https://en.wikipedia.org/wiki/Birthday_problem) tells us that collisions occur long before there are $2^{18}$ clients, and since HLL values are highly non-uniform, we actually start getting collisions around as few as 64 clients! This is a promising start but, as we'll see, there are some issues.

Before we get to the problems with naively using HyperLogLog for counting clients, let’s develop an intuition for how and why HyperLogLog works. Suppose each client draws a random geometric sample value, *i.e.*

- ``0`` has probability $\frac{1}{2}$
- ``1`` has probability $\frac{1}{4}$
- ``2`` has probability $\frac{1}{8}$
- ``k`` has probability $\frac{1}{2^{k+1}}$

The basic intuition of HyperLogLog is that if the maximum geometric sample observed is $k$ then there are around $2^k$ unique values in the set. Of course, this estimator is terrible. It only ever estimates power-of-two cardinalities and has huge variance. HyperLogLog addresses both problems by uniformly splitting the space of clients into buckets, making a geometric estimate for each bucket, and then combining all those estimates into a single improved estimate. With enough buckets this actually produces a high quality estimate. The number of necessary buckets is not even inordinately large—in the common configuration we've mentioned, there are $2^{12} = 4096$ buckets, which gives relative standard error of $1.63\%$, corresponding to a 95% confidence interval of $±3.19\%$. If more accuracy is desired, it's relatively cheap to add more buckets and reduce the error bounds—the error is inversely proportional to the square root of the number of buckets.

Each HyperLogLog sample has two parts:

- `b`: a uniformly distributed **bucket index**
- `k`: a geometrically distributed **sample value**

The common configuration we've referred to above has $2^{12} = 4096$ buckets and caps the geometric sample value at $63$ so that the geometric value can be encoded with $6$ bits and the estimates are capped at $2^{63}$. This is plenty for many purposes, including ours. Each such HyperLogLog value can be encoded with only $18$ bits: $12$ bits for the bucket and $6$ bits for the geometric sample. If we want to splurge and allow estimating cardinalities up to $2^{127}$ we need only a single extra bit per HLL value. If we want to make our estimates twice as accurate, we have to quadruple the number of buckets, which only costs two extra bits per HLL sample. HyperLogLog's extreme memory efficiency means every additional bit has great leverage. We will use the "typical" 18-bit configuration as a reasonable reference point throughout.

A straightforward implementation of generating a random HLL pair is:
```julia
b = rand(UInt16) >> 4
k = min(63, trailing_zeros(rand(UInt64)))
```
Suppose each client generates an 18-bit HLL sample value like this once, saves it, and then includes that value in the header of each Pkg request it makes. Then we can estimate the number of unique clients appearing in any subset of request logs by applying HyperLogLog aggregation and estimation to the HLL values for that set of logs. Since each client's HLL value is only 18 bits, it's impossible to uniquely identify most users—there just aren't enough bits to distinguish more than $2^{18} = 262144$ different users—and since they generate values independently, collisions will start to happen long before that many users.

This seems pretty good. Can we just do this? I have, frankly, considered it. But there are two main issues that have given me pause:

1. Most clients get small geometric sample values that are very common and thus anonymous, but some clients will get rare geometric values and be highly identifiable. HLL samples are anonymous *on average*, but an inherent aspect of how the estimator works is that some values are rare and thus not anonymous.
2. HyperLogLog values are trivially forgeable. If a malicious client were to send requests with intentionally large geometric values, they could inflate the client count estimates arbitrarily high, rendering them useless. Worse: a malicious client only needs to send one forged request per bucket to push the estimate to the maximum.

Using HyperLogLog by itself for unique client estimation, while tantalizing, is not as anonymous as we’d like and requires a trivial amount of work to sabotage. Fortunately, both of these issues can be addressed.
