# Ragged holes

It is possible — and with larger shot groups and/or more precise rifles increasingly likely — that the target will exhibit a "ragged hole" into which one or more shots has disappeared.  In statistics this is known as a "center-censored" sample.  [A detailed explanation and discussion of the problem is here](http://stats.stackexchange.com/questions/76533/estimating-variance-of-center-censored-normal-samples).  This frustrates the direct calculation of [Invariant Measures](describing-precision.md#invariant-measures) like Mean Radius and Variance, which require the location of each and every shot on a target.  Alternatives:

- [Range Statistics](describing-precision.md#range-statistics) can still be measured directly from a center-censored target.
- [Order Statistics](order-statistics.md) can also provide good estimates of precision in most cases.
- Another alternative is to [use CEP](closed-form-precision.md#summary-probabilities) to back into an estimate of *σ*:
- # Find the radius *r* of the smallest circle that covers the ragged hole.
- # From *n*, the number of shots fired, and *c* the number outside the circle, calculate the proportion of shots inside that covering circle *p=(n-c)/n*.  This turns the target into a sample with CEP(*r*).
- Avoid ragged holes by shooting fewer shots per target (meaning per "point of aim" since some paper targets have more than one POA), and then aggregate them into a data set.  Software to automate this process exists, e.g. [OnTarget Target Data System](http://ontargetshooting.com/tds/) offers coded target sheets that can be scanned to automatically aggregate shots across multiple aiming points.
