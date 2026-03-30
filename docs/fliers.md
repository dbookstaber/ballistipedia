# Fliers

In target shooting a "flier" is a shot that flies wide of the target, or appears to be an outlier.  Every shooter experiences fliers.  Sometimes the cause of a flier is known or can be found.

But in statistics we should be careful to discriminate between "outliers" and "fliers."

**A priori we should reserve *flier* to describe a sample that is not characteristic of the underlying process.**  E.g.,
- For benchrest-level shooting we traditionally discard the first round after cleaning a barrel as a "fouling" shot.  The difference between a clean and a fouled bore in terms of friction and fit are enough to significantly alter the point of impact.
- We can allow a shooter to "call a flier" if he knows he committed an error that is not characteristic of his shooting.
- If we are focusing on a gun’s inherent precision we might "call a flier" on any shot that chronographs less than 90% of the mean velocity of that batch of ammo.

**But not every outlier is a flier.**  We have accepted unbounded distributions as models of the shooting process, and so we have to also accept that outliers are part of both the model and the real world, and that our model can correctly account for them if they are part of the modeled process.  Granted, if I had a rail gun on an indoor range and had triple-checked every component of every round I sent downrange I may not accept an unbounded normal distribution as a model of my shot dispersion.  But once we allow for outdoor conditions and normal ammunition, not to mention a shooter operating the gun, then in the normal course of events we will get outliers, and they *are* representative of the underlying normally-distributed process.

It is not unreasonable to accept a model that says 1 round in 100 is going to miss the target entirely.  If we are recording statistically significant samples and using robust estimators then including such outliers will not ruin our estimates.  And in a way our metrics for "statistical significance" will tell us whether an outlier is valid.  E.g., if in my first three shots after sighting in one shot nicks the edge of the target backer then I know right away I need more samples because so far my confidence interval is wider than my target!  If I take another 20 shots and they cluster into a single hole then perhaps I can decide whether to exclude the outlier as a "flier" or incorporate it as a sample from the "true" model of my precision.
 
Ideally maybe we do want to clip our unbounded distribution models, or maybe we want to overlay our shot distribution model with a Poisson dispersion model that allows us to exclude samples that may be due to a defective round, wind gust, etc.  But practically we are already pushing the bounds on the sample size needed just to determine covariance for a general bivariate normal model, so adding a fourth parameter to the [models of dispersion](precision-models.md#models-of-dispersion) may be a bridge too far.

  


---


Note on spelling: *Flier* vs *flyer* has not been well established.  We use the former spelling here because [*flyer* seems to be more commonly used to refer to leaflets and architectural features, as opposed to "things that fly"](http://www.dailywritingtips.com/flier-vs-flyer/).
