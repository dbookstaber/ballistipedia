# Home

This site explains and demonstrates statistics for analyzing the precision of projectile weapon systems.[^1]

High level topics, which are good places to start exploring the site, include:

- [What is Precision?](what-is-precision.md): An important explanation of the difference between precision and accuracy as the terms are used in statistics
- [Why BPC](why-bpc.md): How people go wrong evaluating precision
- [Describing Precision](describing-precision.md): Units, terms, and relationships
- [Precision Models](precision-models.md): Statistical approaches for efficient estimation and inference of precision
- [Prior Art](prior-art.md): Reviews of past efforts to address this question
- [FAQ](faq.md)
- [Ballistic Precision Classification](ballistic-precision-classification.md): A proposed industry standard for determining and describing precision

# Synopsis

When testing a gun, shooter, and/or ammunition the most popular measure is [Extreme Spread](describing-precision.md#extreme-spread) or "group size" of a sample of target shots.  However Extreme Spread must be used with care since it is frequently and easily abused.[^2] As with all measures, the single best measurement is meaningless in isolation. The proper statistical estimator is an "average" (a.k.a., *expected value*) of the measurement.

Another consideration is that some measures, such as Extreme Spread, change value when there are more shots in a group.  Measures that have such a dependency will be referred to here as *variant measures*.  There are *[invariant measures](describing-precision.md#invariant-measures)*, like [Circular Error Probable](circular-error-probable.md) or Mean Radius, for which the expected values do not change with the number of shots in a group. Instead [having more shots only increases the confidence in the measure's value](closed-form-precision.md#how-large-a-sample-do-we-need). Of course the experimental error of either type of measurement can also be decreased by increasing the sample size (i.e., shooting more groups). 

Furthermore, by first making assumptions about the inherent shot dispersion, then it is possible to use theoretical models to estimate measurements and their precision.  The distributions are of two basic types:  If the expected values and the expected precision factor for the measurements depend on distributions which have an [explicit solution](closed-form-precision.md) then the values can be calculated formulaically.  If the values don't have a distribution with a closed form expression then they can be estimated via Monte Carlo approaches.
 
Examples of the application of these [methods](precision-models.md) and [tools](measuring-tools.md) include:

- Determining how many sighter shots you should take.
- Determining the likelihood of a hit on a particular target by a zeroed shooting system.
- Comparing the inherent precision of different shooting systems.
- Determining which ammunition shoots better in a particular gun.

----


[^1]: <small>Typical examples would be target shooting with a rifle or pistol. Such weapons as shotguns, mortars, and ballistic missiles would have some similar characteristics, but also have factors that are neglected in the discussions and measurements.  The wiki will discuss some factors of ballistics, but it is not intended to address all the nuances of internal, external, or terminal ballistics. Rather, the focus is primarily on the analysis of the precision of the whole weapon system which can be observed directly by the relative impact points on a target. Some effort will be made to explore the precision of weapons subsystems.</small>

[^2]: <small>Many references are worth reading for further background on how Extreme Spread is broadly misunderstood. Recommended include [Understanding Rifle Precision](https://www.thetruthaboutguns.com/understanding-rifle-precision/), [Thinking Statistically](https://www.autotrickler.com/blog/thinking-statistically), and [Rifle Ammunition Load Workup](http://www.bisonops.com/2019/08/17/rifle-ammunition-load-workup/).</small>
