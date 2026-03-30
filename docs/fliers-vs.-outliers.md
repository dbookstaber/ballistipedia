# Fliers vs. Outliers

**Title page "Fliers" to replace existing page of that name**

# Fliers - Defects or Outliers?

To target shooters a shot is colloquially labeled as a "flier" if the shot flies wide of the target, or if has a POI too far from the COI of the other shots on the target. Every shooter experiences fliers. The shooter may, or may not know the cause of the flier. 

# Definitions

In using statistics to analyze target precision it is very necessary to formalize the definition of a "flier" and to separate fliers into two types. 

**Definitions:**

; Flier
> A shot which is either a defective shot or an outlier.

; Defect 
> A shot that is the result of an atypical factor, the defect, that affects one or more of the shooting processes. The shooter may or may not be aware of the atypical factor.
 
; Outlier
> A shot that is at an extreme value for the distribution. In other words a shot with such a large position difference from the average of the other shots that it seems improbable.

Thus there is a subtle but significant difference between formal definitions of defect and outlier. The gist is that there is no way to guarantee that a defect can be modeled since we don't know all of the atypical factors that might cause a defect, or what sort of distribution most of those aberrations would create even if we could identify them. Worse given the nature shooting, we might not be able to directly measure the factor which is atypical. In that case the only objective evidence would be from the target analysis. However outliers can be readily mathematically modeled using target analysis.   

# Defects

In essence a defective shot is a shot which has a different distribution than "regular" shots. For the purposes of this discussion let's define the distribution of shots by their Center Of Impact (COI) and their Mean Radius (MR). So for regular shots assume:  


$$
COI_{reg}, MR_{reg}
$$


but let's assume that we have 6 different possible defects each with its own different distribution:  

    $COI_{defect(1)}, MR_{defect(1)}$  

    $COI_{defect(2)}, MR_{defect(2)}$  

    ...  

    $COI_{defect(6)}, MR_{defect(6)}$  


When we were considering the overall system error we added the various error sources as:  

    $\sigma_{System}^2 = \sigma_{Weapon}^2 + \sigma_{Ammunition}^2 + \sigma_{Shooter}^2$  

because the weapon contributes some error, then the ammunition contributes some more error, and the shooter even more. However defects are handled differently. 

you and buddy...

Let's simulate what happens with defects. Consider that 1/6 of the shots are fliers. Get a blank regular d6 die from the hobby store and mark one face with a black circle, the other faces being plain. When you roll the die, you get a regular shot if you roll a plain face, and a flier if you roll the face with the black spot. Now for the nasty rub. A hole is a hole is a hole. There is no way to tell if a hole is a "true flier shot" or a "true regular shot"!

Remember that we have used infinite distributions to model the shot distribution. So if we are given the parameters for the horizontal and vertical normal distributions, then we can readily calculate a pistol shot being 5 miles wide even though that is impossible by physics. Consider being able to "fairly" flip a coin heads a thousand times in a row. It is so impossible that you'd have to suspect cheating, but it is *possible*. So to a statistician improbable certainly does not mean impossible. To prove that a pistol shot won't travel 5 miles we couldn't use the probability models, but rather external ballistics using muzzle velocity and a drag coefficient for the projectile. This of course requires a more sophisticated distribution model, and more assumptions. 

You might think that we might be able to get the fliers to go to one group and the regular shots to another. It might work for different ammunition types, but there is no guarantee. It could be that $COI_{reg} = COI_{flier}$ but that $MeanRadius_{flier} \gg MeanRadius_{reg}$ in which case the fliers would be sprayed all around the COI. 

Obviously another problem with trying to "pattern" the two types is sampling. With 1/6 fliers (an appalling bad defect probability!) to get 5 flier shots we'd shoot on average 30 shots. To most shooters that is an exorbitant amount of shooting. 

A flier might have a known cause before the target is examined, for example:

- Benchrest-level shooters traditionally discard the first round(s) after cleaning a barrel as a "fouling" shot(s).  The friction difference between a clean and a fouled bore are enough to significantly alter the point of impact.
- A shooter may "call a flier" if he knows he committed an error that is not characteristic of his shooting.
- If the shots are being chronographed, then the shooter might "call a flier" on any shot that chronographs outside of the 95% confidence interval around the mean muzzle velocity.

However a defect (or fliers) might have an unknown cause, and might not be suspected until the shot pattern on the target is observed. 

- If the shots are not being chronographed, but one shot is very low, then the shooter would suspect an excessive muzzle velocity variation.   
  

> So to obtain objective evidence that muzzle velocity variation is the problem, the shooter would need to design an experiment to test for that process difference - for instance using a chronograph. In this situation the chronograph adds objective evidence which is not available looking at the bullet holes alone.

- A projectile might be off-balance in the distribution of mass.   
  

> So in this case the shot is a true defect, but the shooter doesn't have any means to measure the mass balance of a spinning projectile. So there would be now way to know that the projectile was a defect. (Making an example, rather than challenging experimentalists!) The only objective evidence is the bullet hole such a projectile makes - and a bullet hole is a bullet hole is a bullet hole.

The salient point is that some objective evidence of process variation must exist to be able to label a shot a defect. If the only evidence is the position of the holes in the target, then a shot can't be labeled a defect. The only way to analyze a flyer on the target when just the relative holes positions of the shots are known, is through the consideration of outliers. An outlier just being an improbable shot given the target analysis model.

## Outliers

**Not every outlier is a defect.**  

Unbounded distributions have been accepted as models for the shooting process, and so outliers are part of both the model and real world, and that our model can correctly account for them if they are part of the modeled process.  Granted, if I had a rail gun on an indoor range and had triple-checked every component of every round I sent downrange I may not accept an unbounded normal distribution as a model of my shot dispersion.  But once we allow for outdoor conditions and normal ammunition, not to mention a shooter operating the gun, then in the normal course of events we will get outliers, and they *are* representative of the underlying normally-distributed process.

It is not unreasonable to accept a model that says 1 round in 100 is going to miss the target entirely.  If we are recording statistically significant samples and using robust estimators then including such outliers will not ruin our estimates.  And in a way our metrics for "statistical significance" will tell us whether an outlier is valid.  E.g., if in my first three shots after sighting in one shot nicks the edge of the target backer then I know right away I need more samples because so far my confidence interval is wider than my target!  If I take another 20 shots and they cluster into a single hole then perhaps I can decide whether to exclude the outlier as a "flier" or incorporate it as a sample from the "true" model of my precision.
 
Ideally maybe we do want to clip our unbounded distribution models, or maybe we want to overlay our shot distribution model with a Poisson dispersion model that allows us to exclude samples that may be due to a defective round, wind gust, etc.  But practically we are already pushing the bounds on the sample size needed just to determine covariance for a general bivariate normal model, so adding a fourth parameter to the [models of dispersion](precision-models.md#models-of-dispersion) may be a bridge too far.

# Counting Fliers

** A caveat - This section isn't meant to give a detailed or reasonable statistical analysis design model, but rather is a contrived situation to illustrate some points about what in process control is called "Acceptance Testing." **

Consider two different types of ammunition, A and B, which you want to use to hunt. So you decide to test 25 shots of each on the range using the mean radius measurement. Obviously one characteristic to compare would be the precision of the two types of ammunition. When measuring precision we decided to set a clip level at the 95^th^ percentile of the mean radius. So we are going to, on average, throw out the worst 5% of the shots. 

But here is perhaps a unknown twist. You can also measure a quality factor based on what percentage of the shots are **fliers due to the ammunition.** In the above sections we took great care to separate fliers into defects and outliers and here we've combined the two again, but notice the very important qualification phrase **due to the ammunition.** So if type B ammunition has a muzzle velocity problem where some shots will have a abnormally low velocity, then such a shot is a defect. The shot would be low compared to the COI for type B ammo, so that shout wouldn't be counted in measuring the precision of B. But it is a **flier due to the ammunition.**

Labeling a shot a **flier due to the ammunition** is a non-parametric measurement. You just get a count. But that count can be useful as a quality factor. There are two problems with statistical tests based on non-parametric measurements. First, in general they aren't as precise as parametric measurements which means that larger samples are needed. Second we'd expect commercial ammunition to have a relatively low percentage of fliers. So to use non-parametric methods will probably require a lot of shooting. 

Type A has 2 fliers. One was a called flier due to shooter error. The other was an outlier. So overall Type A has 24 shots in the sample with 1 shot being a flier flier due to the ammunition. 

Type B has 3 fliers, 2 outliers and 1 defect due to low muzzle velocity. The parametric measurement (radial distance for the shot to the center) of the 2 bad shots in B radius distance allowed them to be labeled defects. So overall Type B has 25 shots in the sample with 3 shots being a flier flier due to the ammunition.

wiki table...

Based on this data (admittedly small sample) you can statistically compare the flier count of the two types of ammunition. (Based on our clip level we expect a 5% defect rate.) Type A has an expected rate of 4.2%, and type B has an expected rate of 12%. For full details details on the statistical procedure see: **yada yada.**

# Examples

To perhaps belabor the difference between defects and outliers consider the following examples. 

**Example 1** - Ten shots are shot at a paper target with ten bulls-eyes.  The cartridge cases are lined after each is shot. After shooting it is discovered that nine of the shots are ammunition type-A  and one is type-B. Shot 7 is the type-B ammunition shot.  

> *Shot 7 is a defect and should be ignored in the measurement(s). Note here that it doesn't matter where the shot hit. The only reason to allow type-A ammunition and type-B ammunition to mix would be if the two types had already been comparatively tested and found to have the same performance. Here performance doesn't just mean the same precision since the two types of ammunition could have the same precision, but have different average COI positions.*

**Example 2** - Ten shots are shot at a target (single bulls-eye). After shooting it is discovered that nine of the shots are ammunition type-A  and one is type-B. It is unknown which shot used the type-B ammunition. There is one shot which is wide of the group of the other nine. 

> *In this situation the shot with the type-B ammunition is a really defect since it isn't of the same type as the other shots. Since it is unknown which shot used the type-B ammunition, it is invalid to just throw out "wide" shot and assumed it was the shot with type-B ammunition. The shot with the type-B ammunition may in fact be the closet shot to the center of the group!*

> *The wide shot can only be labeled as an "outlier" if it falls outside of some set confidence interval. Ideally the confidence interval for acceptance should be decided upon before the experiment, and then data outside of the confidence interval would be properly rejected.*

> *So here some ad hoc judgment may be required. The best option is probably to throw out the group/measurement entirely. This would be especially true if using the measurement Extreme Spread. However if we're using the mean radius measurement then the one Type-B shot probably won't perturb the mean radius measurement too much. Thus for the mean radius measurement the solution to the predicament might be to consider the confidence interval about the measurement to decide if the wide shot should be thrown out as an outlier, and use the resulting 9 or 10 shot measurement. Such a group measure gets "an asterisk", but it would be very useful to estimate the sample size needed to get a mean radius measurement of specific precision.*

**Example 3** - Ten shots are shot at a target (single bulls-eye). There is one shot which is wide of the group of the other nine. The shooter has no idea why there is one wide shot.

> *In this case the wide shot would be an outlier if it was rejected based on some confidence limit.*

> *The nasty part here is that the wide shot might, unknown to the shooter, truly be a defect. For example in the manufacture of the cartridge, this particular cartridge might have had bad primer which was outside of the normal process variations. After shooting this would of course be impossible to determine. Even if this sort of quality problem had been suspected, it would be virtually impossible to measure for a commercial cartridge. So some ammunition manufacturing problems can not be isolated by independent measurement, but rather only a nebulous judgement is possible that the "quality" of the ammunition is "poor" based on the fact that the system variance was much larger than for other ammunition types, or the particular ammunition had a large number of outliers.*

  


---

Note on spelling: *Flier* vs *flyer* has not been well established.  We use the former spelling here because [*flyer* seems to be more commonly used to refer to leaflets and architectural features, as opposed to "things that fly"](http://www.dailywritingtips.com/flier-vs-flyer/).
