# Ballistic Precision Classification

## Introduction

Ballistic Precision Classification™ is a mathematically rigorous system for describing and understanding the precision of ballistic tools like rifles.  It provides information that is:

- Useful and easy for consumers to understand
- Straightforward for enthusiasts and builders to calculate
- Statistically sound enough for experts to validate
*For more background see [Why BPC](why-bpc.md).*

All firearms and components can be assigned a Ballistic Precision Class™ (BPC™), which fully characterizes their accuracy potential.  Lower numbers are better.  In practice, the lowest possible BPC is 1.  (A theoretically perfect rifle system that always puts every shot through the same hole would be BPC 0.)

Behind the scenes, BPC is defined by a single statistical parameter known as sigma (σ).  Using [the associated statistical model](closed-form-precision.md#symmetric-bivariate-normal-shots-imply-rayleigh-distributed-distances) (known as the Rayleigh distribution), statisticians can not only determine the BPC for a particular gun or component, but also compute its expected shooting precision.[^1]

| Shots within |  |  |  |  |
| --- | --- | --- | --- | --- |
| Class 1 | < 0.1MOA | Rail guns | 96% | 100% |
| Class 2 | < 0.2MOA | Benchrest guns | 54% | 96% |
| Class 3 | < 0.3MOA | Mil-spec for PSR | 29% | 75% |
| Class 4 | < 0.4MOA | Competitive auto-loaders | 18% | 54% |
| Class 5 | < 0.5MOA | Mil-spec for M110 and M24 | 12% | 39% |
| Class 6 | < 0.6MOA | Mil-spec for infantry rifles and ammo | 8% | 29% |


We can also generate the expected values of more familiar measures, like the extreme spread of a 3- or 5-shot group:

| Median 3-shot |  |  |  |  |
| --- | --- | --- | --- | --- |
| Class 1 | 100% | 0.3MOA | 100% | 0.2MOA |
| Class 2 | 98% | 0.6MOA | 99% | 0.5MOA |
| Class 3 | 65% | 0.9MOA | 85% | 0.7MOA |
| Class 4 | 26% | 1.2MOA | 57% | 0.9MOA |
| Class 5 | 9% | 1.5MOA | 35% | 1.2MOA |
| Class 6 | 3% | 1.8MOA | 21% | 1.4MOA |


### Understanding MOA
Accuracy is described in [angular terms](angular-size.md), most commonly using the unit "Minute of Arc" (MOA).  One arc minute spans 1.047" at 100 yards.  Rifle shooters often practice on 100 yard targets, and so they often think in terms of how wide their groups are at 100 yards.  People often just round it off and think of 1 MOA as "one inch at 100 yards."

In the absence of an atmosphere, the angular precision measured at one distance would be valid at all other distances.  I.e., a 1" group at 100 yards would measure 5" at 500 yards.  However, in reality the effects of wind and drag (which, together with gravity, accentuates variations in muzzle velocity) will only increase the angular spread of ballistic groups as distance increases.  Therefore, **in practice one should expect worse-than-advertised accuracy when shooting at longer distances**.  The distance at which atmosphere begins to significantly affect precision depends on a bullet’s muzzle velocity and ballistic coefficient.  For high-power rifles this is usually beyond 100 yards.  Guns shooting subsonic projectiles can begin to suffer after just 25 yards. 

### Nomenclature

BPC™ is only meaningful when it conforms to established terms and conventions.  BPC must be determined by testing in accordance with [the BPC Protocol](#protocol) described in this document.

Ballistic Precision Classification™ must be supported by the following descriptive parameters:

1. Product tested.  (Any component, or group of components, associated with accuracy can be tested.)
1. Configuration tested.  This must include the following details:
    1. Barrel length, material, profile, rifling.  
(E.g., *20" stainless 1" bull contour with 6-land 1:10"-twist cut rifling*.)
    1. Receiver, action, and feed mechanism.  
(E.g., *AR-10 magazine-fed semi-automatic*.)
    1. Ammunition.
    1. # If commercial this must include brand, model, and lot.  
(E.g., *Federal GM308M Lot#214374H077*.)
    1. # If custom, this must list component and load formula.  
(E.g., *Lapua .308 full-sized brass, WLR primer, 168gr SMK, 44gr Varget, 2.800" COAL*.)

1. Confidence in BPC.  When not conspicuously mentioned, it is assumed that the upper **90%** confidence value of sigma is referenced for BPC.

BPC measures are intentionally kept somewhat coarse.  Care should be taken to avoid suggesting more than two significant digits of precision.  For example, even after shooting 20 rounds through a gun, [the 80% confidence interval on the precision estimate typically spans 0.9-1.2 times the estimated value](closed-form-precision.md#how-large-a-sample-do-we-need).

### Trademarks
The following terms are trademarks of Boniface LLC.  They are free to use so long as their use complies with the Nomenclature and Protocol outlined here.

- Ballistic Precision Classification™
- Ballistic Precision Class™
- BPC™
The trademarks are claimed solely for the purpose of maintaining the integrity of the system and avoiding market confusion.

## Theory

### Statistical Model

BPC™ assumes that the impact of ballistic shots on a target are normally distributed with the same variance along any axis.  (Empirical data validate this assumption, and it should be true as long as atmospheric effects are negligible.[^2])  Therefore, [we use the Rayleigh distribution to model the radius r](closed-form-precision.md#symmetric-bivariate-normal-shots-imply-rayleigh-distributed-distances), or dispersion of each shot, from the center of impact.  When the coordinates of the shots have independent $N(0,\sigma)$ distributions along orthogonal axes, the radius of each shot is described by the Rayleigh probability density function:

$$
f(r,\sigma)=\frac{r}{\sigma^2}e^{−r^2/2\sigma^2}
$$

The [unbiased estimator for the parameter σ](closed-form-precision.md#rayleigh-estimates) comes from $\widehat{\sigma^2} = \frac{\sum r_i^2}{2(n−1)}$, with confidence $\widehat{\sigma^2} \sim \frac{\sum r^2}{\chi_{2n-2}^2}$.

### Simulation
Monte Carlo simulation is adequate for studying and characterizing precision.  In fact, many of the results associated with BPC, like the distribution of the extreme spread of a particular number of shots, can only be produced through simulation.

For simulation purposes random shots should be generated as (x, y) coordinates, where $X,Y \sim N(0,\sigma)$.  It is critical to "forget" the known center when using simulated data.  When shooting a real gun we never get to know the true center of impact, and instead have to use the sample center.  Likewise, Monte Carlo simulations must not reference the known 0 center, and should instead only reference the sample center of whatever group size is being studied. 

## Protocol
The Ballistic Precision Classification™ shall be the upper bound of the 90% confidence range on estimated sigma, in units of MOA, multiplied by 10 and rounded to the nearest integer.  For example, if the 90% confidence value for sigma on a tested gun is 0.47MOA, then the BPC value is 10 ✕ 0.47 = 4.7, rounded = 5.  I.e., in this example we are saying with 90% confidence that the tested gun’s accuracy is no worse than Class 5.

### Classifying a Specimen
It is not realistic to assign a BPC with fewer than 10 shots.  The 90% confidence range with just 10 shots will typically extend to 1.4 times the estimated accuracy value.  It will typically take 20 shots to get the outside bound of the 90% confidence interval to within 20% of the estimated value.[^3]

You must not discard data points during testing except for an unrelated failure.  (E.g., if you are testing a barrel and you encounter a squib load, that shot may be excluded.  But "fliers" should not generally be excluded.)

**Data:** Target distance, and (*x,y*) coordinates of the center of each shot impact.  All shots with the same point of aim must be grouped, but multiple groups can be used.

**Calculations:** The formulas to transform these data into a confidence interval for sigma, and the corresponding BPC, are shown in [BallisticAccuracyClassification.xlsx](media/BallisticAccuracyClassification.xlsx).  Given a sample of *n* shots, over *g* groups, at a distance *d*:
1. For each measurement, convert to units of MOA.  For example, if measurements are taken in inches, and the target was shot at a distance of *d* yards, then divide each measurement by 0.01047*d*
1. For each group *g*, calculate the center of the group as $(\bar{x}_{i \in g},\bar{y}_{i \in g})$
1. For each shot *i*, find its radius squared relative to the center of the group as $r_i^2=(x_i−\bar{x}_g)^2+(y_i−\bar{y}_g)^2$
1. Calculate the upper 90% confidence value for sigma as `σ~U~=SQRT[SUM(r2)/CHIINV(0.95,2n-2g)]`
1. Calculate the Ballistic Precision Class as `=ROUND(10 × σ~U~,0)`

## Notes


[^1]: The proportion of shots expected to fall within radius ''rσ'' of the center of impact is <math>1-e^{-r^2/2}</math>.

[^2]: There are two atmospheric effects that can finally create excess variance in one axis:  Variable wind will increase horizontal variance.  This is a function of the wind speed and the projectile's time-of-flight to the target.  Also, as time-of-flight increases, variations in muzzle velocity will appear as excess vertical variance.  Since both of these effects depend on time-of-flight, they are typically negligible for high-velocity rifles before 100 yards, or for subsonic projectiles before 25 yards.

[^3]: The U.S. Army Marksmanship Unit (AMU) has long used a minimum of 3 consecutive 10-shot groups fired from a machine rest to test the accuracy of service rifles.
