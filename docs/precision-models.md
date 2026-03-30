# Models of Dispersion

We present five models for measuring and analyzing precision:

1. [Closed Form Precision](closed-form-precision.md)
1. [Circular Error Probable](circular-error-probable.md)
1. [Elliptic Error Probable](elliptic-error-probable.md)
1. [Range Statistics](range-statistics.md)
1. [Order Statistics](order-statistics.md)

Before selecting one consider the following background:

## General Bivariate Normal
[The normal, a.k.a. Gaussian, distribution](http://en.wikipedia.org/wiki/Normal_distribution) is the broadly accepted model of a random variable like the dispersion of a physical gunshot from its center point.  The normal distribution is parameterized by its mean and standard deviation, or $(\mu, \sigma)$.  As explained in *[What is Precision?](what-is-precision.md)* we are only interested in the dispersion component, since the center point of impact is controlled by [sighting in the gun](faq.md#how-many-shots-do-i-need-to-sight-in) (i.e., adjusting its aiming device).  Therefore we will assume that a gunner can dial $\mu \approx 0$ and leave that parameter out of the question in what follows.

Since we are interested in shot dispersion on a two-dimensional target we will look at a [bivariate normal distribution](http://en.wikipedia.org/wiki/Bivariate_normal_distribution), which has separate parameters for the standard deviation in each dimension, $\sigma_x, \sigma_y$, as well as a correlation parameter *ρ*.

## Uncorrelated Bivariate Normal
We don't have any evidence that there is, or should be, correlation between the horizontal and vertical dispersion of gunshots.  Therefore, throughout our analysis we will assume *ρ* = 0.

We do know that targets can often exhibit vertical or horizontal stringing, and therefore $\sigma_x \neq \sigma_y$.  To the extent these parameters are not equal they produce [elliptic](elliptic-error-probable.md) instead of circular shot groups.

However, we know some of the significant sources of stringing and can potentially factor them out:

1. The primary source of x-specific variance is crosswind.  If we measure the wind while shooting we can bound and remove a “wind variance” term from that axis.  E.g., *Suppose the orthogonal component of wind is ranging at random from 0-10mph during the shooting.  Given lag-time *t* this will expand the no-wind horizontal dispersion at the target by $\sigma_w$.*[^1]  Since variances are additive we could adjust $\sigma_x$ via the equation ${\sigma'}_x^2 = \sigma_x^2 - \sigma_w^2$.
1. The primary source of y-specific variance is muzzle velocity, which we can actually measure with a chronograph (or assert) and then remove from that axis.  E.g., "If standard deviation of muzzle velocity is $\sigma_{mv}$ then, given the bullet's ballistic model for the given target distance, the vertical spread attributable to that is some $\sigma_v$.  Here too we can remove this known source of dispersion from our samples via the equation ${\sigma'}_y^2 = \sigma_y^2 - \sigma_v^2$.  This adjustment is shown in several of the examples:
    - [22LR CCI 40gr HV 40-shot 100-yard Example](22lr-cci-40gr-hv-40-shot-100-yard-example.md)
    - [300BLK Subsonic 20-shot 100-yard Example](300blk-subsonic-20-shot-100-yard-example.md)

# Statistical Analysis of Dispersion

In view of the preceding:

1. The [Closed Form Precision](closed-form-precision.md) model requires that we assume the shot group is, or can be normalized to be, a fairly symmetric bivariate Gaussian process.  This assumption is the most amenable to statistical analysis.
1. [Order Statistics](order-statistics.md) are slightly less efficient and amenable to abstract analysis, but are both more robust and easier to apply "in the field."
1. [Circular Error Probable](circular-error-probable.md) disregards any ellipticity in the actual shot process in order to characterize precision using a single parameter.  Since most of precision estimation is for the purposes of comparing loads, rifles, and shooters, we need a single number and we don't care if the dispersion is elliptic: tighter is always better.
1. [Elliptic Error Probable](elliptic-error-probable.md) allows for a full characterization of the General Bivariate Normal model.  For some applications – e.g., computing hit probabilities on non-circular targets – we want to preserve statistically significant ellipticity.
1. Extreme Spread and the other [Range Statistics](range-statistics.md), which increase with number of shots per group *n*, do not have any useful functional forms.  The characteristics of these measures have to be derived from Monte Carlo simulation.  They are the least efficient statistics but are also the most commonly used because they are so easy to measure in the field and so familiar to shooters.

One practical question that many shooters raise is what to do with outliers, known in the sport as "fliers."  We address [fliers here](fliers.md).

# Tools
See [Measuring Tools](measuring-tools.md) for convenient ways of measuring and analyzing precision.

# References


[^1]: Wind deflection is a function of the ballistic curve and distance, but can be expressed as a simple product of the cross-wind velocity and lag time.  For more information on the "lag rule" see Bryan Litz, ''Applied Ballistics for Long Range Shooting, 2<sup>nd</sup> Edition'' (2011) A4; or Robert McCoy, ''Modern Exterior Ballistics, 2<sup>nd</sup> Edition'' (2012) 7.27.
