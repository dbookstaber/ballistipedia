# Projectile Dispersion Classifications

Before considering the measurements that will be used for the actual statistical analysis, let's consider the assumptions about projectile dispersion about the Center of Impact (COI) and how sets of those assumptions might be grouped into different classifications. The various classifications will offer insight as to the fundamental patterns expected for shots and insights to the interactions of various measures. Thus an understanding of the basic assumptions about projectile dispersion is key in being able to effectively use the measures. 

The COI is the only true point of reference which can be calculated from the pattern of shots on a target. Thus the COI is the reference point for precision measurements. The overall error that we are interested in measuring is the sum of all the various interactions that make multiple projectiles shot to the same point of aim (POA) disperse about the COI.  

Since we are primarily interested in the dispersion relative to the COI, the overall assumption is that the weapon could be properly sighted so that the COI would be the same as the POA. In practice this is achieved by [ adjusting the weapon's sights](faq.md#how-many-shots-do-i-need-to-sight-in). Thus in order to isolate projectile dispersion, all of the factors of internal and external ballistics that cause a bias to the COI on a target will be ignored. For example, for the purposes of classifying projectile dispersion, accuracy errors due to POA errors will be ignored.  

The [Normal distribution](http://en.wikipedia.org/wiki/Normal_distribution) is the broadly assumed probability model used for a single random variable and it is characterized by its mean $(\bar{x})$ and standard deviation $(\sigma)$. The [central limit theorem](http://en.wikipedia.org/wiki/Central_limit_theorem) shows that for measures for the "average" shot, or averages of multiple targets are used, then for "large" samples the averages will conform to Normal distribution even if the fundamental distribution is not a normal distribution. 


<figure class="figure-right" style="max-width: 400px" markdown="span">
  ![Distribution of samples from a symmetric bivariate normal distribution.  Axis units are multiples of σ.](images/Bivariate.png)
  <figcaption markdown="span">Distribution of samples from a symmetric bivariate normal distribution.  Axis units are multiples of σ.</figcaption>
</figure>


Since we are interested in shot dispersion on a two-dimensional target we will assume that the horizontal and vertical dispersions of the population of shots are each Normal distributions. Thus the horizontal dispersion will have mean $\mu_H$ and standard deviation $\sigma_H$. The vertical dispersion will have mean $\mu_V$ and standard deviation $\sigma_V$. Then a further assumption is made by assuming that the two dimensional expansion of the Normal distribution the [Bivariate Normal distribution](http://en.wikipedia.org/wiki/Multivariate_normal_distribution#Non-degenerate_case), applies. This adds an additional term the  [correlation parameter *ρ*](http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient). (See also: [What is ρ in the Bivariate Normal distribution?](what-is-rho-in-the-bivariate-normal-distribution.md)) Thus the expectation is that distribution should then describe, the dispersion of a gunshots about the COI, ($\mu_H$ and $\mu_V$). The full bivariate normal distribution is thus:  

    $f(H,V; \mu_H, \mu_V, \sigma_H, \sigma_V, \rho) =       \frac{1}{2 \pi  \sigma_H \sigma_V \sqrt{1-\rho^2}}       \exp\left(         -\frac{1}{2(1-\rho^2)}\left[           \frac{(H-\mu_H)^2}{\sigma_H^2} +           \frac{(V-\mu_V)^2}{\sigma_V^2} -           \frac{2\rho(H-\mu_H)(V-\mu_V)}{\sigma_H \sigma_V}         \right]       \right)$

where:  

    $-1 &le; \rho &le; 1$  

    $\sigma_H>0$ and $\sigma_V>0$

Note that the above restrictions are not additional restrictions on the model, but rather simply pointing out how the mathematics works. Thus they are more analogous to the mathematical notion that a person can't have a negative age. 

> !!! warning
    An ancillary point worth mentioning is that the assuming the Normal distribution in three dimensions leads to the [Maxwell–Boltzmann distribution](http://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution) which is the foundation of the ideal gas laws.


# Simplification of the Hoyt distribution into Special Cases

To eliminate the COI ($\mu_H$, $\mu_V$) which makes the equations "messier", a translation of the coordinate system to the COI is desired. Thus:  


    $\mu_h = 0$   and    $\mu_v = 0$  


    $h = H - \mu_H$   and    $v = V - \mu_V$  


    and  $\sigma_h = \sigma_H$   and   $\sigma_v = \sigma_V$  


This is a very pragmatic and justifiable consideration since the COI can be measured on the target, and the dispersion about the COI is the aspect of interest when measuring precision. As noted before, by adjusting the weapon's sights the POA can be made to coincide with the COI. Thus this simplification of the dispersion equations is strictly for ease of understanding as is not a limitation on the nature of the dispersion classifications. With the translation of the coordinate system to the COI, then the general bivariate normal equation becomes the Hoyt distribution:  

    $f(h,v; \sigma_h, \sigma_v, \rho) =       \frac{1}{2 \pi  \sigma_h \sigma_v \sqrt{1-\rho^2}}       \exp\left(         -\frac{1}{2(1-\rho^2)}\left[           \frac{h^2}{\sigma_h^2} +           \frac{v^2}{\sigma_v^2} -           \frac{2\rho hv}{\sigma_h \sigma_v}         \right]       \right)$

Looking at this equation two other different mutually exclusive simplifications can be readily seen:

- **Either** $\sigma_h = \sigma_v$ (equal standard deviations) **or** $\sigma_h \neq \sigma_v$ (unequal standard deviations).
> Obviously if we could measure both $\sigma_h$ and $\sigma_v$ with a very high precision (e.g 6 significant figures), then the two quantities would never really be equal. But in many cases the assumption is reasonable. In reality since shooters typically collect only a small amount of data, statistical tests will fail to detect a difference unless the difference is great. In such cases the shot pattern would be noticeably elliptical, not round.

- **Either** $\rho = 0$ (uncorrelated) **or** $\rho \neq 0$ (correlated).

> !!! warning
> !! CAREFUL !! **[Correlation does not imply causation](http://en.wikipedia.org/wiki/Correlation_does_not_imply_causation)**

> There is somewhat famous example. A researcher gathered statistics for stork sightings and births in a particular county over a twenty year period. Analysis of the data showed that over the twenty year period both stork sightings and births had increased with a very significant linear correlation. From the data you might erroneously infer that storks do bring babies!


The pair of mutually exclusive assumptions thus results in four cases for analytical evaluation as shown in the Table below. There is one case that results in circular groups, and three that result in elliptical groups. As the different in variances gets greater, or the further $\rho$ is from 0, then the ellipse will be more pronounced. 

|  |  |  |
| --- | --- | --- |
| + Group Shape vs. Assumptions (COI at Origin) |  |  |
|  | $\sigma_h \approx \sigma_v$ | $\sigma_h \neq \sigma_v$ |
| $\rho \approx 0$ | Case 1 - Circular Groups | Case 3 - Elliptical Groups |
| $\rho \neq 0$ | Case 2- Elliptical Groups | Case 4 - Elliptical Groups |


## Experimental reality of Comparing $s_h$ and $s_v$

The table above uses *approximately equal to* $(\approx)$ rather than *strictly  equal to* $( = )$. This is an acknowledgement that we are dividing the cases into ones that are close enough to be useful, even though they most certainly are not exact. To be overly persnickety there are two considerations. 

First we can only get experimental estimates from calculations based on sample data for the factors $\sigma_h$, $\sigma_v$, $\rho$ and these estimates are at best only good to a scant few significant figures. Thus even though the difference between *approximately equal to* and *strictly equal to* is under some experimental control there are practical limits. In other words, we can theoretically make the measurements as precise as we want by collecting more data, but it is quickly impractical to do so. (Assume that to double the precision that we have to quadruple the sample size. This exponential increase quickly becomes unmanageable. It is easy to pontificate about averaging over a million targets, but no one is going to shoot that many.) Thus even if $\sigma_h \equiv \sigma_v$ we'd never expect that we'd experimentally get $s_h = s_v$ due to experimental error. 

Second there is the good enough. Shooting by definition is going to have fairly small sample sizes. So if $0.66s_h < s_v < 1.5 s_h$ then, as a rule of thumb, that is probably good enough. Of course for large sample we would want to tighten the window. The harsh reality is that if $s_h$ and $s_v$ could be measured with great precision (e.g. to ten significant figures), then two values would always be statistically significantly different. 

Thus the approximation that $\sigma_h \approx \sigma_v$ will be used unless the variances are known to be statistically significantly different. On the experimental data it is possible to test for a statistically significant difference by using a ratio of $s_h^2$ and $s_v^2$ via the F-Test. The "catch" in using the F test is that the variance has poor precision for small samples. Thus the difference must be great for the F-Test to detect that the two variances are not equal.

## Simplifications Reduce Number of Coefficients to Fit

The Hoyt distribution is general enough to be able to fit all four of the special cases in the table above. The point in making special cases of the Hoyt distribution is to reduce the number of coefficients to fit to the data. In general the more coefficients to be fit, the more data is required. Also when fitting multiple coefficients some of the coefficients are determined with greater precision than others. Thus to get a "good" fit for  multiple coefficients a lot more data is required not just the minimum. 

Thus to fit the COI at least two shots are required. To fit the constant for the Rayleigh equation another shot would be required for a total of three shots. To fit the Hoyt distribution an additional five shots would be required for a total of seven shots. Probably 10 shots would be required to get a "decent" fit for the Rayleigh distribution, and at least 25 for the Hoyt distribution.

## Notation in Simplified Cases

The formulas for the distributions in the cases detailed in subsequent parts of this page are given in terms of the population parameters (i.e. $\mu_h, \mu_v, \sigma_h, \mbox{and } \sigma_v$) rather than the experimentally determined factors (i.e. $\bar{h}, \bar{v}, s_h, \mbox{and }  s_v$) on purpose to emphasize the theoretical nature of the assumptions. Of course the "true" population parameters are unknown, and they could only be estimated with the corresponding experimentally fitted values about which there is some error.

# Conformance Testing

## $\rho \approx 0$

only way linear least squares

## $\sigma_h \approx \sigma_v$

1. F-Test $\frac{s_h^2}{s_v^2}$     if    $s_h < s_v$    else    $\frac{s_v^2}{s_h^2}$
1. Studentized Ranges
1. Chi-Squared $(n-1) \frac{s^2}{\hat{\sigma}^2}$

# Circular Shot Distribution about COI

## Case 1, Rayleigh Distribution 

<figure class="figure-right" style="max-width: 250px" markdown="span">
  ![Shots dispersed about the COI. A circular dispersion is the Rayleigh distribution.](images/raleigh.jpg)
  <figcaption markdown="span">Shots dispersed about the COI. A circular dispersion is the Rayleigh distribution.</figcaption>
</figure>

Given:   

1. $\sigma_h \approx \sigma_v$  

1. $\rho \approx 0$  

then the  mathematical formula for the dispersion distribution would be the Rayleigh distribution:  

    $f(r) = \frac{r}{\Re^2} e^{-r^2/(2\Re^2)}, \quad r \geq 0,$ and $\Re$ is the shape factor of the Rayleigh distribution.  


This is really the best case for shot dispersion. Shot groups would be round. 

Strictly, for the Rayleigh distribution to apply, then  $\sigma_h = \sigma_v$, in which case $\Re = \sigma_h = \sigma_v$. For the "loose" application of the Rayleigh distribution to apply, then $\Re \approx (\sigma_h + \sigma_v)/2 \approx \sqrt{\frac{\sigma_h^2 + \sigma_v^2}{2}}$.

The following statistical measurements are appropriate:

- Circular Error Probable (CEP)
- Covering Circle Radius (CCR)
- Diagonal
- Extreme Spread (ES)
- Figure of Merit (FOM)
- Mean Radius (MR)
- Radial Standard deviation

**Notes:**
1. The Diagonal, the Extreme Spread and the FOM are different measurements, even though they conceivable could be based on the same two shots! The Extreme Spread would only depend on the two shots most distant in separation. The the Diagonal and the FOM would depend on two to four shots. For a large number of shots we'd typically expect four different shots to define the extremes for horizontal and vertical deflection.
1. For the measures for the CCR, the Diagonal, the GS and the FOM measurements a target would a ragged hole would be acceptable, but for the rest of the measures the {*h,v*} positions of each shot must be known.
1. Experimentally the radial distance for each shot, *i*, is $r_i = \sqrt{(h_i - \bar{h})^2 + (v_i - \bar{v})^2}$
1. The conversion to polar coordinates results in each shot having coordinates $(r, \theta)$. (a) The conversion implicitly assumes that the polar coordinates have been translated so that the center is at the Cartesian Coordinate of the true center of the population $(\bar{h}, \bar{v})$. (b) The distribution of $\theta$ is assumed to be entirely random and hence irrelevant. This assumption is testable. (c) The distribution is thus converted from a two-variable distribution to a one-variable distribution.

!!! warning
    Note that there is a conundrum in how we are "averaging" the horizontal and vertical standard deviations to get $\sigma_{\Re}$. Look at the two expressions. They lead to two choices, either of which may casually seem valid.
- $\sigma_h = \sigma_v$
- $\sigma_h^2 = \sigma_v^2$

In general if we look at the first formula "averaging" it leads to using:  


   $\Re = \frac{\sigma_h + \sigma_v}{2} = \sigma_h$   (with substituting $\sigma_h$ for $\sigma_v$)

However in statistics standard deviations are "averaged" (pooled) by taking the square root of the average of their variances:  

   $\Re^2 = {\frac{\sigma_h^2 + \sigma_v^2}{2}}$  

   $\Re = \sqrt{\frac{\sigma_h^2 + \sigma_v^2}{2}} = \sigma_h$  (with substituting $\sigma_h$ for $\sigma_v$)

but:


$$
\frac{\sigma_h + \sigma_v}{2} = \sqrt{\frac{\sigma_h^2 + \sigma_v^2}{2}}
$$

  if and only if $\sigma_h \equiv \sigma_v$

Thus we should take the extent that:


$$
\frac{\sigma_h + \sigma_v}{2} \neq \sqrt{\frac{\sigma_h^2 + \sigma_v^2}{2}}
$$


as a severe warning that we can not push the assumption that $\sigma_h \approx \sigma_v$ too far if we expect the simplification of the general Hoyt distribution to the Rayleigh distribution to give meaningful results. 

The situation is even more tenuous given the small samples that shooters typically use. In general the relative precision of the variance value about a mean is much less precise than the relative precision of the mean value. The statistical test to compare two experimental variance values (i.e. $\sigma_h^2, \text{and} \sigma_v^2$ in our case) is the F-Test which uses the ratios of the variances. For small samples a large difference would need to be observed before the ratio would be statistically significantly because of the imprecision in the individual experimental variance values.


# Elliptical Shot Distribution about COI

## Case 2, Equal variances and correlated 
Given:  

1. $\sigma_h \approx \sigma_v$  

1. $\rho \neq 0$  

1. The {*h,v*} position of each shot must be known.

The following statistical measurement is appropriate:

- Elliptic Error Probable

## Case 3, Unequal variances and uncorrelated (Orthogonal Elliptical Distribution) 
Given:  
 
1. $\sigma_h \neq \sigma_v$  

1. $\rho \approx 0$  

1. The {*h,v*} position of each shot must be known.
then the mathematical formula for the dispersion distribution would be:  

    $f(h,v) =       \frac{1}{2 \pi \sigma_h \sigma_v}       \exp\left(         -\frac{1}{2}\left[           \frac{h^2}{\sigma_h^2} +           \frac{v^2}{\sigma_v^2}          \right]       \right)$
For the purposes of this wiki, this distribution will be called the **Orthogonal Elliptical Distribution**. It is obviously a special case of the Hoyt distribution which in turn is a special case of the bivariate normal distribution. 

In order of the model complexity, the following statistical measurements are appropriate:

- Individual Horizontal and Vertical variances
- Elliptic Error Probable

In this case the horizontal and vertical standard deviations could be determined independently from the horizontal and vertical measurements respectively.

## Case 4, Unequal variances and correlated (Hoyt Distribution) 

<figure class="figure-right" style="max-width: 250px" markdown="span">
  ![Hoyt Distribution - Shots dispersed about COI in an elliptical pattern which has its major axis at an angle to the coordinate axes.](images/Hoyt.jpg)
  <figcaption markdown="span">Hoyt Distribution - Shots dispersed about COI in an elliptical pattern which has its major axis at an angle to the coordinate axes.</figcaption>
</figure>


Given:  
 
1. $\sigma_h \neq \sigma_v$  

1. $\rho \neq 0$  

1. The {*h,v*} position of each shot must be known.
then the mathematical formula for the dispersion distribution would be the Hoyt distribution with no simplifications:  

    $f(h,v) =       \frac{1}{2 \pi \sigma_h \sigma_v \sqrt{1-\rho^2}}       \exp\left(         -\frac{1}{2(1-\rho^2)}\left[           \frac{h^2}{\sigma_h^2} +           \frac{v^2}{\sigma_v^2} -           \frac{2\rho h v}{\sigma_h \sigma_v}         \right]       \right)$

Shot groups would be elliptical or egg-shaped if either the horizontal range or vertical range were large. The following statistical measurements are appropriate:

- Elliptic Error Probable

# Related topics

See also the following topics which are closely related:

- [Error Propagation](error-propagation.md) - A basic discussion of how errors propagate when making measurements.
- [Stringing](stringing.md) - Definition of stringing and how it can be handled

# References
