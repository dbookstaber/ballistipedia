# Circular Error Probable (CEP)

The Circular Error Probable $CEP(p)$ for $p \in [0,1)$ is the estimated radius of the smallest circle that is expected to cover proportion $p$ of the shot group. Some authors restrict the name "CEP" to the case of $p = 0.5$, and refer to, e.g., $R95$ for $p = 0.95$. We make no such distinction here.

How $CEP(p)$ should be estimated depends on what assumptions are made regarding the distribution of radial errors, i.e., the distribution of miss distances of shots to the point of aim (POA). In turn, the distribution of radial error depends on the bivariate distribution of ***x***- and ***y***-coordinates of the shots. If the ***x***- and ***y***-coordinates of the shots follow a bivariate normal distribution, the radial error around the POA can follow one of several distributions, depending on the cirumstances (Beckmann [1962](cep-literature.md#beckmann1962); [1964](cep-literature.md#beckmann1964)):


<figure class="figure-right" style="max-width: 400px" markdown="span">
  ![Distribution of radial error (red arrow) in different kinds of bivariate normal distribution. POA = point of aim, POI = mean point of impact](images/radialErrorDistributions.jpg)
  <figcaption markdown="span">Distribution of radial error (red arrow) in different kinds of bivariate normal distribution. POA = point of aim, POI = mean point of impact</figcaption>
</figure>


1. **Rayleigh**: When the true center of the coordinates and the POA coincide, the radial error around the POA in a bivariate uncorrelated normal random variable with equal variances follows a [Rayleigh distribution](http://reference.wolfram.com/language/ref/RayleighDistribution.html). This distribution is described in the [Closed Form Precision](closed-form-precision.md) section. In three dimensions (spherical error probable, SEP), the radial error follows a [Maxwell-Boltzmann distribution](http://reference.wolfram.com/language/ref/MaxwellDistribution.html).
1. **Rice**: When the true center of the coordinates and the POA are not identical, the radial error around the POA in a bivariate uncorrelated normal random variable with equal variances follows a [Rice distribution](http://reference.wolfram.com/language/ref/RiceDistribution.html). The probability density function, the cumulative distribution function, and the quantile function are defined in closed form. The Rice distribution reduces to the Rayleigh distribution if the mean coincides with the POA.
1. **Hoyt**: When the true center of the coordinates and the POA coincide, the radius around the POA in a bivariate correlated normal random variable with unequal variances follows a [Hoyt distribution](http://reference.wolfram.com/language/ref/HoytDistribution.html). The probability density function and the cumulative distribution function are defined in closed form, whereas numerical methods are required to find the quantile function. The Hoyt distribution reduces to the Rayleigh distribution if the correlation is 0 and the variances are equal.
1. The **general case** obtains if the true center of the coordinates and the POA are not identical, and the shots have a bivariate correlated normal distribution with unequal variances. The cumulative distribution function of radial error is equal to the integral of the bivariate normal distribution over an offset disc. Numerical methods are required to evaluate the distribution. The resulting distribution reduces to the Rice distribution if the correlation is 0 and the variances are equal. The resulting distribution reduces to the Hoyt distribution if the mean has no offset.

## Systematic Accuracy Bias
Some approaches to estimating CEP conflate the question of precision with the question of accuracy, or "sighting in."

The simpler case only tries to estimate precision, and computes CEP about the sample center.

The general case allows that the point-of-aim is offset from the true center point-of-impact.  In the literature this is referred to as *systematic accuracy bias*.  Including systematic accuracy bias sets the center of the circle to the point of aim, which means the sample center will probably be offset from that and CEP will be correspondingly larger.

## Estimators
Several different methods for estimating $CEP(p)$ have been proposed which are based on the different assumptions about the underlying distribution of coordinates outlined above. See the [CEP literature overview](cep-literature.md) for references and the [shotGroups](measuring-tools.md#shotgroups-analysis-package) package for a free open source implementation:

- The general correlated normal estimator ([DiDonato & Jarnagin, 1961a](cep-literature.md#didonato1961a); [Evans, 1985](cep-literature.md#evans1985)) is based on the assumption of [bivariate normality](precision-models.md#general-bivariate-normal) of the shot coordinates. It allows the ***x***- and ***y***-coordinates to be correlated and have different variances.
    - Without taking systematic bias into account, this estimate can be based on the closed-form solution for the Hoyt distribution of radial error ([Hoyt, 1947](cep-literature.md#hoyt1947); [Paris, 2009](cep-literature.md#paris2009)). The [Krempasky (2003)](cep-literature.md#krempasky2003) estimate is based on a nearly correct closed-form solution for the 50% quantile of the Hoyt distribution. The [Ignani (2010)](cep-literature.md#ignani2010) estimate is based on a polynomial approximation for the 50%, 90%, 95%, and 99% quantiles of the Hoyt distribution.
    - The modified [RAND R-234](cep-literature.md#rand1952) estimator ([Pesapane & Irvine, 1977](cep-literature.md#pesapane1977)) is an early example of CEP and is based on lookup tables for the 50% quantile of the Hoyt distribution. The tables were later cast into an algebraic form that is essentially the Rayleigh estimator with a weighted average of the variances of the de-correlated data to estimate the true standard deviation. It works best for a mostly circular distribution of $(x,y)$-coordinates (aspect ratio of data ellipse $\leq 3$).
    - The Valstar estimate ([Puhek, 1992](cep-literature.md#puhek1992)) for the 50% quantile of the Hoyt distribution differs from the RAND-estimate only for highly elliptical distributions. This estimate does not generalize to three dimensions.
    - Other old, and less relevant approximations to the 50% quantile of the Hoyt distribution include [Bell (1973)](cep-literature.md#bell1973), [Nicholson (1974)](cep-literature.md#nicholson1974) and [Siouris (1993)](cep-literature.md#siouris1993).
- If systematic accuracy bias is taken into account, [numerical integration](cep-literature.md#algos) of the multivariate normal distribution around an offset circle is required for an exact solution. The calculation of the correlated normal estimator is difficult and requires numerical approaches only available in [specialized software](measuring-tools.md#shotgroups-analysis-package).
    - An approximation for the 50% and 90% quantile when there is systematic bias comes from [Shultz (1963)](cep-literature.md#shultz1963), later modified by [Ager (2004)](cep-literature.md#ager2004). The RAND-tables have also been fitted with a regression model to accommodate systematic accuracy bias in the 50% quantile ([Pesapane & Irvine, 1977](cep-literature.md#pesapane1977)). The Valstar estimate ([Puhek, 1992](cep-literature.md#puhek1992)) is similar but differs in its method of correcting for systematic accuracy bias.
- The Grubbs-Pearson estimator ([Grubbs, 1964](cep-literature.md#grubbs1964)) shares its assumptions with the general correlated normal estimator. It is based on the Pearson three-moment central $\chi^{2}$-approximation ([Imhof, 1961](cep-literature.md#imhof1961); [Pearson, 1959](cep-literature.md#pearson1959)) of the cumulative distribution function of radial error in bivariate normal variables. This approach has the advantage that its calculation is much easier than the exact distribution and does not require special software. For $p \geq 0.25$, the approximation to the true cumulative distribution function is very close but can diverge from it for $p < 0.25$ with some distribution shapes.

- The Grubbs-Patnaik estimator ([Grubbs, 1964](cep-literature.md#grubbs1964)) differs from the Grubbs-Pearson estimator insofar as it is based on the Patnaik two-moment central $\chi^{2}$-approximation ([Patnaik, 1949](cep-literature.md#patnaik1949)) of the true cumulative distribution function of radial error. For $p < 0.5$ with some distribution shapes, the approximation can diverge significantly from the true cumulative distribution function.

- The Grubbs-Liu estimate was not proposed by Grubbs but can be constructed following the same principle as his original estimators. It differs from them insofar as it is based on the recent [Liu, Tang, and Zhang (2009)](cep-literature.md#liu2009) four-moment non-central $\chi^{2}$-approximation of the true cumulative distribution function of radial error.

- The [Rayleigh estimator](closed-form-precision.md) uses the Rayleigh quantile function for radial error ([Culpepper, 1978](cep-literature.md#culpepper1978); [Singh, 1992](cep-literature.md#singh1992)). It assumes an uncorrelated bivariate normal process with equal variances and zero mean. Note that this estimator is essentially the same as the RMSE estimator often described in the GPS literature when using centered data for calculating MSE.[^1] [^2][^3] The only difference is that RMSE uses a biased estimate of the true standard deviation.
    - If systematic accuracy bias is taken into account, this estimator becomes the Rice estimator. Note that for small bias, this estimator is similar to the RMSE estimator often described in the GPS literature when using the original, non-centered data for calculating MSE. With large bias however, the RMSE estimator becomes seriously wrong.

- The [Ethridge (1983)](cep-literature.md#ethridge1983) estimator is not based on the assumption of bivariate normality of $(x,y)$-coordinates but uses a robust unbiased estimator for the median radius ([Hogg, 1967](cep-literature.md#hogg1967)). This estimator "assumes that the square root of the radial miss distances follows the logarithmic generalized exponential power distribution." ([Williams, 1997](cep-literature.md#williams1997)). It is only available for $p = 0.5$.

# Comparing CEP estimators

If the true variances of ***x***- and ***y***-coordinates as well as their covariance is known then the closed-form general correlated normal estimator is ideal.

When we are confident in asserting a bivariate normal model for shot dispersion the Grubbs estimators are excellent approximations for reasonable values of ***p*** and ellipticity.

- To date most [comparison studies](cep-literature.md#compstudies) have only used the Grubbs-Patnaik estimator.
- The Grubbs-Pearson estimator has the theoretical advantage over the Grubbs-Patnaik estimator that the approximating distribution matches the true distribution not only in mean and variance but also in skewness.  Both the Grubbs-Pearson and Grubbs-Patnaik estimators are easy to calculate with standard software as long as the central $\chi^{2}$-distribution is available (as it is, for example, in [spreadsheets](measuring-tools.md#spreadsheet-analysis)).
- If systematic accuracy bias is taken into account, the Grubbs-Liu estimator has the theoretical advantage over the Grubbs-Pearson estimator that the approximating distribution matches the true distribution not only in mean, variance, and skewness but also in kurtosis. If systematic accuracy bias is ignored, the Grubbs-Liu estimator is equivalent to the Grubbs-Pearson estimator. Its calculation is less complicated than the exact correlated normal estimator but requires the *non-central* $\chi^{2}$-distribution. This distribution might not be available in general tools like spreadsheets, but it is implemented in all statistics packages.

One shortcoming of the Grubbs estimators is that it is not possible to incorporate the confidence intervals of the variance estimates into the CEP estimate.  This is particularly relevant to small samples where the variance estimates themselves are subject to considerable error.

In the special case where we assume uncorrelated bivariate normal data with equal variances the [Rayleigh estimator](closed-form-precision.md) does provide true confidence intervals, and it is easy to calculate using [spreadsheets](measuring-tools.md#spreadsheet-analysis).

The Ethridge estimator stands out because it does not require bivariate normality of the $(x,y)$-coordinates. It generalizes to three-dimensional data and can accommodate systematic accuracy bias, but it is limited to the 50% CEP.

## Small Samples

For small samples we are more sensitive to which estimator is least biased and most efficient. This question has been studied, e.g., by [Williams (1997)](cep-literature.md#williams1997).

A related question is which estimator is most robust to a very small number of outliers (or [Fliers](fliers.md)) that may result from clear operator error. See the literature overview for more [comparison studies](cep-literature.md#compstudies).

# References


[^1]: [''GPS Accuracy: Lies, Damn Lies, and Statistics'', Frank van Diggelen, GPS World, 1998](http://gpsworld.com/gps-accuracy-lies-damn-lies-and-statistics/)

[^2]: [''Update: GNSS Accuracy: Lies, Damn Lies, and Statistics'', Frank van Diggelen, GPS World, 2007](http://gpsworld.com/gpsgnss-accuracy-lies-damn-lies-and-statistics-1134/)

[^3]: [''The Mathematics of GPS'', Richard B Langley, GPS World, 1991](http://gauss.gge.unb.ca/gpsworld/EarlyInnovationColumns/Innov.1991.07-08.pdf)
