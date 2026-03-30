# Terms
Overall the definitions of terms, in and out of the glossary, strives to correspond to typical usage, in the shooting sports, in mathematics, or in statistics. Any term used with an atypical meaning will be so noted in the glossary.  If a term is specifically defined in this glossary, then the desire is to use that term consistently in the wiki. 

; Accuracy
> The accuracy of a weapon is the distance between the point of aim (POA) and point of impact (POI) for the projectile. Smaller differences reflect greater accuracy.  See [What is Precision?](what-is-precision.md#precision-vs-accuracy)
> > **Shooters** colloquially interchange *accuracy* with *precision*, or conflate the two concepts.

; [Bivariate Distribution](http://en.wikipedia.org/wiki/Joint_probability_distribution)

;c-t-c
> Center-to-center

; [Circular Error Probable](circular-error-probable.md) (CEP)
> CEP(*p*), for *p* ∈ [0, 1), is the expected radius of the smallest circle, whose center is the COI, that covers proportion *p* of the shot group. When *p* is not indicated it is assumed to be 50%, so by default CEP is the estimated true median shot radius.

; [Chi-squared ($\chi^{2}$) distribution](http://en.wikipedia.org/wiki/Chi-squared_distribution)

; [Chronograph](http://en.wikipedia.org/wiki/Gun_chronograph)
> An instrument which measures projectile velocity
> > **Shooters:** Use the term as defined.
> > **Scientists:** Use the term for a [timing device](http://en.wikipedia.org/wiki/Chronograph)

; [Closed form](http://en.wikipedia.org/wiki/Closed-form_expression)

; COI
> Center of Impact. The COI is the point on the target defined by the mean horizontal and mean vertical measurements of the sample of shots.

; Confidence Interval
> A range around an **Estimate** associated with a **Confidence Level *K***.  The correct way to read these values is to say, "If we repeatedly sampled new data using the same method and computed this Estimate each time, we would expect the estimated true parameter value to fall within the Confidence Interval *K%* of the time."  Therefore, the smaller the Confidence Interval the higher the precision of our Estimate.
 
; Confidence Level *K*
> The probability associated with a **Confidence Interval**.  ***K*** ∈ (0, 1), and we expect that in ***n*** estimates, the true parameter value will fall within the **Confidence Interval** ***n / K*** times.

; [Covering Circle Radius (CCR)](describing-precision.md#covering-circle-radius)

; [Degrees of freedom](http://en.wikipedia.org/wiki/Degrees_of_freedom_(statistics))

; Dispersion
> A term to denote that shots are spread around the *Center Of Impact* (COI) without regard to any particular statistical model or measure. Shot groups that have a small spread around the COI have low dispersion, and shots groups that have a large spread around the COI have a high dispersion. High precision implies low dispersion, thus in a mathematical sense precision and dispersion are inversely related.

; Estimator
> A formula or algorithm for estimating the true value of a population parameter from samples of the population.

; [Extreme Spread](describing-precision.md#extreme-spread)
> When discussing precision, and throughout this site, this refers to *Group Size*, a [Range Statistic](range-statistics.md) defined as the maximum distance between any two shots on a target.
> > **Mathematicians:** Also call the target measure the *bivariate range*.
> > **Shooters:** In shooting sports the availability of cheap and accurate chronographs has proliferated their use. Shooters (the intended audience) typically use *extreme spread* " to mean the difference between the highest and lowest velocities measured by a chronograph. To mathematicians this difference would be the *range* of the measurement. Typically shooters will call the target measure *group size*.

; [Figure of Merit (FOM)](describing-precision.md#figure-of-merit)

; Gaussian Distribution
> Synonym for **Normal Distribution** which is the preferred term.

; Group size
> Center-to-center distances between the two widest shots on a sample target.  Because it is measured from center points it is independent of the caliber of the projectile.  In practice this measurement is typically produced from the outside edge of projectile holes in a target and then the projectile caliber is subtracted.
> > **Mathematicians:** Call this measure the *bivariate range* or the *extreme spread*.

; Mean 

- the mathematical average.
- The mean of a sample of *n* values from the population of values (i.e. the set of all values), $\lbrace x_1, x_2, x_3, ..., x_n\rbrace$, is defined as $\bar{x} \equiv \frac{1}{n} \sum_{i=1}^n x_i$. In the limit as *n* approaches infinity, then $\bar{x}$ approaches $\mu_x$.
- The mean value for the population, as opposed to the mean of a sample, is denoted by $\mu_x$.

; Mean Diameter (MD)
> Twice the **Mean Radius**

; [Mean Radius (MR)](describing-precision.md#mean-radius-mr)
> Mean value of the **Radius** of shots on target.

; Median
> The median is the value for which half the measures are smaller and half are larger. For example for a discrete set of five values ranked from samllest to largest, then the 3^rd^ value would be the median. For a probability distribution function the median would be the value at the 50^th^ percentile.

; [Minute of Angle (MOA)](angular-size.md)
> Synonym for Minute of Arc, i.e., one arc minute

; Mode
> For a continuous distribution, the mode is the peak value. For a sample of values Mode is the value which occurs most frequently (typically from a histogram of the data).

; [Mil](describing-precision.md#units)
> One thousandth.

; [Milliradian](angular-size.md), a.k.a. mrad or milrad

; [Normal distribution](http://en.wikipedia.org/wiki/Normal_distribution)

; POA
> Point of Aim

; POI
> Point of Impact

; [Probability density function (PDF)](http://en.wikipedia.org/wiki/Probability_density_function)

; Precision
> A term to denote how shots are spread around the *Center of Impact* (CoI), generally in reference to a particular statistical model or measure. Shot groups that have a small spread around the CoI have high precision, and shots groups that have a large spread around the CoI have a low precision. Thus low precision implies high dispersion. In this mathematical sense precision and dispersion are inversely related.
> > **Shooters** colloquially interchange *accuracy* with *precision*, or conflate the two concepts.

; [Radial Standard Deviation (RSD)](describing-precision.md#radial-standard-deviation-rsd)

; [Rayleigh Distribution](http://en.wikipedia.org/wiki/Rayleigh_distribution)

; Sample
> *Sample* is used in two different contexts, so careful reading may be required. First a sample could mean the set of measurements made on the shots on a single target. Second it could mean the set of measurements made on a set of such targets.

; [Shooter's MOA (SMOA)](angular-size.md)
> Angular measure defined as one inch at a hundred yards.

; Target
> In this wiki target has the implicit notion of one or more shots taken at the same POA. There are of course *paper targets* with multiple bulls-eyes per paper sheet. For such a paper target each bulls-eye would be a different target.

# Mathematical Notation

| Formula |  |  |  |
| --- | --- | --- | --- |
| $\bar{ }$ as $\bar{h}$ | "Bar" | Denotes a sample average |  |
| $\hat{ }$ as $\hat{\sigma}$ | "Hat" | Denotes an statistical estimator. For example $\sigma$ would be calculated directly from the set of data in the usual way, but $\hat{\sigma}$ would calculated some other unusual way. |  |
| $\alpha$ | Probability of Type I error | A Type I error, also known as an error of the first kind, occurs when the null hypothesis (H~0~) is true, but is rejected. |  |
| $\beta$ | Probability of Type II error | A type II error, also known as an error of the second kind, occurs when the null hypothesis is false, but erroneously fails to be rejected. |  |
| $\mu$ | Population Mean | In essence the population mean is a theoretical abstraction which is unknowable. (eg. Mean position for all possible shots as opposed to a sample of *n* shots.) | For a set of *n* values, $\lbrace x_1, x_2, x_3, ..., x_n\rbrace$ |
| σ | Greek letter sigma, [Population Standard Deviation](http://en.wikipedia.org/wiki/Standard_deviation) | The true population's standard deviation. In essence the population's standard deviation is a theoretical abstraction which is unknowable. (eg. Standard deviation for all possible shots as opposed to a sample of *n* shots.) | For a set of *n* values, $\lbrace x_1, x_2, x_3, ..., x_n\rbrace$ |
| $\rho$ | [Pearson correlation coefficient](http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient) | Measures the amount of liner correlation between the two variables. | :$\rho = \rho_{hv} =\frac{\sum ^n _{i=1}(h_i - \bar{h})(v_i - \bar{v})}{\sqrt{\sum ^n _{i=1}(h_i - \bar{h})^2} \sqrt{\sum ^n _{i=1}(v_i - \bar{v})^2}}$ |
| $c_B(n)$ | [ Bessel Correction Factor](closed-form-precision.md#bessel-correction-factor) | Multiplicative factor which adjusts the sample variance to use (n-1) instead of (n) which reduces the bias of the estimate. | $c_B(n) = \frac{n}{n-1}$ |
| $c_G(n)$ | [ Gaussian Correction Factor](closed-form-precision.md#gaussian-correction-factor) | Multiplicative factor which adjusts the standard deviation to correct bias due to taking the square root of the variance. Thus the Bessel correction corrects the variance itself, but taking the square root introduces a different bias. | $\frac{1}{c_{G}(n)} = \sqrt{\frac{2}{n-1}}\,\frac{\Gamma\left(\frac{n}{2}\right)}{\Gamma\left(\frac{n-1}{2}\right)} \, \approx \, 1 - \frac{1}{4n} - \frac{7}{32n^2} - \frac{19}{128n^3}$ |
| $c_R(n)$ | [Rayleigh Correction Factor](closed-form-precision.md#rayleigh-correction-factor) | Multiplicative factor which adjusts the Rayleigh shape parameter to reduce the bias of the estimate. | $\frac {4^n n!(n-1)!\sqrt{n}} {(2n)!\sqrt{\pi}}$ |
| $h$ | Horizontal position or axis of a target. | Typically a target would be mounted vertically and the line of fire would be horizontal. Thus the horizontal axis of the target would be synonymous with the *X* axis when using X-Y Cartesian coordinates for the target. | $\bar{h} = \frac{1}{n} \sum_{i=1}^n h_i$ |
| $\mathcal{N}(\mu,\,\sigma^2)$ | Normal Distribution | This denotes the normal distribution with mean, $\mu$, and variance $\sigma^2$. |  |
| $r$ | Radius | Here this almost always refers to the distance of a shot on target from the center, or sample center, of a target group. |  |
| $\bar{r}$ | Mean Radius | Mean radius of a sample of *n* shots. | $\bar{r} = \frac{1}{n} \sum_{i=1}^n r_i$ |
| $r_i$ | Radius of *i*^th^ shot | Unless otherwise stated it is implicit that the polar coordinate center is at the mean POI. |  |
| %RSD | [% Relative Standard Deviation](http://en.wikipedia.org/wiki/Relative_standard_deviation) | The relative standard deviation which is a ratio of the standard deviation to the mean as a percentage. | %RSD = $100 \frac{\sigma_x}{\bar{x}}$ |
| $s$ | [Sample Standard deviation](http://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation) | An experimental estimate of the population's standard deviation which is calculated from a of a set of *n* values, $\lbrace x_1, x_2, x_3, ..., x_n\rbrace$. In the limit as *n* approaches $\infty$, then *s* approaches the population standard deviation, *σ*. | $s = \sqrt{\frac{1}{n-1} \sum_{i=1}^n (x_i - \bar{x})^2}$, |
| $s^2$ | [Sample Variance](http://en.wikipedia.org/wiki/Variance#Sample_variance) | The square of the sample standard deviation, typically as done with Bessel's correction of (n-1) instead of (n) to get a more unbiased estimator. |  |
| $v$ | Vertical position or axis of a target. | Typically a target would be mounted vertically and the line of fire would be horizontal. Thus the vertical axis of the target would be synonymous with the *Y* axis when using X-Y Cartesian coordinates for the target. | $\bar{v} = \frac{1}{n} \sum_{i=1}^n v_i$ |
| $x$ | variable $x$ | probably a undefined variable. *X* will be used for the X-axis. |  |
| $X$ | X-axis | $X_i$ will be used for X position of the *i*^th^ shot. |  |
| $y$ | variable $y$ | probably a undefined variable. $Y$ will be used for the Y-axis. |  |
| $Y$ | Y-axis | $Y_i$ will be used for Y position of the *i*^th^ shot. |  |
