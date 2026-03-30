# References

**[Prior Art](prior-art.md)** details previous work on the problem of estimating shooting statistics.

**[CEP literature](cep-literature.md)** focuses on the broader body of work related to characterizing Circular Error Probable, which is applicable not only to ballistics but also to fields like navigation and signal processing.

Following is a complete list of useful References and Prior Art:

- Bookstaber, David (2014).  [**Understanding Rifle Precision**](https://web.archive.org/web/20171120094420/http://www.thetruthaboutguns.com/2014/12/daniel-zimmerman/understanding-rifle-precision/).

- Danielson, Brent J. (2005).  [**Testing Loads** – *detailed in Prior Art*](prior-art.md#danielson-2005-testing-loads).

- Gammon, W. J. (2017), [**Shot Group Statistics for Small Arms Applications** – ''detailed in Prior Art](prior-art.md#gammon-2017-shot-group-statistics-for-small-arms-applications).

- Grubbs, Frank E. (1964).  [**Statistical Measures of Accuracy for Riflemen and Missile Engineers** – *detailed in Prior Art*](prior-art.md#grubbs-1964-statistical-measures-of-accuracy-for-riflemen-and-missile-engineers).

- Hogema, Jeroen (2005).  [**Shot group statistics** – *detailed in Prior Art*](prior-art.md#hogema-2005-shot-group-statistics).

- Hogema, Jeroen (2006).  [**Measuring Precision** – *detailed in Prior Art*](prior-art.md#hogema-2006-measuring-precision).

- Hornady Podcast (2022).  [**Episode 050 – Your Groups Are Too Small**](https://www.youtube.com/watch?v=QwumAGRmz2I).  [Summary here](https://www.snipershide.com/shooting/threads/ill-post-this-here-hornadys-podcast-50-i-thought-it-was-one-of-their-best-but-some-reloaders-might-not-like-what-they-see.7152297/post-10620383).

- Kolbe, Geoffrey (2010).  [**Group Statistics** – *detailed in Prior Art*](prior-art.md#kolbe-2010-group-statistics).

- [mailto:jleslieiii@icloud.com Leslie, John E. III] (1993).  [**Is "Group Size" the Best Measure of Accuracy?** – *detailed in Prior Art*](prior-art.md#leslie-1993-is-group-size-the-best-measure-of-accuracy).

- MacDonald, Adam (2017). [**Thinking Statistically**](https://www.autotrickler.com/blog/thinking-statistically).

- Molon (2006). [**The Trouble With 3-Shot Groups** – *detailed in Prior Art*](prior-art.md#molon-2006-the-trouble-with-3-shot-groups).

- Precision Rifle Blog (2020). [**Statistics for Shooters**](https://precisionrifleblog.com/2020/12/12/measuring-group-size-statistics-for-shooters/).

- Rifleslinger (2014). [**On Zeroing**](http://artoftherifleblog.com/on-zeroing/2014/02/on-zeroing.html).

- Saleh, A. K. Md. Ehsanes (1967). [**Determination of the Exact Optimum Order Statistics for Estimating the Parameters of the Exponential Distribution from Censored Samples**](media/Order_Statistics_of_Exponential_Distribution_in_Censored_Samples.pdf). Technometrics 9, no. 2.

- Sarhan, A. E., Greenberg, B. G., & Ogawa, J. (1963). [**Simplified Estimates for the Exponential Distribution**](media/Simplified_Estimates_for_the_Exponential_Distribution.pdf). The Annals of Mathematical Statistics, 34(1), 102–116.

- Siddiqui, M. M. (1961). [**Some Problems Connected With Rayleigh Distributions**](media/Some%20Problems%20Connected%20With%20Rayleigh%20Distributions%20-%20Siddiqui%201961.pdf).  The Journal of Research of the National Bureau of Standards, Sec. D: Radio Science, Vol. 68D, No. 9.

- Siddiqui, M. M. (1964). [**Statistical Inference for Rayleigh Distributions**](media/Statistical%20Inference%20for%20Rayleigh%20Distributions%20-%20Siddiqui,%201964.pdf).  The Journal of Research of the National Bureau of Standards, Sec. D: Radio Propagation, Vol. 66D, No. 2.  (*Summarizes and extends Siddiqui, 1961.*)

***Important Note on Siddiqui**: Siddiqui parameterizes the Rayleigh distribution with $\frac{\sigma}{\sqrt{2}}$.  Therefore, should you endeavor to relate Siddiqui's work to that referenced here and in more modern usage, remember that $\sigma_{modern} = \sqrt{2} \sigma_{Siddiqui}$.*

- Strohm, Luke (2013).  [**An Introduction to the Sources of Delivery Error for Direct-Fire Ballistic Projectiles** (ARL-TR-6494)](https://apps.dtic.mil/sti/tr/pdf/ADA588846.pdf).

- Taylor, M. S. & Grubbs, Frank E. (1975).  [**Approximate Probability Distributions for the Extreme Spread** – *detailed in Prior Art*](prior-art.md#taylor-grubbs-1975-approximate-probability-distributions-for-the-extreme-spread).

- Triplett, Ben (2019). [**Rifle Ammunition Load Workup**](http://www.bisonops.com/2019/08/17/rifle-ammunition-load-workup).


# Reference Data

- [Confidence Interval Convergence.xlsx](media/Confidence%20Interval%20Convergence.xlsx): Shows how precision confidence intervals shrink as sample size increases.

- [Sigma1RangeStatistics.xls](media/Sigma1RangeStatistics.xls): Simulated median, 50%, 80%, and 95% quantiles, plus first four sample moments, for shot groups containing 2 to 100 shots, of: Extreme Spread, Diagonal, Figure of Merit.

- [SymmetricBivariateSigma1.xls](media/SymmetricBivariateSigma1.xls): Monte Carlo simulation results validating the [Closed Form Precision](closed-form-precision.md) math.

[BallisticSimulations.ipynb](https://github.com/dbookstaber/ballistipedia/blob/main/BallisticSimulations.ipynb) is a Jupyter notebook containing extensive illustrations and validation of the math used throughout this site.
