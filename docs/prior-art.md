# [Danielson, 2005, *Testing loads*](http://www.public.iastate.edu/~jessie/PPB/Stats/Testing%20loads.htm)
[Brent J. Danielson](http://www.public.iastate.edu/~jessie/PPB/Stats/Testing%20loads.htm) suggests that a practical way to assess and compare precision is to shoot many 2-shot groups and measure the extreme spread of each.

- When one simply wants to assess whether one sample is more precise than another that probability is given by the one-tailed T-test for two samples with unequal variance -- i.e., the spreadsheet function `=1-TTEST({Sample1},{Sample2},1,3)`.
- Using the same data it is possible to determine the precision parameter.  It is also noteworthy that the square of the 2-tailed T-Test gives the exact confidence range for which the precision parameters of two samples do not overlap -- i.e., the probability that two samples have different precision parameters is `=POWER(1-TTEST({Sample1},{Sample2},2,3), 2)`.  This is all illustrated using Danielson's own data set in [DanielsonExample.xlsx](media/DanielsonExample.xlsx).

# Gammon, 2017, *Shot Group Statistics for Small Arms Applications*
[*Shot Group Statistics for Small Arms Applications*](media/Gammon_2017_Shot_Group_Statistics.pdf) gives distributions from Monte Carlo simulation and shows numerous examples of statistical inference from measurements of extreme spread, mean radius, and radial standard deviation.

# Grubbs, 1964, *Statistical Measures of Accuracy for Riflemen and Missile Engineers*
[*Statistical Measures of Accuracy for Riflemen and Missile Engineers*, Frank E. Grubbs, 1964](media/Statistical%20Measures%20for%20Riflemen%20and%20Missile%20Engineers%20-%20Grubbs%201964.pdf).

# [Hogema, 2005, *Shot group statistics*](http://sport-shooters-ots.com/various/ballist/sgs/sgs.htm)
[Jeroen Hogema](http://sport-shooters-ots.com/various/ballist/sgs/sgs.htm) provides an accessible proof of the equivalence between the symmetric bivariate normal and Rayleigh distributions.  He provides extensive examples, simulations, and applications to scoring and load selection, and begins to address the problem of estimating the Rayleigh parameter.

# [Hogema, 2006, *Measuring Precision*](http://sport-shooters-ots.com/various/ballist/precision/Measuring_precision.htm)
[*Picking the most precise ammo, probably.*](http://sport-shooters-ots.com/various/ballist/precision/Measuring_precision.htm)  Jeroen Hogema:

- Reproduces Leslie’s 1993 results.
- Confirms that for radius measures it is preferable to incorporate all data at once, not to break them into separate groups.
- For FOM and ES it is best to generate many groups so as to preserve more data points.
- Looks at T-tests for significance and shows very large groups are needed to detect statistically meaningful differences.

# Leslie, 1993, *Is "Group Size" the Best Measure of Accuracy?*
[*Is "Group Size" the Best Measure of Accuracy?*, John "Jack" E. Leslie III, 1993](media/Is_Group_Size_the_Best_Measure_of_Accuracy_by_J.E._Leslie_III.pdf).  Notes:

- Extreme Spread: Maximum distance between any two shots in group.  Note that this effectively only uses two data points.
- Figure of Merit (FoM): Average of the maximum horizontal group spread and the maximum vertical group spread.  This uses only 2-4 data points depending on the group.  Like Diagonal, FoM becomes more efficient than Extreme Spread for larger group sizes.
- Mean Radius: Average distance to center of group for all shots.
- Radial Standard Deviation: Sqrt (Horizontal Variance + Vertical Variance).
- Found military using RSD and Mean Radius as early as 1918.

His Monte Carlo analysis shows sample RSD to be most efficient predictor of precision, followed closely by Mean Radius.  I.e., they can distinguish between loads of different inherent precision more accurately and using fewer sample shots than the other measures.

Note that, like Grubbs, Leslie estimates MR by sampling the mean of radii.  This is less efficient than using the Rayleigh estimator on the radii, and then [computing MR based on the sample Rayleigh parameter](closed-form-precision.md#mean-radius-mr).  The latter process is equally and maximally efficient for all invariant measures that are products of the Rayleigh parameter σ.

# [Kolbe, 2010, *Group Statistics*](http://www.geoffrey-kolbe.com/articles/rimfire_accuracy/group_statistics.htm)
Attributing the work to [Sitton et. al., 1990](media/Sitton%201990.pdf), [Geoffrey Kolbe](http://www.geoffrey-kolbe.com/articles/rimfire_accuracy/group_statistics.htm) walks through applied statistics to show how many groups of how many shots are required to estimate Extreme Spread to within 10% of the true value with 90% confidence.  Using the values from [Grubb's](prior-art.md#grubbs-1964-statistical-measures-of-accuracy-for-riflemen-and-missile-engineers) thousand-iteration simulations he determines that 7-shot groups are the most efficient for estimating Extreme Spread.

Running [the same analysis](range-statistics.md#estimation) with our [million-iteration simulation values](media/Sigma1RangeStatistics.xls) reveals that 6-shot groups are actually optimal, and not significantly more so than 5-shot groups.

# [Molon, 2006, *The Trouble With 3-Shot Groups*](http://www.ar15.com/mobile/topic.html?b3&f=118&t=279218) =
[Through an extended forum thread](http://www.ar15.com/mobile/topic.html?b=3&f=118&t=279218) Molon offers intuitive explanations and illustrations of the problems with Extreme Spread samples.

*We have not been able to identify the real person behind that now-inactive user account.  If anyone knows please contact us so we can give well-deserved credit!*

# Taylor & Grubbs, 1975, *Approximate Probability Distributions for the Extreme Spread*
[*Approximate Probability Distributions for the Extreme Spread*, Taylor & Grubbs, 1975](media/Approximate%20Distributions%20for%20Extreme%20Spread%20-%20Taylor%201975.pdf): notes that there is no closed-form expression for [Extreme Spread](describing-precision.md#extreme-spread).  Uses Monte Carlo simulations of 10,000 iterations to estimate the first four moments of the distribution of extreme spread for shot groups of size 3 to 34.  Checks fit against the Chi, LogNormal, and Weibull distributions.

Note that in [Sigma1RangeStatistics.xls](media/Sigma1RangeStatistics.xls) we have generated quantiles and the first four moments for the Extreme Spread using 1,000,000 iteration simulations for groups of size 2 to 100.
