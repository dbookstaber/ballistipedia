# Leslie, 1993, *Is "Group Size" the Best Measure of Accuracy?*
[*Is "Group Size" the Best Measure of Accuracy?*, John "Jack" E. Leslie III, 1993](media/Is_Group_Size_the_Best_Measure_of_Accuracy_by_J.E._Leslie_III.pdf).  

**ABSTRACT:**

Compares the following measures as a function of the number of shots per group.

- Extreme Spread: Maximum distance between any two shots in group.
- Figure of Merit (FoM): Average of the maximum horizontal group spread and the maximum vertical group spread.  This uses only 2-4 data points depending on the group.  Like Diagonal, FoM becomes more efficient than Extreme Spread for larger group sizes.
- Mean Radius: Average distance to center of group for all shots.
- Radial Standard Deviation: Sqrt (Horizontal Variance + Vertical Variance).

Found military using RSD and Mean Radius as early a 1918.

His Monte Carlo analysis shows sample RSD to be most efficient predictor of precision, followed closely by Mean Radius.  I.e., they can distinguish between loads of different inherent precision more accurately and using fewer sample shots than the other measures.

**Ballistipedia Notes:** 
<ol>
<li> In discussing Extreme spread measurement Leslie makes the following statement: *Also, by only using data from two shots within the group, it ignores the data represented by the other, more likely to be repeated, shots.*

<p>This is the right notion, but not quite correct from the point of view of statistics. From a statistical point of view there is the sample size, the number of shots in a group, and an "effective" sample size also known as the degrees of freedom. On average the ES will increase as the number of shots in a group increases. But with each increase in average ES, the next shot is less likely to increase the size of the ES. Thus subsequent shots don't increase the degrees of freedom by 1, but only by a fraction, and that fraction gets smaller and smaller as the number of shots increases. </p> 

<p>The point  would apply to the Diagonal and the FOM measurements as well.</p>
<li> Leslie, like Grubbs, estimates MR by sampling the mean of radii.  This is less efficient than using [the MR computed from the Rayleigh parameter fitted to the sample](closed-form-precision.md#mean-radius-mr).  The latter process is equally and maximally efficient for all invariant measures that are products of the Rayleigh distribution parameter $\Re$.
<li> Leslie notes "The RSD was the most accurate measure I examined for determining the tightest grouping load."
> But of course, since the simulation is based on the Rayleigh model for which the RSD measurement is essentially the fitted parameter.
<li> Leslie compares "load differences" to get a notion of the relative performance of the statistics. Take the analysis with a grain of salt. There are additional considerations to such an analysis. First, different "loads" would probably have different COI's as well as different dispersions. Second, since 5-shot groups are about optimal for ES, then 20 total shots should be shot as four 5-shot groups since the average of 4 groups would have a much better precision than the ES of one 20-shot group. Third, this analysis ignores fliers. 
</ol>
