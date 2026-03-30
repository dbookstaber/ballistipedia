# 300BLK Subsonic 20-shot 100-yard Example

- **Shots:** 20
- **Distance:** 100 yards
- **Ammunition:** Remington Express 220gr Subsonic
- **Barrel:** Noveske 16" 1:7 stainless
- **Action:** Noveske AR-15 with pistol-length gas tube
- **Conditions:** 35F, 40%RH, slight wind
- **Precision:** Raw precision is 1.4MOA.  Statistically adjusted for muzzle velocity variance precision is 1.0MOA.  Rank-adjusted for measured muzzle velocity precision is 0.70MOA.

This example shows how critical it is to minimize muzzle velocity variance for subsonic rifle precision.  Average muzzle velocity was 960fps, but the values observed for this sample ran from 883 to 1013 fps, which represents a 14% extreme spread.  At 100 yards that spread equates to over 5 inches of drop!

[300BLK 220gr Factory Subsonic 100yd.xlsx](media/300BLK%20220gr%20Factory%20Subsonic%20100yd.xlsx) shows the raw target data as well as two methods of adjustment.  It contains 5 sheets:

1. Standard precision analysis of the Raw Shot Data
1. Analysis of the chronograph data to determine vertical variance attributable to the dispersion in muzzle velocity
1. **Statistically Adjusted** Shot Data, correcting each shot uniformly for muzzle velocity variance
1. **Rank Adjusted** Shot Data, correcting each shot by correlating each muzzle velocity measurement with the vertical point of impact
1. [Chart of both the raw and rank-adjusted shot data](#target-data)

The **Statistical Adjustment** based on variances is a *lower bound* on the inherent precision of this system after controlling for the spread in muzzle velocity.  The **Rank Adjustment** is an *upper bound*, because it assumes that the vertical point of impact is 100% correlated with muzzle velocity, which is probably not true.  In order to determine the *true precision controlling for muzzle velocity* we would have to keep track of the velocity associated with each impact, which is a more tedious exercise.

Nevertheless, the Rank Adjustment is probably close to the true value because it reduces the vertical variance to the same magnitude as the horizontal variance, and also because the resulting sigma is in line with [typical values of light factory rifle precision](closed-form-precision.md#typical-values-of).

# Target Data

![300BLK Subsonic Target w MV Adjustment.png](images/300BLK%20Subsonic%20Target%20w%20MV%20Adjustment.png)
