# 22LR CCI 40gr HV 40-shot 100-yard Example

- **Shots:** 40
- **Distance:** 100 yards
- **Ammunition:** CCI 40gr HV .22LR
- **Barrel:** EAB 18" 1:9, Target chamber, with Gem-Tech Outback II suppressor
- **Action:** PWS T3 toggle bolt
- **Conditions:** 82F, 60%RH, no significant wind
- **Precision:** 0.75MOA, 95% confidence interval (0.65, 0.89)

In this example I shot two strings of 20 shots of CCI 40gr HV .22LR ammunition at 100 yards.  I was collecting data to estimate the ballistic coefficient, so precision wasn't the primary objective.

This example is notable because I had a chronograph on each string, so I have a measure of the dispersion of muzzle velocity.  The raw targets show significantly more spread on the Y axis: vertical standard deviation is .9" vs. .8" horizontal standard deviation.  However when we factor out the vertical spread attributable to muzzle velocity the variance in each axis is practically identical.

The precision analysis, as well as the data and process for adjusting the sample for muzzle velocity dispersion, is documented in [CCI 40gr HV 100yd.xlsx](media/CCI%2040gr%20HV%20100yd.xlsx), which has four sheets:

1. Standard precision analysis of the raw data
1. Analysis of the chronograph data to determine vertical variance attributable to the dispersion in muzzle velocity
1. Precision analysis of the shot data adjusted for muzzle velocity variance
1. [Chart of the adjusted shot data](#combined-adjusted-data)

# Raw Targets
![40grCCI HV 100Yards String1.jpg](images/40grCCI%20HV%20100Yards%20String1.jpg) 
![40grCCI HV 100Yards String2.jpg](images/40grCCI%20HV%20100Yards%20String2.jpg)

# Combined Adjusted Data
The combined, adjusted data look like this.  The estimated CEP is indicated, and it happens to fully encompasses exactly 50% of the sample shots, although two more (for a total of 55%) are technically inside its radius:

![22LR CCI 40gr HV 40shot.png](images/22LR%20CCI%2040gr%20HV%2040shot.png)


---

For a case in which velocity adjustment is applied by rank see [300BLK Subsonic 20-shot 100-yard Example](300blk-subsonic-20-shot-100-yard-example.md).
