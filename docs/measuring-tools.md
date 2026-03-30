# Measuring Tools

Tools and practical methods for measuring and analyzing precision.

# [Ballistipedia Spreadsheet](https://docs.google.com/spreadsheets/d/1i_trin4mHuTJI4HnAnPiVwJHMc0m0X09E7Ya7_aWGFQ/)
Paste your data; get [the statistics](https://colab.research.google.com/drive/1FN4Nq-14N-JhYEtwYcSydJgJG5o0io6t#scrollTo=IVslsqfF-WGL).  [Also available as an Excel spreadsheet](https://github.com/dbookstaber/ballistipedia/blob/main/Ballistipedia.xlsx).

# [2-Shot Method](prior-art.md#danielson-2005-testing-loads)
If you're willing to sacrifice statistical efficiency, [Brent Danielson](prior-art.md#danielson-2005-testing-loads) noted that you can get away with a single measurement for successive 2-shot groups, instead of the coordinate (*x*, *y*) measurements of each shot required for the [maximally efficient estimator](closed-form-precision.md#rayleigh-estimates).  I.e.,

1. Fire two shots at a single point of aim.
1. Measure their center-to-center distance.

This produces two sample radii, both equal to half the measured distance.

The benefit of this approach is that it only requires one measurement in one dimension for every two shots.  The drawback is that, to get the same statistical confidence, this requires almost double the number of shots as if measuring the coordinates of all the shots in a single group.  (With *n* shots split over *g* groups, [the statistical formulas](closed-form-precision.md#confidence-intervals) show that confidence is an increasing function of (*n-g*), so going from 1 group to *n*/2 groups requires 2*n*-1 shots.)

Calculation of sigma from Danielson's sample data, as well as confidence intervals and hypothesis testing, are shown in [DanielsonExample.xlsx](media/DanielsonExample.xlsx).

# [OnTarget](http://ontargetshooting.com/)
Jeffrey Block's [OnTarget Precision Calculator](http://ontargetshooting.com/) is the most convenient package for converting a target image into data points for analysis. It accounts for scale and distance and automatically calculates [Mean Radius](describing-precision.md#mean-radius-mr) (called "Average to Center" in the software) and [Extreme Spread](describing-precision.md#extreme-spread) (called "Max Spread").

The more expensive [Target Data System](http://ontargetshooting.com/tds/) can automatically identify and aggregate shots on scans of its specially-coded targets.

# [Taran](http://taran.ptosis.ch/taran.html)
[Taran](http://taran.ptosis.ch/taran.html) (target analysis and shooting precision calculator) is a free online application to upload a target image, mark the points of impact, and download the coordinates of the points. Among others, it also calculates the Rayleigh CEP.

# [shotGroups Analysis Package](https://github.com/dwoll/shotGroups)
The free [shotGroups](http://cran.fhcrc.org/web/packages/shotGroups/index.html) package for the open-source [statistical environment R](http://www.r-project.org/) provides functions to analyze target groups with respect to their shape, location (accuracy) and spread (precision). Among others, it provides implementions for many [CEP estimators](circular-error-probable.md) and descriptive [precision measures](describing-precision.md). The package works with point data exported from [OnTarget](#ontarget) or [Taran](http://taran.ptosis.ch/taran.html) and includes functions to plot the group with precision indicators like the bounding box, maximum spread or minimum covering circle.

**The main functionality of the package is also available as a set of web applications that do not require installing R or using R syntax.:**
- [Comprehensive shot group analysis](http://dwoll.shinyapps.io/shotGroupsApp)
- [Absolute <-> angular size conversion](http://dwoll.shinyapps.io/shotGroupsAngular)
- [Region <-> hit probability calculations](http://dwoll.shinyapps.io/shotGroupsHitProb)
- [Estimate Rayleigh *σ* parameter from range statistics](http://dwoll.shinyapps.io/shotGroupsRangeStat)

For more information, see the [package description](http://cran.fhcrc.org/web/packages/shotGroups/vignettes/shotGroups.pdf) including a walk-through with sample diagrams and the [complete manual for all functions](http://cran.fhcrc.org/web/packages/shotGroups/shotGroups.pdf). After having installed [R](http://cran.fhcrc.org/bin/windows/base/) and [RStudio](http://www.rstudio.com/ide/download/desktop), open RStudio and install shotGroups by running: `install.packages("shotGroups")` . For a first introduction to R, see:

- http://tryr.codeschool.com/
- http://www.statmethods.net/

# Apps

- [TargetScan](http://www.targetshootingapp.com/) (iOS) computes the unbiased estimate of Mean Radius for supported targets.

# Spreadsheet Analysis
Given a target data set, whether compiled using the [2-Shot Method](#2-shot-method) or [OnTarget](#ontarget), [Closed Form Precision](closed-form-precision.md) analysis can be performed using standard spreadsheet functions.  See, for example [CCI 40gr HV 100yd.xlsx](media/CCI%2040gr%20HV%20100yd.xlsx) or any of the other workbooks linked in the Examples.

[RangeStatisticEstimation.xls](media/RangeStatisticEstimation.xls) is spreadsheet for calculating the statistical significance of [Range Statistics](range-statistics.md) estimates.
