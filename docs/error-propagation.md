# Error Propagation

## Basic Assumptions

In general shooting has a number of different error sources. The "basic" assumption is that the error sources are each normally distributed and all are independent of each other. Thus for example the vertical error component would be represented by the equation:  

    $\sigma_v^2 = \sigma_{v,Weapon}^2 + \sigma_{v,Ammunition}^2 + \sigma_{v,Shooter}^2 + \sigma_{v,Measure}^2$  

and there would similar equations for $\sigma_h^2$, $\sigma_{RSD}^2$ and other measurements.

Notes:

1. The error sources can't be measured independently. Only the total error is observable. This is a very important consideration since the total error thus limits the precision which any of the individual error factors can be determined.
1. If any one of the error sources is more than about three times the sum of the others, then that error source essentially controls the overall error of the measurement.
1. These error sources could be decomposed further into "finer" errors. For instance consider the multitude of variables that a handloader controls when loading a cartridge.
1. The fact that a very small amount of data is typically collected by a shooter greatly limits the overall precision of the mathematical analysis. In statistics this is known as the "small sample" problem. For example to get "good" estimates for the normal distribution parameters $\bar{x}$ and particularly $\sigma$, 30 measurements are desired!
1. In general for a given sample size of *n* measurements the mean, $\bar{x}$, is determined relatively (i.e as a percentage) much more precisely than the standard deviation, $\sigma$ (or other dispersion factor).
1. For practical purposes the Measurement error is assumed to be insignificant.

## Random Error Scenario

## Random Defect Scenario

# Error Handling Examples


## Examples Using System Errors Of The Populations

Let's assume the following for an Extreme Spread measurements for the **population values**:  


$\sigma_{ES}^2 = \sigma_{ES,Weapon}^2 + \sigma_{ES,Ammunition}^2 + \sigma_{ES,shooter}^2$  


$\sigma_{ES,Weapon}^2 = 0.30$ inches

$\sigma_{ES,Ammunition(1)}^2 = 0.50$ inches

$\sigma_{ES,Ammunition(2)}^2 = 0.70$ inches

$\sigma_{ES,Shooter}^2 = 0.40$ inches

### What is the expected precision, $\sigma_{ES(1)}$, for type 1 ammo?

Substituting for the values into the general equation:  

 
$\sigma_{ES(1)}^2 = 0.30^2 + 0.50^2 + 0.40^2 = 0.50$  


thus $\sigma_{ES(1)} = \sqrt{0.5} = 0.707$ inches  


### What is the expected precision, $\sigma_{ES(2)}$, for type 2 ammo?

Substituting for the values into the general equation:  


$\sigma_{ES(2)}^2 = 0.30^2 + 0.70^2 + 0.40^2 = 0.740$  


thus $\sigma_{ES(2)} = \sqrt{0.5} = 0.860$ inches  


### From the population values, what is the expected difference between Type(1) and Type(2) ammo?


## Examples Using Experimental System Errors Of The Populations

Given:  


$s_{ES}^2 = s_{ES,Weapon}^2 + s_{ES,Ammunition}^2 + s_{ES,shooter}^2$  


$s_{ES(1)} = 0.707$ inches, i.e. ammunition type 1  


$s_{ES(2)} = 0.860$ inches, i.e. ammunition type 2  


### What is the precision of the difference measurement between Type(1) and Type(2) ammo?

**There is a nasty surprise in this.** 

The precision with which this measurement can be made is very poor.   

$\sigma_{diff}^2 = \sigma_{ES(1)}^2 + \sigma_{ES(2)}^2$ inches  

or:  

$\sigma_{diff} = \sqrt{0.707^2 + 0.860^2} = 1.113$ inches  


The nasty surprise is that while variances can be subtracted to get differences, the standard deviation of the difference is equal to the sum of all of the variances. Thus the trend is that since both measurements have some random error, that the random errors don't negate each other but combine to give an error bigger than any of the individual errors.  

If multiplying or dividing, then it is the relative error (%) which is propagated. Given $\sigma_A = 10$% and $\sigma_B = 18$% then for the product A*B the error would be:  


$\sigma_{A*B} = \sqrt{0.1^2 + 0.18^2} = 20.6%$  


## Error Dominance

### Measurement Error

Measurement error was specifically included in the overall system error to acknowledge that it does indeed exist. There must be some error in locating the (h,v) positions of holes in a target. However in general the measurement error is assumed to very small. Mathematically this can be expressed as:

    $\sigma_{v,Measure}^2 << \sigma_{v,Weapon}^2 + \sigma_{v,Ammunition}^2 + \sigma_{v,Shooter}^2$  


Thus the overall error is dominated by the errors other than the measurement error.

### Large Ammunition Error

Let's again assume the basic equation for the Extreme Spread measurement, but in addition let's assume that:  

  $s_{ES,Ammunition} = 3\sqrt{s_{ES,Weapon}^2 + s_{ES,shooter}^2}$.

**What is the effect of this situation?**

Substituting into the basic equation we can see that:  

  $s_{ES}^2 = s_{ES,Ammunition}^2 + \frac{1}{9}s_{ES,Ammunition}^2 = \frac{10}{9}s_{ES,Ammunition}^2$  


Thus in this case the ammunition error dominates the total error of the measurement. Experiments using different kinds of ammunition would tend to readily show differences between ammunition types. Shooting a few groups of each type should be sufficient.

### Small Ammunition Error
As before let's assume the basic equation for the Extreme Spread measurement, but in let's flip the assumption about ammunition controlling the experimental error around. Now let's assume that:  

  $s_{ES,Ammunition} = \frac{1}{3}\sqrt{s_{ES,Weapon}^2 + s_{ES,shooter}^2}$.

**What is the effect of this situation?**

Now the ammunition contribution is only 1/10th of the total error in variance. In order to precisely measure differences between ammunition types a lot of measurements are going to be needed. Thus a large number of groups for each ammunition types will be needed.
