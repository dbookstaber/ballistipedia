# Sighting a Weapon

"Sighting" a weapon involves trying to align its sighting device so that the Point Of Aim (POA) coincides with the Center Of Impact (COI) as measured for some sample of shots on a given target. Let $COI_{\infty}$ represent the COI for an infinite number of possible shots, the population of shots.  The statistical challenge with sighting a weapon is that given any shot dispersion, and every weapon has some, then we can only estimate $COI_{\infty}$ with the COI from the target. It is also implicit that increasing the number of shots used to determine the COI for the target will improve our estimate of the $COI_{\infty}$. 

<br clear=all>

# Sighting Process

As Einstein said, "Everything is relative." So for the purpose of this wiki we need to develop some frame of reference around which we can analyze situations. The sighting process nails down the theoretical perspective. 

Assuming that "perfect" horizontal and vertical sight adjustments can be made for the fixed target distance used for sighting, the process is to:

- Take $n$ shots at target
- Measure the horizontal and vertical distance from each shot to the center of the target.
- Calculate the mean horizontal, $\bar{h}$, and vertical, $\bar{v}$, distances from the sample COI to the target center
- Adjust the Horizontal and vertical sights to perfectly correct for the error measured.

It is also assumed that the weapon can maintain this "zero" until the weapon is sighted again. So we're ignoring the real world problems of the weapon getting banged around, undergoing temperature changes, and a zillion other real world details that effect the ballistics of the weapon over time. 

The only thing that we can measure is the COI for the shot pattern, and hence this is generally our point of reference since the wiki is primarily focused on precision not accuracy. 


# Sighting Analysis Assuming Rayleigh Distribution

To analyze the sighting bias we continue to use the symmetric bivariate normal model of shot impacts. 

First using Cartesian coordinates, we assume that both the $n$ measurements along horizontal and vertical axis have a Normal Distribution with the same variance of $\sigma^2$ and that the distributions are both independent and uncorrelated. Thus:   


    $h \sim \mathcal{N}(\mu_h, \sigma^2)$ and $v \sim \mathcal{N}(\mu_v, \sigma^2)$

From the $n$ "sighter" shots and the mean horizontal, $\bar{h}$, and vertical, $\bar{v}$, distances have been calculated which gives the COI for the target. 

## Sighting Accuracy


### Sighting Bias

Breaking the sighting error into parts assuming each error source is independent and normally distributed we have:  


## Sighting Precision

In the Section on [Precision in Shooting](http://ballistipedia.com/index.php?title=What_is_Precision%3F#Precision_in_Shooting) we define the precision of the shooting system by the equation:  


    $\sigma_{System}^2 = \sigma_{Weapon}^2 + \sigma_{Ammunition}^2 + \sigma_{Shooter}^2$

However this equation assumes that the weapon has been sighted once and that the sighting holds for the weapon. Let's consider some different scenarios with multiple sighting adjustments. We will take as a given:

- The value of $\sigma$ for the Rayleigh Distribution for how individual shots are dispersed around the COI from previous experiments. (eg we shot 100 shots one one target to get the "right" number for $\sigma$).

- We'll use multiple targets in the experiment. The left target will be for our sighter groups and the right target to record shots after each sighting adjustment. In other words we'll shoot a group of some sort, adjust the sights, and change the left target, then shoot one shot at the right target.

- We'll shoot $m$ sighter targets with $n$ shots per target.

- After measuring each sighter target we can ask the oracle Carnac the Magnificent for help. Carnac is fickle however. He'll just give us one magic correction factor at a time.

- After each sighter target we'll take one shot at the right target.

- The overall precision of the shots on the right target will thus be:
    $\sigma_{System}^2 = \sigma_{Weapon}^2 + \sigma_{Ammunition}^2 + \sigma_{Shooter}^2 + \sigma_{Sighting}^2$

What we want to understand is how \sigma_{Sighting}^2 changes with different sighting process. 

**Scenario 1** We shot one sighting shot per target.  

For a value of  \sigma_{Sighting}^2 we start with XXXX which is the relationship for the variance of individual shots about one shot. can use the value from the Rayleigh Distribution for the variance based on the sigma value. 


where $\sigma$is the term fit from the Raleigh distribution

**Scenario 2** We shoot </math>n</math> shots per sighter target. 

**Scenario 3** We shoot </math>n</math> shots per sighter target and consult Carnac for the true average $POA^*$ for each group. In other words, the intent was to shoot at the center of the target on each shot. But if we shot an infinite number of shots where would our average aim point have been? 

**Scenario 4** We shoot </math>n</math> shots per sighter target and consult Carnac for the true average $COI^*$ for each group. So Carnac gives us the COI not for the sample of shots which we can measure, but gives the position of the COI if an infinite number of shots had been made.
