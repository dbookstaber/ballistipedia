# Extreme Spread

!!! warning
    This page is a draft and needs review!


# Experimental Summary

|  |  |
| --- | --- |
| Given |  |
| Assumptions |  |
| Data transformation | Identify two holes, $i, j$ which are the farthest apart and measure $ES$. |
| Experimental Measure | $ES$ |


## Given

The requirements for this test are very basic. Just a target with $n$ shots, and some measuring device. Assuming an Extreme spread of under 6 inches then a vernier caliper is used. A measurement is possible to a few thousandths of an inch which is vastly more precision than is usually required. From longer distance a ruler, or perhaps a tape measure.

## Assumptions

None are needed to make measurement. However some points are worth considering.

- The same ES measurement could result from a vertical group to a round group. If the shooting process can vary that much then the ES measurement won't give any indication of the change.

> > If the shot patterns aren't "fairly" round, then using the measurement makes little sense. For instance if muzzle velocity variations are severe, then the vertical range will dominate the ES measurement. Muzzle velocity variations would correlate better with vertical range than with ES.

- Making assumptions about the dispersion will enable theoretical predictions about the ES measurement. It must be realized that the theoretical solution, assuming the Rayleigh distribution and using Monte Carlo simulation, isn't some arbitrary goal, it is the best case scenario.

## Data transformation

The data transformation for a human has simple requirements, just the ability to locate the holes which are the furthest apart and measure the distance between them. If the target has a ragged hole it can be a bit tricky, but the edges of the hole should have enough curvature to make shot location possible.

If measuring on the range, then the center of the hole is difficult to locate. Typically a vernier caliper (cheap is fine!) would be used to measure the distance from the outside edges of the holes, then the bullet caliber subtracted to get a c-t-c measurement. 

> !!! warning
    A cheap ($10-$20) vernier caliper works fine. There is no need for a $2,000 one that measures to 1/10,000^th^ of an inch and has National Bureau of Standards calibration. The vernier caliper is nice for the c-t-c measurement because the knife edges will be parallel and won't obscure the edges of the bullet hole. Thus it is easy to accurately place both of the knife edges on a tangent to the curved bullet holes.


If using a computer then the center location would be a matter programming. For example a mouse might be used simply to point out the holes, or to drop a dot at the center of the hole, or to drag a circle over the hole. The computer would then make the c-t-c measurement.

## Experimental Measure

No calculation needs to be done to get the measurement. The single physical measurement is the data sought for the target.

# See Also

[Projectile Dispersion Classifications](projectile-dispersion-classifications.md) - A discussion of the different cases for projectile dispersion
