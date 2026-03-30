# Sighter Distribution

<p style="text-align:right">... back to [Closed Form Precision](closed-form-precision.md#how-many-sighter-shots-do-you-need)</p>

"Sighting in" a gun involves trying to align its sighting device so that the point of aim coincides with the center of the point of impact (POI).  The statistical challenge with sighting in a gun that has any shot dispersion is that we can only estimate the POI.  In principle we know that increasing the number of samples will improve our estimate, but so long as there is some shot dispersion our POI estimate will be some distance from the true POI center.

To analyze the characteristics of sighting error we continue to use the symmetric bivariate normal model of shot impacts that allowed us to perform [Closed Form Precision](closed-form-precision.md) analysis.  Using Cartesian coordinates we have

$$X, Y \sim N(0, \sigma)$$

I.e., we designate the (unknown) true POI as the origin (0, 0).  Now we take some number *n* of "sighter" shots and calculate their center.  The distance of this sample center from the true POI is

$$R(n) =  \sqrt{\overline{x_i}^2 + \overline{y_i}^2}$$

### What is the distribution of Gun Zero?

*Thanks to Alecos Papadopoulos for the following solution:*

Given the preceding definitions, the sample means are also mean-zero normal random variables with
$$\bar X, \bar Y \sim N(0,\sigma^2/n)$$

Define random variables
$$Z_x = \left(\frac{\bar X}{\sigma/\sqrt n}\right)^2 = \frac n{\sigma^2}\bar X^2 \sim \chi^2(1),\;\;Z_y = \left(\frac{\bar Y}{\sigma/\sqrt n}\right)^2 =\frac n{\sigma^2}\bar Y^2\sim \chi^2(1),$$

Then the random variable 
$$W = Z_x + Z_y =\frac n{\sigma^2}\left(\bar X^2+\bar Y^2\right)\sim \chi^2(2)$$

By the properties of a chi-square random variable, we have
$$W_n=\frac {\sigma^2}nW \sim \text {Gamma}(k=1, \theta = 2\sigma^2/n) = \text{Exp}(2\sigma^2/n)$$
i.e.
$$f_{W_n}(w_n) = \frac {n}{2\sigma^2}\cdot \exp\Big \{-\frac {n}{2\sigma^2} w_n\Big\}$$

Define $R_n = \sqrt {W_n}$. By the change-of-variable formula we have

$$W_n = R_n^2 \Rightarrow \frac {dW_n}{dR_n} = 2R_n$$ and so

$$f_{R_n}(r_n) = 2r_n\frac {n}{2\sigma^2}\cdot \exp\Big \{-\frac {n}{2\sigma^2} r_n^2\Big\} = \frac {r_n}{\alpha^2} \exp\Big \{-\frac {r_n^2}{2\alpha^2} \Big\},\;\;\alpha \equiv \sigma/\sqrt n$$

So $R_n$ follows a Rayleigh distribution with parameter $\alpha = \sigma / \sqrt{n}$.
