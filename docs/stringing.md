# Stringing

For a "large" number of shots, to the extent that either $\sigma_h \neq \sigma_v$ or $\rho \neq 0$ then elliptical shot groups will result instead of circular shot groups. A "small" number of shots may have any pattern due to chance. So at least 10 shots should be used before any visual impression is taken. 

If the shot groups are not round then there are three options. 

> (1) Isolate the error source experimentally and remove it (for instance weigh gunpowder carefully to remove vertical stringing).
> > Obviously the experimental reason for stinging may not be obvious and easy to remove. Experimental designs to isolate and quantify the source of the stringing are beyond this basic discussion at this point, but possible.
> (2) Use a mathematical model for target analysis that allows for stringing (i.e. fit an ellipse, not a circle).
> (3) Use exterior ballistics to measure the factor causing the elliptical grouping independently, and correct the target data for that factor.


# Horizontal Stringing

Given that $\sigma_h > \sigma_v$ then the shot positions on the target would be elliptical in shape with the longer axis of the ellipse along the horizontal axis. 

We do know that targets can often exhibit vertical or horizontal stringing as evidenced by an elliptical shaped group along the vertical or horizontal axis respectively. Obviously in such cases $\sigma_h \neq \sigma_v$. 

The primary source of horizontal stringing is crosswind.  

> Given that we are measuring the horizontal and vertical variances, if we measure the wind while shooting we can bound and remove a “wind correction” term from that axis.  E.g., "Suppose the orthogonal component of wind is ranging at random from 0-10mph during the shooting.  Given lag-time *t* this will expand the no-wind horizontal dispersion at the target by $f(\sigma_{wind}^2)$. Wind deflection is a function of the ballistic curve and distance, but can be expressed as a simple product of the cross-wind velocity and lag time known as the "lag rule" [^1][^2]

> Since variances are additive we could adjust $\sigma_h$ via the equation ${\sigma'}_h^2 = \sigma_h^2 - f(\sigma_{wind}^2)$.

# Vertical Stringing

Given that $\sigma_h < \sigma_v$ then the shot positions on the target would be elliptical in shape with the longer axis of the ellipse along the vertical axis. 

A typical source of vertical stringing is muzzle velocity. 

> We can actually measure the muzzle velocity for each shot with a chronograph and then correct for the muzzle velocity dispersion. E.g., If standard deviation of muzzle velocity is $s_{mv}$ then, given the bullet's ballistic model for the given target distance, the vertical spread attributable to that is some $s_v$.  Here too we can remove this known source of dispersion from our samples via the equation ${\sigma'}_v^2 = \sigma_v^2 - f(\sigma_{mv}^2)$.  and/or  appropriately rescale the horizontal component of the shot positions, depending on the measurement being used. This adjustment is shown in several of the examples:
> > * [22LR CCI 40gr HV 40-shot 100-yard Example](22lr-cci-40gr-hv-40-shot-100-yard-example.md)
> > * [300BLK Subsonic 20-shot 100-yard Example](300blk-subsonic-20-shot-100-yard-example.md)

# Stringing due to Horizontal-Vertical Correlation

Given that $\rho \neq 0$ then the shot positions on the target would be elliptical in shape with the longer axis at some angle to both the horizontal and vertical axes.

Even though exterior ballistics does indicate an interaction between the horizontal and vertical dispersion of gunshots, at a "short" distance the interaction is assumed negligible, or at least consistent enough so that the differences between shots are negligible. Therefore, most statistical measures documented in this wiki implicitly assume $\rho = 0$. In general if $\rho \neq 0$, then there would be an elliptical shaped group with the major axis oriented at some angle to the horizontal or vertical axis.

# References


[^1]: Bryan Litz, ''Applied Ballistics for Long Range Shooting, 2<sup>nd</sup> Edition'' (2011) A4

[^2]: Robert McCoy, ''Modern Exterior Ballistics, 2<sup>nd</sup> Edition'' (2012) 7.27.
