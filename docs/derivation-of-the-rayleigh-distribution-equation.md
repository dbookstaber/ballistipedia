# Bivariate Normal Distribution
Starting only with the assumptions that the horizontal and vertical measurements are normally distributed as notated by:  


$$
h \sim \mathcal{N}(\mu_h,\sigma_h^2), \, \, v \sim \mathcal{N}(\mu_v,\sigma_v^2)
$$


then the horizontal and vertical measures follow the general bivariate normal distribution which is given by the following equation:

$f(h,v) =       \frac{1}{2 \pi  \sigma_h \sigma_v \sqrt{1-\rho^2}}       \exp\left(         -\frac{1}{2(1-\rho^2)}\left[           \frac{(h-\mu_h)^2}{\sigma_h^2} +           \frac{(v-\mu_v)^2}{\sigma_v^2} -           \frac{2\rho(h-\mu_h)(v-\mu_v)}{\sigma_h \sigma_v}         \right]       \right)$

# The Hoyt Distribution
Simplification of the Bivariate Normal Distribution to the Hoyt Distribution

## Derivation From the Bivariate Normal distribution
A simple translation of the Cartesian Coordinate System converts the Bivariate Normal distribution to the Hoyt distribution. This translation will not affect measurements about COI, but it would of course affect measurements which are measured about POA. 

Given a translation to point $(\mu_h, \mu_v)$ then let:


$$
h_* =  h - \mu_h, \, \, v_* =  v - \mu_v
$$


Since the derivative of $h_*$ with respect to $(h - \mu_h)$ is 1, (and similarity for $v_*$) then no change results to the integration constant of the function. Thus $h_*$ can be substituted for $(h - \mu_h)$ and $v_*$ for $(v - \mu_v)$. At this point the asterisk subscript is superfluous and will be dropped, giving the Hoyt distribution: 

$f(h,v) =       \frac{1}{2 \pi  \sigma_h \sigma_v \sqrt{1-\rho^2}}       \exp\left(         -\frac{1}{2(1-\rho^2)}\left[           \frac{h^2}{\sigma_h^2} +           \frac{v^2}{\sigma_v^2} -           \frac{2\rho h v }{\sigma_h \sigma_v}         \right]       \right)$


# The Rayleigh distribution

The Rayleigh Distribution makes the following simplifying assumptions to the general bivariate normal distribution:

- Horizontal and vertical dispersion are independent.
- $\sigma_h = \sigma_v$ (realistically $\sigma_h \approx \sigma_v$)
- $\rho = 0$
- No Fliers
for which the PDF for any shot, $i$, around the horizontal and vertical point $(\mu_h, \mu_v)$ is given by:

$PDF(r; \sigma_{\Re}) = \frac{r}{\sigma_{\Re}^2 }       \exp\left(         - \frac{r^2}{2\sigma_{\Re}^2}        \right)$

where


$$
\sigma_{\Re} = \sigma_h = \sigma_v, \, \, r = \sqrt{h_i - \mu_h)^2 + sqrt(v_i - \mu_v)^2}
$$


## Proof

Using the assumptions in the first section, the distribution of an individual shot is easily simplified from the Bivariate Normal Distribution which has the equation:

$f(h,v) =       \frac{1}{2 \pi  \sigma_h \sigma_v \sqrt{1-\rho^2}}       \exp\left(         -\frac{1}{2(1-\rho^2)}\left[           \frac{(h-\mu_h)^2}{\sigma_h^2} +           \frac{(v-\mu_v)^2}{\sigma_v^2} -           \frac{2\rho(h-\mu_h)(v-\mu_v)}{\sigma_h \sigma_v}         \right]       \right)$

By substituting $\rho = 0$ the equation reduces to:

$f(h,v) =       \frac{1}{2 \pi  \sigma_h \sigma_v }       \exp\left(         -\frac{1}{2}\left[           \frac{(h-\mu_h)^2}{\sigma_h^2} +           \frac{(v-\mu_v)^2}{\sigma_v^2}          \right]       \right)$

Since $\sigma_h$ and $\sigma_v$ are equal, substitute $\sigma$ for each, then collect terms in the exponential, after which the equation reduces to:

$f(h,v) =       \frac{1}{2 \pi  \sigma^2 }       \exp\left(         -\left[           \frac{(h-\mu_h)^2 + (v-\mu_h)^2}{2\sigma^2}          \right]       \right)$

Letting $r^2 = (h-\mu_h)^2 + (v-\mu_v)^2$ the equation becomes: 

$f(h,v) =       \frac{1}{2 \pi  \sigma^2 }       \exp\left(         - \frac{r^2}{2\sigma^2}        \right)$

Now transforming to the polar coordinate system:

    - ok, here I'm lost **

and finally:

$f(r) =       \frac{r}{\sigma^2 }       \exp\left(         - \frac{r^2}{2\sigma^2}        \right)$

### Accuracy of $n$ Sighting Shots

Given the assumptions in the starting section we again substitute $\sigma$ for both $\sigma_h$ and $\sigma_v$. This simplifies the distributions of $h$ and $v$ to:  


$$
h \sim \mathcal{N}(\mu_h,\sigma^2), \, \, v \sim \mathcal{N}(\mu_v,\sigma^2)
$$


Now we take some number $n$ of shots $( n \geq 1)$and calculate their centers $\bar{h}$ and $\bar{v}$ which will be normal distributions as well. 


$$
\bar{h} \sim \mathcal{N}(\mu_h,\sigma^2/n),  \, \, bar{v} \sim \mathcal{N}(\mu_v,\sigma^2/n)
$$


Let $r_n$ be the distance of this sample center $(\bar{h}, \bar{v})$ from the true distribution center $(\mu_h, \mu_v)$ as:  


$$
r_n = \sqrt{(\bar{h}-\mu_h)^2 + (\bar{v}-\mu_v)^2}
$$


Define random variables $Z_h$ and $Z_v$ as the squared *Studentized* horizontal and vertical errors by dividing by the respective standard deviations. Each of these variables with have a Chi-Squared Distribution with one degree of freedom.

$Z_h = \left(\frac{(\bar{h}-\mu_h)}{\sigma/\sqrt n}\right)^2 = \frac n{\sigma^2}(\bar{h}-\mu_h)^2 \sim \chi^2(1)$


$Z_v = \left(\frac{(\bar{v}-\mu_v)}{\sigma/\sqrt n}\right)^2 = \frac n{\sigma^2}(\bar{v}-\mu_v)^2\sim \chi^2(1)$


Define random the variable $W$ which will have a Chi-Squared Distribution with two degrees of freedom as:  
 


$$
W = Z_x + Z_y =\frac n{\sigma^2}\left((\bar{v}-\mu_v)^2+(\bar{v}-\mu_v)^2\right)\sim \chi^2(2)
$$


Rescale the variable $W$ by $\frac {\sigma^2}{n}$ and denote the new variable $w_n$:  


$w_n=\frac {\sigma^2}nW$ and note that $w_n=r_n^2$

By the properties of a chi-square random variable, we have:  


$$
w_n \sim \text {Gamma}(k=1, \theta = 2\sigma^2/n) = \text{Exp}(2\sigma^2/n)
$$


so:  


$$
PDF(w_n) = \frac {n}{2\sigma^2}\cdot \exp\Big \{-\frac {n}{2\sigma^2} w_n\Big\}
$$


But from above $r_n = \sqrt {w_n}$. By the change-of-variable formula we have

$w_n = r_n^2 \Rightarrow \frac {dw_n}{dr_n} = 2r_n$


and so:  


$$
PDF(r_n) = 2r_n\frac {n}{2\sigma^2}\cdot \exp\Big \{-\frac {n}{2\sigma^2} r_n^2\Big\} = \frac {r_n}{\alpha^2} \exp\Big \{-\frac {r_n^2}{2\alpha^2} \Big\},\;\;\alpha \equiv \sigma/\sqrt n
$$


So for any number of shots $n$, the expected accuracy is given by $r_n$ follows a Rayleigh distribution with parameter $\alpha = \sigma / \sqrt{n}$ where $\sigma$ is the Rayleigh shape factor for one shot. 

**Thanks to Alecos Papadopoulos for the solution.**
