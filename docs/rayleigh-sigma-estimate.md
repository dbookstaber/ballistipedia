# Estimating parameter $\sigma$ for the Rayleigh distribution

Let $X, Y$ be uncorrelated, jointly normally distributed cartesian coordinates with means $\mu_{X}, \mu_{Y}$, and equal variances $\sigma_{X}^{2} = \sigma_{Y}^{2}$. The radius of an $(X,Y)$-coordinate pair around the true mean is:

$R := \sqrt{(X-\mu_{X})^{2} + (Y-\mu_{Y})^{2}}$

With $N$ observations of $(x,y)$-coordinates, $\text{SSR}$ is the sum of squared radii:

$\begin{array}{rcl} \text{SSR} &:=& \sum\limits_{i=1}^{N} r_{i}^{2} = \sum\limits_{i=1}^{N} ((x_{i}-\mu_{X})^{2} + (y_{i}-\mu_{Y})^{2})\\  &=& \sum\limits_{i=1}^{N} (x_{i}-\mu_{X})^{2} +  \sum\limits_{i=1}^{N} (y_{i}-\mu_{Y})^{2} \end{array}$

The Maximum-Likelihood-estimate of the variance of $r$ is the total variance of all $2N$ separate $x$- and $y$-coordinates:

$\widehat{\sigma^{2}_{ML}} = \frac{1}{2N} \cdot \text{SSR}$

## Unknown true center
### Deriving Singh's (1992) $C_{2}$ estimator

When the $(\mu_{X}, \mu_{Y})$-center is estimated by $(\bar{x}, \bar{y})$, the [Bessel correction](closed-form-precision.md#bessel-correction-factor) is required to make the ML-estimate of the variance unbiased:

$\widehat{\sigma^{2}_{ub}} = \frac{N}{N-1} \cdot \frac{1}{2N} \cdot \text{SSR} = \frac{1}{2(N-1)} \cdot \text{SSR}$

However, taking the square root of the unbiased variance estimate makes it biased ([Jensen's inequality](http://en.wikipedia.org/wiki/Jensen%27s_inequality)). Specifically, since the square root is concave, the bias is negative and $\sqrt{\widehat{\sigma^{2}_{ub}}}$ underestimates $\sigma$. Another correction factor is needed to counter this effect. To this end, let $Q$ be a $\chi$-distributed variable with $k$ degrees of freedom. Then its mean is:

$E(Q) = \sqrt{2} \cdot \frac{\Gamma\left(\frac{k+1}{2}\right)}{\Gamma\left(\frac{k}{2}\right)}$

When $X$ is a normally distributed random variable with $n$ observations, and $s^{2} := \frac{1}{n-1} \sum\limits_{i=1}^{n}(x-\bar{x})^{2}$ is the Bessel-corrected variance estimate, $E(s^{2}) = \sigma^{2}$. Then [$c_{4}(n)$ is the correction factor](closed-form-precision.md#gaussian-correction-factor) such that $E(s) = c_{4}(n) \cdot \sigma$. With $Q$ as given above, $c_{4}(n) = \frac{1}{\sqrt{k}} \cdot E(Q)$ with $k := n-1$:

$c_{4}(n) = \sqrt{\frac{2}{n-1}} \cdot \frac{\Gamma\left(\frac{n}{2}\right)}{\Gamma\left(\frac{n-1}{2}\right)}$

We lose 2 degrees of freedom estimating the center.  Therefore we really only have 2N - 2 degrees of freedom, so we set $n := 2N-1$ so that $c_{4}(n)$ is the scaled mean of a $\chi$-distributed variable with $k = n-1 = 2N-1-1 = 2N-2$ degrees of freedom. This gives:

$\begin{array}{rcl} \widehat{\sigma_{ub}} &=& \frac{1}{c_{4}(2N-1)} \cdot \sqrt{\widehat{\sigma^{2}_{ub}}} \\  &=& \frac{1}{\sqrt{\frac{2}{2N-1-1}} \cdot \frac{\Gamma\left(\frac{2N-1}{2}\right)}{\Gamma\left(\frac{2N-1-1}{2}\right)}} \cdot \sqrt{\frac{N}{N-1} \cdot \frac{1}{2N} \cdot \text{SSR}} \\  &=& \sqrt{\frac{2(N-1)}{2}} \cdot \frac{\Gamma\left(\frac{2(N-1)}{2}\right)}{\Gamma\left(\frac{2N-1}{2}\right)} \cdot \sqrt{\frac{N}{N-1} \cdot \frac{1}{2N} \cdot \text{SSR}} \\   &=& \sqrt{N} \cdot \frac{\Gamma(N-1)}{\Gamma\left(\frac{2N-1}{2}\right)} \cdot \sqrt{\frac{1}{2N} \cdot \text{SSR}} \\  &=& \frac{\Gamma(N-1)}{\Gamma\left(N - \frac{1}{2}\right)} \cdot \sqrt{\frac{1}{2} \cdot \text{SSR}} \end{array}$

The penultimate expression is [Singh's 1992](cep-literature.md#singh1992) definition 3.3 for estimator $C_{2}$, also given by [Moranda (1959)](cep-literature.md#moranda1959).

## Known true center
### Deriving Singh's (1992) $C_{1}$ estimator

When the $(\mu_{X}, \mu_{Y})$-center is known, the Bessel correction is **not** required to make the Maximum-Likelihood-estimate of the variance unbiased.  In this case we simply set $n := 2N+1$ for the $c_{4}(n)$ correction factor so that $c_{4}(n)$ is the scaled mean of a $\chi$-distributed variable with $k = n-1 = 2N+1-1 = 2N$ degrees of freedom. This gives:

$\begin{array}{rcl} \widehat{\sigma_{ub}} &=& \frac{1}{c_{4}(2N+1)} \cdot \sqrt{\widehat{\sigma^{2}_{ub}}} \\  &=& \frac{1}{\sqrt{\frac{2}{2N+1-1}} \cdot \frac{\Gamma\left(\frac{2N+1}{2}\right)}{\Gamma\left(\frac{2N+1-1}{2}\right)}} \cdot \sqrt{\frac{1}{2N} \cdot \text{SSR}} \\  &=& \sqrt{N} \cdot \frac{\Gamma(N)}{\Gamma\left(\frac{2N+1}{2}\right)} \cdot \sqrt{\frac{1}{2N} \cdot \text{SSR}} \\   &=& \frac{\Gamma(N)}{\Gamma\left(N + \frac{1}{2}\right)} \cdot \sqrt{\frac{1}{2} \cdot \text{SSR}} \end{array}$

The penultimate expression is [Singh's 1992](cep-literature.md#singh1992) definition 2.2 for estimator $C_{1}$.
