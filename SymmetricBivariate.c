/* 
compile line 
gcc -o cor cor.c -lgsl -lgslcblas -lm

Based on simulation from Charles McMillan 10/25/2013
Modified by David Bookstaber 11/17/2013
*/
/* Simulate shots from a symmetric bivariate to see:
* 1. What is distance between sample center and true center?  A: Significant! See output columns 2 & 3. 
* 2. Do sample Rs from sample center underestimate sigma (with the Rayleigh correction)?  A: Yes!
* 3. If not, can we do a Bessel-like correction to sample Rs to get unbiased estimate of sigma?  A: Yes:
*    - Apply Bessel's correction to the R2 estimator.  Multiply sqrt(R2) by c4(2n-1) to get unbiased R estimate.
*    - Apply sqrt(Bessel) to correct [mean radius to sample center] to be [mean radius to actual center]
* 4. Can we use average sample variance in x and y to improve estimate of sigma?  A: Yes, when bivariate is symmetric this gives 2n sample points with same estimate and power as Rayleigh.
* 5. What is best estimator?  A: In decreasing order:
*    - Rayleigh estimator from true center
*    - Rayleigh estimator from sample center = univariate normal estimators (i.e., using sample values from both x and y)
*    - Univariate normal from one dimension
* 6. What is the distribution of N-linked spread statistics?  (Rx is range of sample x, Ry is range of sample y)
*    A. Diagonal = Sqrt(Rx^2 + Ry^2)
*    B. Figure of Merit = (Rx + Ry)/2
*    C. Extreme Spread = max(sqrt((x_i - x_j)^2 - (y_i - y_j)^2))
*/

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_gamma.h>

#define ITERATIONS 1000000
#define PI 3.14159265358979
#define MAX_N 100 // Max group size to simulate
#define SIGMA 1.0

#define SQR(x) ((x)*(x))

double c4 (long n);
int shoot (long n, gsl_rng *r, double x0, double y0, double sigma1, double sigma2, double rho, double *x, double *y);

int main (void){
	/* Initialize GSL random number generator */
	const gsl_rng_type *T;
	gsl_rng *rng;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rng = gsl_rng_alloc(T);

	long i, j, n = ITERATIONS;
	unsigned group; // Group size to simulate -- iterates [2 to MAX_N]
	double sigma = SIGMA;
	double cBessel, cRayleigh, c4_n, c4_2n1;
	// Samples for each shot
	double *x, *y;
	// Variables with one value per group
	double *centerDistance; // Distance from sample center to true center (the origin {0,0})
	double *centerR_bar; // Average of radii from sample center
	double *originR_bar; // Average of radii from true center (origin)
	double *centerR2_estimate; // Rayleigh s^2 estimator using radii to sample center
	double *originR2_estimate; // Rayleigh s^2 estimator using radii to true center
	double *centerR_estimate; // Rayleigh s estimator using radii to sample center (biased)
	double *originR_estimate; // Rayleigh s estimator using radii to true center (biased)
	double *varX_bar, *varY_bar, *var_bar; // Sample variances
	double *sigmaX_bar, *sigmaY_bar, *sigma_bar; // Sample sigmas of X and Y (biased)
	double *diagonal, *FoM, *extremeSpread;
	// Variables for averaging over all samples in a group size
	double x_min, x_max, x_bar, x_range, y_min, y_max, y_bar, y_range, radius;
	unsigned left, right, top, bottom; // Index of extreme shots in each group

	centerDistance = (double *) malloc(ITERATIONS*sizeof(double));
	centerR_bar = (double *) malloc(ITERATIONS*sizeof(double));
	originR_bar = (double *) malloc(ITERATIONS*sizeof(double));
	centerR2_estimate = (double *) malloc(ITERATIONS*sizeof(double));
	originR2_estimate = (double *) malloc(ITERATIONS*sizeof(double));
	centerR_estimate = (double *) malloc(ITERATIONS*sizeof(double));
	originR_estimate = (double *) malloc(ITERATIONS*sizeof(double));
	varX_bar = (double *) malloc(ITERATIONS*sizeof(double));
	varY_bar = (double *) malloc(ITERATIONS*sizeof(double));
	var_bar = (double *) malloc(ITERATIONS*sizeof(double));
	sigmaX_bar = (double *) malloc(ITERATIONS*sizeof(double));
	sigmaY_bar = (double *) malloc(ITERATIONS*sizeof(double));
	sigma_bar = (double *) malloc(ITERATIONS*sizeof(double));
	diagonal = (double *) malloc(ITERATIONS*sizeof(double));
	FoM = (double *) malloc(ITERATIONS*sizeof(double));
	extremeSpread = (double *) malloc(ITERATIONS*sizeof(double));

	printf("Iterations=,%ld,Sigma=,%lf,Variance=,%lf,MeanRadius=,%lf\n\n", n, sigma, SQR(sigma), sigma*sqrt(PI/2.0));

	/* ************************** */
	/* Outer loop over group size */

	// Header for columns
	printf("Group Size,");
	printf("Avg |Center|,");
	printf("StDev |Center|,");
	printf("Mean Radius to Center *Sqrt(cB),");
	printf("Mean Radius to Origin,");
	printf("CenterR2 Estimate *cB,");
	printf("CenterR2 StDev,");
	printf("CenterR Estimate *cN(2n-1),");
	printf("CenterR StDev,");
	printf("OriginR2 Estimate,");
	printf("OriginR2 StDev,");
	printf("OriginR Estimate *cR,");
	printf("OriginR StDev,");
	printf("VarianceX Estimate,");
	printf("VarianceX StDev,");
	printf("Variance Estimate,");
	printf("Variance StDev,");
	printf("SigmaX Estimate *cN,");
	printf("SigmaX StDev,");
	printf("Sigma Estimate *cN(2n-1),");
	printf("Sigma StDev,");
	printf("Bessel Correction,");
	printf("Normal Correction,");
	printf("Rayleigh Correction,");
	printf("Diagonal,");
	printf("Diagonal StDev,");
	printf("FoM,");
	printf("FoM StDev,");
	printf("Extreme Spread,");
	printf("Extreme Spread Stdev");
	printf("\n");
	for (group = 2; group <= MAX_N; group += 1){
		x = (double *) malloc(group*sizeof(double));
		y = (double *) malloc(group*sizeof(double));

		// Correction factors given sample (group) size
		c4_n = c4(group);
		c4_2n1 = c4(2*group - 1);
		// cRayleigh = sqrt(group / PI) * pow(4, group) * gsl_sf_fact(group) * gsl_sf_fact(group-1) / gsl_sf_fact(2*group); // Overflows factorial for group > 85
		cRayleigh = exp(log(sqrt(group / PI)) + group * log(4) + gsl_sf_lnfact(group) + gsl_sf_lnfact(group-1) - gsl_sf_lnfact(2*group)); // This version avoids overflow
		cBessel = group / (group - 1.0);

		// Run N iterations of group shots
		for (j = 0; j < n; j++){
			// For each iteration:
			shoot(group, rng, 0, 0, sigma, sigma, 0, x, y);

			// Compute group center and range
			left = right = top = bottom = 0;
			x_bar = y_bar = 0;
			x_min = x_max = x[0];
			y_min = y_max = y[0];
			for (i = 0; i < group; i++){
				if(x[i] < x_min){
					x_min = x[i];
					left = i;
				}
				if(x[i] > x_max){
					x_max = x[i];
					right = i;
				}
				x_bar += x[i];

				if(y[i] < y_min){
					y_min = y[i];
					bottom = i;
				}
				if(y[i] > y_max){
					y_max = y[i];
					top = i;
				}
				y_bar += y[i];
			}
			x_bar /= (double)group;
			y_bar /= (double)group;
			x_range = x_max - x_min;
			y_range = y_max - y_min;
			centerDistance[j] = sqrt(SQR(x_bar) + SQR(y_bar));
			diagonal[j] = sqrt(SQR(x_range) + SQR(y_range));
			FoM[j] = (x_range + y_range) / 2.0;
			// Compute extreme spread by checking distance between the four extreme shots
			extremeSpread[j] = SQR(x[top] - x[bottom]) + SQR(y[top] - y[bottom]);
			radius = SQR(x[left] - x[right]) + SQR(y[left] - y[right]);
			if(radius > extremeSpread[j]) extremeSpread[j] = radius;
			radius = SQR(x[top] - x[right]) + SQR(y[top] - y[right]);
			if(radius > extremeSpread[j]) extremeSpread[j] = radius;
			radius = SQR(x[left] - x[top]) + SQR(y[left] - y[top]);
			if(radius > extremeSpread[j]) extremeSpread[j] = radius;
			radius = SQR(x[left] - x[bottom]) + SQR(y[left] - y[bottom]);
			if(radius > extremeSpread[j]) extremeSpread[j] = radius;
			radius = SQR(x[bottom] - x[right]) + SQR(y[bottom] - y[right]);
			if(radius > extremeSpread[j]) extremeSpread[j] = radius;
			extremeSpread[j] = sqrt(extremeSpread[j]);

			// Compute variances, radii and estimators
			varX_bar[j] = varY_bar[j] = 0;
			centerR_bar[j] = centerR2_estimate[j] = 0;
			originR_bar[j] = originR2_estimate[j] = 0;
			for (i = 0; i < group; i++){
				varX_bar[j] += SQR(x[i] - x_bar);
				varY_bar[j] += SQR(y[i] - y_bar);

				radius = sqrt(SQR(x[i] - x_bar) + SQR(y[i] - y_bar));
				centerR_bar[j] += radius;
				centerR2_estimate[j] += SQR(radius);

				radius = sqrt(SQR(x[i]) + SQR(y[i]));
				originR_bar[j] += radius;
				originR2_estimate[j] += SQR(radius);
			}
			varX_bar[j] /= (group - 1.0);
			varY_bar[j] /= (group - 1.0);
			var_bar[j] = (varX_bar[j] + varY_bar[j]) / 2.0;
			sigmaX_bar[j] = sqrt(varX_bar[j]) * c4_n; // Single dimension's sigma requires c4 correction
			sigmaY_bar[j] = sqrt(varY_bar[j]) * c4_n;
			sigma_bar[j] = sqrt(var_bar[j]) * c4_2n1; // Sigma estimate from both dimensions' samples requires c4(2n-1) correction
			centerR_bar[j] /= group / sqrt(cBessel);
			originR_bar[j] /= group;
			centerR2_estimate[j] /= (2.0*group) / cBessel; // Rayleigh R2 Estimate from sample center requires Bessel correction
			originR2_estimate[j] /= (2.0*group);
			centerR_estimate[j] = sqrt(centerR2_estimate[j]) * c4_2n1; // Rayleigh R estimate from sample center requires Bessel and c4(2n-1) correction
			originR_estimate[j] = sqrt(originR2_estimate[j]) * cRayleigh; // Rayleigh R estimate from true center requires Rayleigh correction
		}

		printf("%ld,", group);
		printf("%lf,", gsl_stats_mean(centerDistance, 1, n));
		printf("%lf,", gsl_stats_sd(centerDistance, 1, n));
		printf("%lf,", gsl_stats_mean(centerR_bar, 1, n));
		printf("%lf,", gsl_stats_mean(originR_bar, 1, n));
		printf("%lf,", gsl_stats_mean(centerR2_estimate, 1, n));
		printf("%lf,", gsl_stats_sd(centerR2_estimate, 1, n));
		printf("%lf,", gsl_stats_mean(centerR_estimate, 1, n));
		printf("%lf,", gsl_stats_sd(centerR_estimate, 1, n));
		printf("%lf,", gsl_stats_mean(originR2_estimate, 1, n));
		printf("%lf,", gsl_stats_sd(originR2_estimate, 1, n));
		printf("%lf,", gsl_stats_mean(originR_estimate, 1, n));
		printf("%lf,", gsl_stats_sd(originR_estimate, 1, n));
		printf("%lf,", gsl_stats_mean(varX_bar, 1, n));
		printf("%lf,", gsl_stats_sd(varX_bar, 1, n));
		printf("%lf,", gsl_stats_mean(var_bar, 1, n));
		printf("%lf,", gsl_stats_sd(var_bar, 1, n));
		printf("%lf,", gsl_stats_mean(sigmaX_bar, 1, n));
		printf("%lf,", gsl_stats_sd(sigmaX_bar, 1, n));
		printf("%lf,", gsl_stats_mean(sigma_bar, 1, n));
		printf("%lf,", gsl_stats_sd(sigma_bar, 1, n));
		printf("%lf,", cBessel);
		printf("%lf,", c4_n);
		printf("%lf,", cRayleigh);
		printf("%lf,", gsl_stats_mean(diagonal, 1, n));
		printf("%lf,", gsl_stats_sd(diagonal, 1, n));
		printf("%lf,", gsl_stats_mean(FoM, 1, n));
		printf("%lf,", gsl_stats_sd(FoM, 1, n));
		printf("%lf,", gsl_stats_mean(extremeSpread, 1, n));
		printf("%lf", gsl_stats_sd(extremeSpread, 1, n));
		printf("\n");

		free(x);
		free(y);
	}

	gsl_rng_free(rng);
	free(centerDistance);
	free(centerR_bar);
	free(originR_bar);
	free(centerR2_estimate);
	free(originR2_estimate);
	free(centerR_estimate);
	free(originR_estimate);
	free(varX_bar);
	free(varY_bar);
	free(var_bar);
	free(sigmaX_bar);
	free(sigmaY_bar);
	free(sigma_bar);
	free(diagonal);
	free(FoM);
	free(extremeSpread);

	return 0;
}

// Gaussian correction factor
double c4(long n){
	return 1./(1. - 1./(4.*n) - 7./(32.*n*n) - 19./(128.*n*n*n));
}

/*****************************
* Generate a set of shots 
* n - number of shots to fire
* r - pointer to random number generator
* (x0, y0) - aim point doesn't have to be zero
* (sigma1, sigma2) - Standard deviations in the x and y directions
* rho - covariance
* (*x, *y) - the array containing the shots - have to allocated externally
***************************** */
int shoot (long n, gsl_rng *r, double x0, double y0, double sigma1, double sigma2, double rho, double *x, double *y){
	long i;
	for (i = 0; i < n; i++){
		gsl_ran_bivariate_gaussian (r, sigma1, sigma2, rho, &x[i], &y[i]);
//		x[i] += x0;
//		y[i] += y0; // Caution: this value seems to get changed after first call -- probable memory leak somewhere!
	}
	return 0;
}

/* ******************************* */      
