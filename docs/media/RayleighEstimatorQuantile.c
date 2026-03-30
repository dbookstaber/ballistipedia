/* 
 * compile line 
 * gcc -o cor cor.c -lgsl -lgslcblas -lm
 * 
 * David Bookstaber 3/18/2014
 *
 * Simulate shots from a symmetric bivariate to
 * calculate Rayleigh estimate of sigma
 * and output quantile function (in steps of 0.5%) of estimated sigma
 * for each sample size 2-50.
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
#define MAX_N 50 // Max group size to simulate
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
	double *centerR2_estimate; // Rayleigh s^2 estimator using radii to sample center
	double *originR2_estimate; // Rayleigh s^2 estimator using radii to true center
	double *centerR_estimate; // Rayleigh s estimator using radii to sample center (biased)
	double *originR_estimate; // Rayleigh s estimator using radii to true center (biased)
	// Variables for averaging over all samples in a group size
	double x_bar, y_bar, radius;

	centerR2_estimate = (double *) malloc(ITERATIONS*sizeof(double));
	originR2_estimate = (double *) malloc(ITERATIONS*sizeof(double));
	centerR_estimate = (double *) malloc(ITERATIONS*sizeof(double));
	originR_estimate = (double *) malloc(ITERATIONS*sizeof(double));

	printf("Iterations=,%ld,Sigma=,%lf,Variance=,%lf,MeanRadius=,%lf\n\n", n, sigma, SQR(sigma), sigma*sqrt(PI/2.0));

	/* ************************** */
	/* Outer loop over group size */

	// Header
	printf("Shots,Mean,Stdev,Skew,Kurtosis");
	for(i = 1;i<=200;i++)
		printf(",%lf",i/200.0);
	printf("\n");

	for (group = 2; group <= MAX_N; group += 1){
		x = (double *) malloc(group*sizeof(double));
		y = (double *) malloc(group*sizeof(double));

		// Correction factors given sample (group) size
		c4_2n1 = c4(2*group - 1);
		// cRayleigh = sqrt(group / PI) * pow(4, group) * gsl_sf_fact(group) * gsl_sf_fact(group-1) / gsl_sf_fact(2*group); // Overflows factorial for group > 85
		cRayleigh = exp(log(sqrt(group / PI)) + group * log(4) + gsl_sf_lnfact(group) + gsl_sf_lnfact(group-1) - gsl_sf_lnfact(2*group)); // This version avoids overflow
		cBessel = group / (group - 1.0);

		// Run N iterations of group shots
		for (j = 0; j < n; j++){
			// For each iteration:
			shoot(group, rng, 0, 0, sigma, sigma, 0, x, y);

			x_bar = y_bar = 0;
			for (i = 0; i < group; i++){
				x_bar += x[i];
				y_bar += y[i];
			}
			x_bar /= (double)group;
			y_bar /= (double)group;

			// Compute variances, radii and estimators
			centerR2_estimate[j] = 0;
			originR2_estimate[j] = 0;
			for (i = 0; i < group; i++){
				radius = sqrt(SQR(x[i] - x_bar) + SQR(y[i] - y_bar));
				centerR2_estimate[j] += SQR(radius);

				radius = sqrt(SQR(x[i]) + SQR(y[i]));
				originR2_estimate[j] += SQR(radius);
			}
			centerR2_estimate[j] /= (2.0*group) / cBessel; // Rayleigh R2 Estimate from sample center requires Bessel correction
			originR2_estimate[j] /= (2.0*group);
			centerR_estimate[j] = sqrt(centerR2_estimate[j]) * c4_2n1; // Rayleigh R estimate from sample center requires Bessel and c4(2n-1) correction
			originR_estimate[j] = sqrt(originR2_estimate[j]) * cRayleigh; // Rayleigh R estimate from true center requires Rayleigh correction
		}

		gsl_sort(centerR_estimate, 1, n);

		printf("%ld,", group);
		printf("%lf,", gsl_stats_mean(centerR_estimate, 1, n));
		printf("%lf,", gsl_stats_sd(centerR_estimate, 1, n));
		printf("%lf,", gsl_stats_skew(centerR_estimate, 1, n));
		printf("%lf", gsl_stats_kurtosis(centerR_estimate, 1, n));
		for(i = 1;i <= 200;i++)
			printf(",%lf", centerR_estimate[(int)(i*n/200.0)-1]);
		printf("\n");

		free(x);
		free(y);
	}

	gsl_rng_free(rng);
	free(centerR2_estimate);
	free(originR2_estimate);
	free(centerR_estimate);
	free(originR_estimate);

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
