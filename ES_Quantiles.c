/* 
 * compile line 
 * gcc -o cor cor.c -lgsl -lgslcblas -lm
 * 
 * David Bookstaber 3/18/2014
 *
 * Given n shots, simulate the average extreme spread of s groups of size g, for all s*g = n
 * and where Extreme Spread = max(sqrt((x_i - x_j)^2 - (y_i - y_j)^2)),
 * Output sample moments as well as quantile function in steps of 0.5%.
*/

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics_double.h>

#define ITERATIONS 2000000
#define SHOTS 24
#define SIGMA 1.0

#define SQR(x) ((x)*(x))

int shoot (long n, gsl_rng *r, double x0, double y0, double sigma1, double sigma2, double rho, double *x, double *y);

int main (void){
	/* Initialize GSL random number generator */
	const gsl_rng_type *T;
	gsl_rng *rng;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rng = gsl_rng_alloc(T);

	long i, j, k, n = ITERATIONS;

	unsigned group; // Group size to simulate -- iterates [2 to SHOTS]
	unsigned groups; // Number of groups to average: iterates for all groups [1 to SHOTS/2] where groups * group = SHOTS
	double sigma = SIGMA;
	// Samples for each shot
	double *x, *y;
	// Variables with one value per group
	double *extremeSpread;
	// Variables for averaging over all samples in a group size
	double x_min, x_max, y_min, y_max, radius, ES;
	unsigned left, right, top, bottom; // Index of extreme shots in each group

	extremeSpread = (double *) malloc(ITERATIONS*sizeof(double));

	printf("Iterations=,%ld,Sigma=,%lf\n\n", n, sigma);

	/* ************************** */
	/* Outer loop over group size */

	// Header
	printf("Shots,Groups,Mean,Stdev,Skew,Kurtosis");
	for(i = 1;i<=200;i++)
		printf(",%lf",i/200.0);
	printf("\n");

	for (group = 2; group <= SHOTS; group += 1){
		if(SHOTS % group > 0) continue; // Only run scenarios where group * groups = SHOTS exactly (i.e., we use exactly all our SHOTS in the measure)
		groups = SHOTS / group;

		// Run N iterations of of the average of "groups" groups of "group" shots
		for (j = 0; j < n; j++){
			extremeSpread[j] = 0;
			for(k = 0;k < groups;k++){
				// For each iteration:
				x = (double *) malloc(group*sizeof(double));
				y = (double *) malloc(group*sizeof(double));

				shoot(group, rng, 0, 0, sigma, sigma, 0, x, y);

				// Compute group center and range
				left = right = top = bottom = 0;
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

					if(y[i] < y_min){
						y_min = y[i];
						bottom = i;
					}
					if(y[i] > y_max){
						y_max = y[i];
						top = i;
					}
				}
				// Compute extreme spread by checking distance between the four extreme shots
				ES = SQR(x[top] - x[bottom]) + SQR(y[top] - y[bottom]);
				radius = SQR(x[left] - x[right]) + SQR(y[left] - y[right]);
				if(radius > ES) ES = radius;
				radius = SQR(x[top] - x[right]) + SQR(y[top] - y[right]);
				if(radius > ES) ES = radius;
				radius = SQR(x[left] - x[top]) + SQR(y[left] - y[top]);
				if(radius > ES) ES = radius;
				radius = SQR(x[left] - x[bottom]) + SQR(y[left] - y[bottom]);
				if(radius > ES) ES = radius;
				radius = SQR(x[bottom] - x[right]) + SQR(y[bottom] - y[right]);
				if(radius > ES) ES = radius;
				extremeSpread[j] += sqrt(ES);
			}
			extremeSpread[j] /= groups;
		}
	
		gsl_sort(extremeSpread, 1, n);

		printf("%ld,", group);
		printf("%ld,", groups);
		printf("%lf,", gsl_stats_mean(extremeSpread, 1, n));
		printf("%lf,", gsl_stats_sd(extremeSpread, 1, n));
		printf("%lf,", gsl_stats_skew(extremeSpread, 1, n));
		printf("%lf", gsl_stats_kurtosis(extremeSpread, 1, n));
		for(i = 1;i <= 200;i++)
			printf(",%lf", extremeSpread[(int)(i*n/200.0)-1]);
		printf("\n");

		free(x);
		free(y);
	}

	gsl_rng_free(rng);

	free(extremeSpread);

	return 0;
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
//		y[i] += y0;
	}
	return 0;
}

/* ******************************* */      
