/* 
compile line:
gcc -o SymBiSampleRanges SymBiSampleRanges.c -lgsl -lgslcblas -lm
*/
/* Simulate shots from a symmetric bivariate to see mean, -47.5%, -40%, -25%, +25%, +40%, +47.5% values of:
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
#include <sys/time.h>

#define ITERATIONS 1000000
#define MAX_N 100 // Max group size to simulate
#define SIGMA 1.0

#define SQR(x) ((x)*(x))

int shoot (long n, gsl_rng *r, double sigma1, double sigma2, double rho, double *x, double *y);

int main (void) {
	/* Initialize GSL random number generator */
	const gsl_rng_type *T;
	gsl_rng *rng;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	rng = gsl_rng_alloc(T);
	//gsl_rng_set(rng, random_seed());  // Random seed

	long h, i, j, n = ITERATIONS;
	unsigned group; // Group size to simulate -- iterates [2 to MAX_N]
	double sigma = SIGMA;
	// Samples for each shot
	double *x, *y;
	// Variables with one value per group
	double *diagonal, *FoM, *extremeSpread;
	// Variables for averaging over all samples in a group size
	double x_min, x_max, x_bar, x_range, y_min, y_max, y_bar, y_range, radius, spread, ES;

	diagonal = (double *) malloc(ITERATIONS*sizeof(double));
	FoM = (double *) malloc(ITERATIONS*sizeof(double));
	extremeSpread = (double *) malloc(ITERATIONS*sizeof(double));

	printf("Iterations=,%ld,Sigma=,%lf\n\n", n, sigma);

	/* ************************** */
	/* Outer loop over group size */

	// Header for columns
	printf("Group Size,");
	printf("ES -47.5%,");
	printf("ES -40%,");
	printf("ES -25%,");
	printf("ES Median,");
	printf("ES +25%,");
	printf("ES +40%,");
	printf("ES +47.5%,");
	printf("Diag -47.5%,");
	printf("Diag -40%,");
	printf("Diag -25%,");
	printf("Diagonal Median,");
	printf("Diag +25%,");
	printf("Diag +40%,");
	printf("Diag +47.5%,");
	printf("FoM -47.5%,");
	printf("FoM -40%,");
	printf("FoM -25%,");
	printf("FoM Median,");
	printf("FoM +25%,");
	printf("FoM +40%,");
	printf("FoM +47.5%,");
	printf("Extreme Spread Mean,");
	printf("Extreme Spread Stdev,");
	printf("Extreme Spread Skew,");
	printf("Extreme Spread Kurtosis,");
	printf("Diagonal Mean,");
	printf("Diagonal StDev,");
	printf("Diagonal Skew,");
	printf("Diagonal Kurtosis,");
	printf("FoM Mean,");
	printf("FoM StDev,");
	printf("FoM Skew,");
	printf("FoM Kurtosis");
	printf("\n");
	for (group = 2; group <= MAX_N; group += 1){
		x = (double *) malloc(group*sizeof(double));
		y = (double *) malloc(group*sizeof(double));

		// Run N iterations of group shots
		for (j = 0; j < n; j++) {
			// For each iteration:
			shoot(group, rng, sigma, sigma, 0, x, y);

			// Compute group center and range
			ES = 0;
			x_bar = y_bar = 0;
			x_min = x_max = x[0];
			y_min = y_max = y[0];
			for (i = 0; i < group; i++) {
				if(x[i] < x_min)
					x_min = x[i];
				if(x[i] > x_max)
					x_max = x[i];
				x_bar += x[i];

				if(y[i] < y_min)
					y_min = y[i];
				if(y[i] > y_max)
					y_max = y[i];
				y_bar += y[i];

				// Brute-force extreme spread.
				// (Want to make it faster? Find convex hull, then Shamos's Rotating Calipers algorithm.  E.g., https://code.activestate.com/recipes/117225/)
				for(h = i+1; h < group; h++) {
					spread = SQR(x[h] - x[i]) + SQR(y[h] - y[i]);
					if(spread > ES) ES = spread;
				}
			}
			x_bar /= (double)group;
			y_bar /= (double)group;
			x_range = x_max - x_min;
			y_range = y_max - y_min;
			diagonal[j] = sqrt(SQR(x_range) + SQR(y_range));
			FoM[j] = (x_range + y_range) / 2.0;
			extremeSpread[j] = sqrt(ES);
		}

		gsl_sort(extremeSpread, 1, n);
		gsl_sort(diagonal, 1, n);
		gsl_sort(FoM, 1, n);

		printf("%ld,", group);
		printf("%lf,", extremeSpread[(int)(n/40)]);
		printf("%lf,", extremeSpread[(int)(n/10)]);
		printf("%lf,", extremeSpread[(int)(n/4)]);
		printf("%lf,", extremeSpread[(int)(n/2)]);
		printf("%lf,", extremeSpread[(int)(3*n/4)]);
		printf("%lf,", extremeSpread[(int)(9*n/10)]);
		printf("%lf,", extremeSpread[(int)(39*n/40)]);
		printf("%lf,", diagonal[(int)(n/40)]);
		printf("%lf,", diagonal[(int)(n/10)]);
		printf("%lf,", diagonal[(int)(n/4)]);
		printf("%lf,", diagonal[(int)(n/2)]);
		printf("%lf,", diagonal[(int)(3*n/4)]);
		printf("%lf,", diagonal[(int)(9*n/10)]);
		printf("%lf,", diagonal[(int)(39*n/40)]);
		printf("%lf,", FoM[(int)(n/40)]);
		printf("%lf,", FoM[(int)(n/10)]);
		printf("%lf,", FoM[(int)(n/4)]);
		printf("%lf,", FoM[(int)(n/2)]);
		printf("%lf,", FoM[(int)(3*n/4)]);
		printf("%lf,", FoM[(int)(9*n/10)]);
		printf("%lf,", FoM[(int)(39*n/40)]);
		printf("%lf,", gsl_stats_mean(extremeSpread, 1, n));
		printf("%lf,", gsl_stats_sd(extremeSpread, 1, n));
		printf("%lf,", gsl_stats_skew(extremeSpread, 1, n));
		printf("%lf,", gsl_stats_kurtosis(extremeSpread, 1, n));
		printf("%lf,", gsl_stats_mean(diagonal, 1, n));
		printf("%lf,", gsl_stats_sd(diagonal, 1, n));
		printf("%lf,", gsl_stats_skew(diagonal, 1, n));
		printf("%lf,", gsl_stats_kurtosis(diagonal, 1, n));
		printf("%lf,", gsl_stats_mean(FoM, 1, n));
		printf("%lf,", gsl_stats_sd(FoM, 1, n));
		printf("%lf,", gsl_stats_skew(FoM, 1, n));
		printf("%lf", gsl_stats_kurtosis(FoM, 1, n));
		printf("\n");

		free(x);
		free(y);
	}

	gsl_rng_free(rng);

	free(diagonal);
	free(FoM);
	free(extremeSpread);

	return 0;
}


/*****************************
* Generate a set of shots 
* n - number of shots to fire
* r - pointer to random number generator
* (sigma1, sigma2) - Standard deviations in the x and y directions
* rho - covariance
* (*x, *y) - the array containing the shots - have to allocated externally
***************************** */
int shoot (long n, gsl_rng *r, double sigma1, double sigma2, double rho, double *x, double *y){
	long i;
	for (i = 0; i < n; i++){
		gsl_ran_bivariate_gaussian (r, sigma1, sigma2, rho, &x[i], &y[i]);
	}
	return 0;
}

unsigned long int random_seed() {
	struct timeval tv;
	gettimeofday(&tv,0);
	return (tv.tv_sec + tv.tv_usec);
}
