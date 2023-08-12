#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "poisson.h"

/*
 * generates an artificial ACP(1,1) series
 */

char *Usage = "acp-gen -a alpha\n\
\t-b beta\n\
\t-o omega\n\
\t-c count\n\
\t-s seasonal-period (in lags)\n";

#define ARGS "a:b:o:c:s:"

double Alpha;
double Beta;
double Omega;
double Count;
int Period;


int main(int argc, char *argv[])
{
	int c;
	int i;
	int j;
	double mu_t;
	double mu_t_minus_1;
	double n_t;
	double n_t_minus_1;
	double *history;
	

	while((c = getopt(argc,argv,ARGS)) != EOF)
	{
		switch(c)
		{
			case 'a':
				Alpha = atof(optarg);
				break;
			case 'b':
				Beta = atof(optarg);
				break;
			case 'o':
				Omega = atof(optarg);
				break;
			case 'c':
				Count = atoi(optarg);
				break;
			case 's':
				Period = atoi(optarg);
				break;
			default:
				fprintf(stderr,
			"unrecognized command %c\n",(char)c);
				fprintf(stderr,"usage: %s",Usage);
				exit(1);
		}
	}

	if(Alpha == 0) {
		fprintf(stderr,"must specify alpha\n");
		fprintf(stderr,"%s",Usage);
		exit(1);
	}
	if(Beta == 0) {
		fprintf(stderr,"must specify beta\n");
		fprintf(stderr,"%s",Usage);
		exit(1);
	}
	if(Count == 0) {
		fprintf(stderr,"must specify a count\n");
		fprintf(stderr,"%s",Usage);
		exit(1);
	}

	if(Period > 0) {
		history = (double *)malloc(Period * sizeof(double));
		if(history == NULL) {
			exit(1);
		}
	}

	mu_t_minus_1 = 0;
	n_t_minus_1 = 1;
	mu_t = Omega + (Alpha * n_t_minus_1) + (Beta * mu_t_minus_1);
	n_t = InvertedPoissonCDF(mu_t);
	n_t_minus_1 = n_t;
	mu_t_minus_1 = mu_t;
	for(i=0; i < Count; i++) {
		mu_t = Omega + (Alpha * n_t_minus_1) + (Beta * mu_t_minus_1);
		n_t = InvertedPoissonCDF(mu_t);
		if(Period > 0) {
			if(i < Period) {
				history[i] = n_t;
				printf("%d %f\n",i,n_t);
			} else {
				if(drand48() < Alpha) {
					n_t = n_t + history[0];
				}
				for(j=0; j < Period-1; j++) {
					history[j] = history[j+1];
				}
				history[Period-1] = n_t;
				printf("%d %f\n",i,n_t);
			}
		} else {
			printf("%d %f\n",i,n_t);
		}
				
		n_t_minus_1 = n_t;
		mu_t_minus_1 = mu_t;
	}

	exit(0);
}



