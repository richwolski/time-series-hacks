#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mio.h"
#include "poisson.h"
#include "meanvar.h"
#include "yw-estimate.h"
#include "redblack.h"

/*
 * for WTB interarrival times where cameras take sequences
 * when ever they are triggered
 *
 * input time series is 
 *  ts count
 * where ts are evenly spaced and count is the number fo values
 * including zeros, in each interval
 */

char *Usage = "inar-series -f filename\n\
\t-c count\n\
\t-l lags\n\
\t-V <verbose mode>\n";

#define ARGS "f:l:V"

char Fname[255];
int Verbose;
int Sample_count;
int Lags;

	
double EstimateLambda(MIO *d_mio, unsigned long size, int f, 
		int lags, double* alphas)
{
	int i;
	double mu;
	double var;
	double lambda;
	double sum = 0.0;
	int err;

	err = MeanVar(d_mio,f,&mu,&var);
	if(err < 0) {
		fprintf(stderr,"no mean and var for data\n");
		exit(1);
	}

	sum = 1.0;
	for(i=0; i < lags; i++) {
		sum -= alphas[i];
	} 

	lambda = mu * sum;

	return(lambda);
}

int BernoulliThin(int x, double phi)
{
	double r;
	int i;
	int sum = 0;

	if(phi < 0) {
		return(0);
	}

	for(i=0; i < x; i++) {
		r = drand48();
		if(r < phi) {
			sum += 1;
		}
	}
	return(sum);
}


int main(int argc, char *argv[])
{
	int c;
	int err;
	int i;
	int j;
	int k;
	int f;
	MIO *data_mio;
	MIO *raw_mio;
	double *data;
	int data_fields;
	unsigned long size;
	unsigned long recs;
	double *alpahs;
	double lambda;
	int innovation;
	double *history;
	double y_t;
	double *alphas;

	memset(Fname,0,sizeof(Fname));
	Lags = 0;
	while((c = getopt(argc,argv,ARGS)) != EOF)
	{
		switch(c)
		{
			case 'f':
				strncpy(Fname,optarg,sizeof(Fname));
				break;
			case 'V':
				Verbose = 1;
				break;
			case 'c':
				Sample_count = atoi(optarg);
				break;
			case 'l':
				Lags = atoi(optarg);
				break;
			default:
				fprintf(stderr,
			"unrecognized command %c\n",(char)c);
				fprintf(stderr,"usage: %s",Usage);
				exit(1);
		}
	}

	if(Fname[0] == 0)
	{
		fprintf(stderr,"must specify file name\n");
		fprintf(stderr,"usage: %s",Usage);
		exit(1);
	}

	if(Lags == 0) {
		fprintf(stderr,"musty specify lags\n");
		fprintf(stderr,"usage: %s",Usage);
		exit(1);
	}


	size = MIOSize(Fname);
	if(size <= 0) {
		fprintf(stderr,"couldn't get size from %s\n",Fname);
		exit(1);
	}

	raw_mio = MIOOpenText(Fname,"r",size);
	if(raw_mio == NULL) {
		fprintf(stderr,"couldn't make mio for %s\n",Fname);
		exit(1);
	}

	data_mio = MIODoubleFromText(raw_mio,NULL);
	if(data_mio == NULL) {
		fprintf(stderr,"couldn't make data mio for %s\n",Fname);
		exit(1);
	}

	data_fields = MIOTextFields(raw_mio); /* data in last field */
	f = data_fields - 1;
	recs = MIOTextRecords(raw_mio);

	MIOClose(raw_mio);

	data = MIOAddr(data_mio);

	alphas = YWEstimate(data_mio,data_fields-1,Lags);
	if(alphas == NULL) {
		fprintf(stderr,"YWestimate failed\n");
		exit(1);
	}
	lambda = EstimateLambda(data_mio, recs, data_fields-1, Lags, alphas);


	if(Verbose == 1) {
		printf("CLS lambda: %f\n",lambda);
		for(i=1; i <= Lags; i++) {
			printf("alpha[%d]: %f\n",i,alphas[i-1]);
		}
		exit(1);
	}

	/*
	 * generate an artificial INAR(1)_s series with Poisson innovations
	 *
	 * prime the pump with initial period
	 */
	history = (double *)malloc(Lags * sizeof(double));
	if(history == NULL) {
		exit(1);
	}

	for(i=0; i < Lags; i++) {
		printf("%10.0f %d\n",data[i*data_fields+0],(int)data[i*data_fields+f]);
		history[i] = data[i*data_fields+f];
	}


	for(i=Lags; i < recs; i++) {
		innovation = InvertedPoissonCDF(lambda);
		y_t = 0.0;
		for(j=0; j < Lags; j++) {
			y_t += BernoulliThin(history[Lags-j],alphas[j]);
		}
		y_t += innovation;
		printf("%10.0f %d\n",data[i*data_fields+0],(int)y_t);
		for(j=1; j < Lags; j++) {
			history[j-1] = history[j];
		}
		history[Lags - 1] = y_t;
	}
		

	MIOClose(data_mio);
	free(alphas);

	return(0);
}

	
