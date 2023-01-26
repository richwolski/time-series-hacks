#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mio.h"
#include "poisson.h"
#include "meanvar.h"
#include "autoc.h"
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
\t-p <use parcor coeff>\n\
\t-P period (for seasonal case)\n\
\t-V <verbose mode>\n";

#define ARGS "f:l:VpP:"

char Fname[255];
int Verbose;
int Sample_count;
int Lags;
int Use_parcor;
int Period;

	
double EstimateLambda(MIO *d_mio, unsigned long size, int f, 
		int lags, int period, double* alphas)
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
		if((period != 0) && (period > lags)) {
			sum -= alphas[period-1-i];
		} 
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
	double *pc;
	double *pc2D;
	int length;

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
			case 'p':
				Use_parcor = 1;
				break;
			case 'l':
				Lags = atoi(optarg);
				break;
			case 'P':
				Period = atoi(optarg);
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

	if((Period != 0) && (Period > Lags)) {
		length = Period;
	} else {
		length = Lags;
	}

	alphas = YWEstimate(data_mio,data_fields-1,length+5);
	if(alphas == NULL) {
		fprintf(stderr,"YWestimate failed\n");
		exit(1);
	}
	lambda = EstimateLambda(data_mio, recs, data_fields-1, Lags, Period,
			alphas);

	if(Use_parcor == 1) {
		pc2D = (double *)malloc((length+1) *(length+1)*sizeof(double));
		if(pc2D == NULL) {
			exit(1);
		}
		pc = (double *)malloc(length * sizeof(double));
		if(pc == NULL) {
			exit(1);
		}
		err = ParCor(data_mio,data_fields-1,length,pc2D);
		if(err < 0) {
			fprintf(stderr,"no parcor\n");
			exit(1);
		}
		for(i=1; i <= Lags; i++) {
			pc[i-1] = pc2D[i*(length+1) + i]; // pc[0] must alpha_1
		}
		free(pc2D);
	}


	if(Verbose == 1) {
		printf("CLS lambda: %f\n",lambda);
		for(i=1; i <= length; i++) {
			printf("alpha[%d]: %f\n",i,alphas[i-1]);
		}
		if(Use_parcor == 1) {
			for(i=1; i <= length; i++) {
				printf("pc[%d]: %f\n",i,pc[i-1]);
			}
		}
			
		exit(1);
	}

	/*
	 * generate an artificial INAR(1)_s series with Poisson innovations
	 *
	 * prime the pump with initial period
	 */
	history = (double *)malloc(length * sizeof(double));
	if(history == NULL) {
		exit(1);
	}

	for(i=0; i < length; i++) {
		printf("%10.0f %d\n",data[i*data_fields+0],(int)data[i*data_fields+f]);
		history[i] = data[i*data_fields+f];
	}


	for(i=length; i < recs; i++) {
		innovation = InvertedPoissonCDF(lambda);
		y_t = 0.0;
		if((Period == 0) || (Lags > Period)) {
			for(j=0; j < Lags; j++) {
				if(Use_parcor == 0) {
					y_t += BernoulliThin(history[Lags-j-1],
							alphas[j]);
				} else {
					y_t += BernoulliThin(history[Lags-j-1],
							pc[j]);
				}
			}
			y_t += innovation;
			printf("%10.0f %d\n",data[i*data_fields+0],(int)y_t);
			for(j=1; j < Lags; j++) {
				history[j-1] = history[j];
			}
			history[Lags - 1] = y_t;
		} else {
			for(j=0; j < Lags; j++) {
				if(Use_parcor == 0) {
					y_t += 
					    BernoulliThin(history[length-j-1],
                                                        alphas[j]);
					y_t += BernoulliThin(history[j],
                                                        alphas[length-j-1]);
				} else {
					y_t += 
					    BernoulliThin(history[length-j-1],
                                                        pc[j]);
					y_t += BernoulliThin(history[j],
                                                        pc[length-j-1]);
				}
			}
			y_t += innovation;
			printf("%10.0f %d\n",data[i*data_fields+0],(int)y_t);
			for(j=1; j < length; j++) {
				history[j-1] = history[j];
			}
			history[length - 1] = y_t;
		}
	}
		

	MIOClose(data_mio);
	free(alphas);
	if(Use_parcor == 1) {
		free(pc);
	}

	return(0);
}

	
