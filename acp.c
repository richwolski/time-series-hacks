#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mio.h"
#include "meanvar.h"

/*
 * for WTB interarrival times where cameras take sequences
 * when ever they are triggered
 *
 * input time series is 
 *  ts count
 * where ts are evenly spaced and count is the number fo values
 * including zeros, in each interval
 */

char *Usage = "acp -f filename\n\
\t-I iterations (for gradient descent)\n\
\t-R rate (learning rate)\n\
\t-l ar lags\n\
\t-L lambda lags\n\
\t-P period (for seasonal case)\n\
\t-V <verbose mode>\n";

#define ARGS "f:l:L:VP:I:R:"

char Fname[255];
int Verbose;
int Sample_count;
int ARLags;
int LLags;
int Period;
double Rate;
int Iterations;

double ACPLogLike(double *data, int fields, int f, int t, double *lam, double omega,
                              double *alphas, int arlags, 
			      double *betas, int llags)
{
	int i;
	double sum = 0.0;
	double lgnfac = 0.0;
	double ll;
	double y_t;

	sum = omega;
	/*
	 * compute AR terms
	 */
	for(i=0; i < arlags; i++) {
		sum += (alphas[i] * data[(t-i-1)*fields + f]);
	}

	/*
	 * compute lambda terms
	 */
	for(i=0; i < llags; i++) {
		sum += (betas[i] * lam[t-i-1]);
	}

	/*
	 * update lambda
	 */
	lam[t] = sum;
	if(lam[t] == 0) {
		lam[t] = 0.00001;
	}

	/*
	 * compute log of y_t!
	 */
	y_t = data[t*fields+f];
	sum = 0;
	for(i=1; i <= y_t; i++) {
		sum += log(i);
	}
	lgnfac = sum; 

	/*
	 * lloglike = N_t * log(lam_t) - lam[t] - log(N_t!)
	 */
	ll = (data[t*fields+f] * log(lam[t])) - lam[t] - lgnfac;
	return(ll);
}

double *ACPGrad(double *data, int fields, int f, int t, double *lam, double ll, 
		double *coeff, int count, double rate)
{
	double *new;
	int i;

	new = (double *)malloc(count * sizeof(double));
	if(new == NULL) {
		exit(1);
	}

	for(i=0; i < count; i++) {
		new[i] = (data[t*fields+f] - lam[t]) / lam[t];
	}

	return(new);
}
		
	
int main(int argc, char *argv[])
{
	int c;
	int err;
	int i;
	int j;
	int t;
	int f;
	MIO *data_mio;
	MIO *raw_mio;
	double *data;
	int data_fields;
	unsigned long size;
	unsigned long recs;
	double *alphas;
	double *totalphagrad;
	double *newalphagrad;
	double *betas;
	double *newbetagrad;
	double *totbetagrad;
	double totalgrad;
	double mu;
	double var;
	double omega;
	double *lam;
	int start;
	double ll;
	double sum;

	memset(Fname,0,sizeof(Fname));
	ARLags = 0;
	LLags = 0;
	Iterations = 1;
	Rate = 0.01;
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
				ARLags = atoi(optarg);
				break;
			case 'L':
				LLags = atoi(optarg);
				break;
			case 'P':
				Period = atoi(optarg);
				break;
			case 'I':
				Iterations = atoi(optarg);
				break;
			case 'R':
				Rate = atof(optarg);
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

	if(ARLags == 0) {
		fprintf(stderr,"musty specify ar lags\n");
		fprintf(stderr,"usage: %s",Usage);
		exit(1);
	}

	if(LLags == 0) {
		fprintf(stderr,"musty specify lambda lags\n");
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

	/*
	 * generate an artificial ACP series with Poisson innovations
	 *
	 * prime the pump with initial period
	 */

	/*
	 * compute unconditional mean of ACP(p,q) to use as initial
	 * lambdas
	 *
	 * make sure initial alphas and betas sume to less than 1
	 */

	alphas = (double *)malloc(ARLags * sizeof(double));
	if(alphas == NULL) {
		exit(1);
	}
	betas = (double *)malloc(LLags * sizeof(double));
	if(betas == NULL) {
		exit(1);
	}

	sum = 0;
	for(i=0; i < ARLags; i++) {
		alphas[i] = 0.5 / (double)(ARLags+LLags);
		sum += alphas[i];
	}

	for(i=0; i < LLags; i++) {
		betas[i] = 0.5 / (double)(ARLags+LLags);
		sum += betas[i];
	}

	err = MeanVar(data_mio,f,&mu,&var);
	if(err < 0) {
		exit(1);
	}

	/*
	 * unconditional expectation mu = omega / (1 - sum(alphas and betas)
	 */
	omega = mu * (1 - sum);
	
	lam = (double *)malloc(recs * sizeof(double));
	if(lam == NULL) {
		exit(1);
	}

	/*
	 * initialize first set of lambdas with mu
	 */
	if(ARLags > LLags) {
		start = ARLags;
	} else {
		start = LLags;
	}
	for(i=0; i < start; i++) {
		lam[i] = mu;
	}

	totalphagrad = (double *)malloc(ARLags * sizeof(double));
	if(totalphagrad == NULL) {
		exit(1);
	}
	totbetagrad = (double *)malloc(LLags * sizeof(double));
	if(totbetagrad == NULL) {
		exit(1);
	}
	for(j = 0; j < Iterations; j++) {
		totalgrad = 0;
		memset(totalphagrad,0,ARLags*sizeof(double));
		memset(totbetagrad,0,LLags*sizeof(double));
		for(t=start; t < recs; t++) {
			ll = ACPLogLike(data,data_fields,f,t,lam,
					omega,alphas,ARLags,betas,LLags); 
			newalphagrad = ACPGrad(data,data_fields,f,t,lam,ll,alphas,ARLags,Rate);
			newbetagrad = ACPGrad(data,data_fields,f,t,lam,ll,betas,LLags,Rate);
			for(i=0; i < ARLags; i++) {
				totalphagrad[i] += newalphagrad[i];
			}
			for(i=0; i < LLags; i++) {
				totbetagrad[i] += newbetagrad[i];
			}
			totalgrad++;
			free(newalphagrad);
			free(newbetagrad);
		}
		for(i=0; i < ARLags; i++) {
			alphas[i] = alphas[i] - (Rate*(totalphagrad[i]/totalgrad)*ll);
		}
		for(i=0; i < LLags; i++) {
			betas[i] = betas[i] - (Rate*(totbetagrad[i]/totalgrad)*ll);
		}
		sum = 0;
		for(i=0; i < ARLags; i++) {
			sum += alphas[i];
		}
		for(i=0; i < LLags; i++) {
			sum += betas[i];
		}
		omega = mu * (1 - sum); 
	}

	if(Verbose == 1) {
		printf("alphas:\n");
		for(i=0; i < ARLags; i++) {
			printf("\t%f\n",alphas[i]);
		}
		printf("betas:\n");
		for(i=0; i < LLags; i++) {
			printf("\t%f\n",betas[i]);
		}
		exit(1);
	}
			
	MIOClose(data_mio);
	free(alphas);
	free(betas);
	free(lam);

	return(0);
}

