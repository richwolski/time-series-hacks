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
	double lgnfac = 0.0;
	double ll;
	double y_t;

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

double OldACPLogLike(double *data, int fields, int f, int t, double *lam, double omega,
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

/*
 * returns a verctor of partials for omega, the alphas, and the betas
 * 
 * the partials all use the beta weighted sume of previous partials (each of
 * which is a vector)
 *
 * the partials also requite previous data values (for the alphas) and
 * previous u_t values (for the betas)
 */
double *Partialmu_t(double *data, int fields, int f, int t, double omega, double *alphas, int acount, double *betas,
			int bcount, double *mu_history, double **partial_hist)
{
	double *part;
	double *prev_part;
	int p;
	int i;
	double sum;

	part = (double *)malloc((acount+bcount+1)* sizeof(double));
	if(part == NULL) {
		exit(1);
	}

	
	/*
	 * compute the sum of previous partials
	 *
	 * first element in partial is omega
	 */
	sum = 0;
	for(p = 0; p < bcount; p++) {
		prev_part = partial_hist[p];
		sum += (prev_part[0] * betas[p]);
	}

	/*
	 * 0 is the omega term
	 */
	sum += 1;
	part[0] = sum;

	/*
	 * now do the alphas (after omgea in partials array)
	 */
	for(i=0; i < acount; i++) {
		sum = 0;
		/*
		 * previous partials
		 */
		for(p = 0; p < bcount; p++) {
			prev_part = partial_hist[p];
			sum += (prev_part[i+1] * betas[p]);
		}
		/*
		 * current partial
		 */
		sum += data[t*fields+f - i - 1];
		part[i+1] = sum;
	}

	/*
	 * now do the betas
	 */
	for(i=0; i < bcount; i++) {
		sum = 0;
		for(p = 0; p < bcount; p++) {
			prev_part = partial_hist[p];
			sum += (prev_part[i+1+acount] * betas[p]);
		}
		sum += mu_history[i]; /* first element is mu_(t-1) */
		part[i+1+acount] = sum;
	}

	return(part);
}


	
/*
 * computes the partial of ll with respect to parameters omega, alphas, and
 * betas
 *
 * requires the partial of mu_t at t which requires previous partials of mu_t
 * and previous values of mu
 */

double *ACPGrad(double *data, int fields, int f, int t, double *lam,
		double omega, double *alphas, int acount, double *betas, int bcount, 
		double *lam_history, double **part_hist)
{
	double *new;
	int i;
	double *part;

	/*
	 * get partials for u_t
	 */
	part = Partialmu_t(data,fields,f,t,
			   omega,alphas,acount,betas,bcount,
			   lam_history,part_hist);
	if(part == NULL) {
		exit(1);
	}

	new = (double *)malloc((acount+bcount+1) * sizeof(double));
	if(new == NULL) {
		exit(1);
	}

	for(i=0; i < (acount+bcount+1); i++) {
		new[i] = ((data[t*fields+f] - lam[t]) / lam[t]) * part[i];
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
	double *totgrad;
	double *betas;
	double *newbetagrad;
	double totalgrad;
	double mu;
	double var;
	double omega;
	double *lam;
	int start;
	double ll;
	double sum;
	double *grad;
	double *lam_history;
	double **part_history;
	double *plam;

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
		alphas[i] = 0.01 / (double)(ARLags+LLags);
		sum += alphas[i];
	}

	for(i=0; i < LLags; i++) {
		betas[i] = 0.01 / (double)(ARLags+LLags);
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

	totgrad = (double *)malloc((ARLags+LLags+1) * sizeof(double));
	if(totgrad == NULL) {
		exit(1);
	}
	lam_history = (double *)malloc(LLags * sizeof(double));
	if(lam_history == NULL) {
		exit(1);
	}
	memset(lam_history,0,LLags * sizeof(double));

	part_history = (double **)malloc(LLags * sizeof(double *));
	if(part_history == NULL) {
		exit(1);
	}
	for(i=0; i < LLags; i++) {
		part_history[i] = (double *)malloc((ARLags+LLags+1)*sizeof(double));
		if(part_history == NULL) {
			exit(1);
		}
		memset(part_history[i],0,(ARLags+LLags+1)*sizeof(double));
	}

	for(j=0; j < Iterations; j++) {
		for(i=0; i < LLags; i++) {
			memset(part_history[i],0,(ARLags+LLags+1)*sizeof(double));
		}
		memset(lam_history,0,LLags * sizeof(double));
		totalgrad = 0;
		memset(totgrad,0,(ARLags+LLags+1)*sizeof(double));
		for(t=start; t < recs; t++) {
			/*
			 * compute current lam_t
			 */
			lam[t] = omega;
			for(i=0; i < ARLags; i++) {
				lam[t] += (alphas[i] * data[(t-i-1)*data_fields+f]);
			}
			for(i=0; i < LLags; i++) {
				lam[t] += (betas[i] * lam_history[i]);
			}
			grad = ACPGrad(data,data_fields,f,t,lam,
					omega, alphas,ARLags, betas, LLags,
					lam_history,part_history);
			for(i=0; i < (ARLags+LLags+1); i++) {
				totgrad[i] += grad[i];
			}
			totalgrad++;
			free(grad);
			plam = Partialmu_t(data,data_fields,f,t,
					   omega,alphas,ARLags,betas,LLags,
					   lam_history,part_history);
			if(plam == NULL) {
				exit(1);
			}
			/*
			 * age the histories
			 * most recent values are in [0] location
			 * so we age off the end
			 */
			free(part_history[LLags-1]);
			for(i=LLags-2; i >= 0; i--) {
				lam_history[i+1] = lam_history[i];
				part_history[i+1] = part_history[i];
			}
			lam_history[0] = lam[t];
			part_history[0] = plam;
		}

		/*
		 * update the parameters using LL and average gradient
	 	 */
		ll = ACPLogLike(data,data_fields,f,t,lam,
					omega,alphas,ARLags,betas,LLags); 
printf("ll: %f\n",ll);

		omega = omega - (Rate * (totgrad[0] / totalgrad)* -1.0 * ll);
		for(i=0; i < ARLags; i++) {
			alphas[i] = alphas[i] - (Rate*(totgrad[i+1]/totalgrad) * -1.0 * ll);
printf("a[%d]: %f, %f\n",i,alphas[i],Rate*(totgrad[i+1]/totalgrad));
		}
		for(i=0; i < LLags; i++) {
			betas[i] = betas[i] - (Rate*(totgrad[i+ARLags+1]/totalgrad) * -1.0 * ll);
printf("b[%d]: %f, %f\n",i,betas[i],Rate*(totgrad[i+ARLags+1]/totalgrad));
		}

		if(Verbose == 1) {
			printf("iteration: %d\n",j);
			printf("\tomega: %f\n",omega);
			for(i=0; i < ARLags; i++) {
				printf("\ta[%d]: %f\n",i,alphas[i]);
			}
			for(i=0; i < LLags; i++) {
				printf("\tb[%d]: %f\n",i,betas[i]);
			}
		}

	}

	for(i=0; i < LLags; i++) {
		free(part_history[i]);
	}
	free(part_history);
	free(lam_history);
			
	MIOClose(data_mio);
	free(alphas);
	free(betas);
	free(lam);

	return(0);
}

