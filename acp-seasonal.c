#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mio.h"
#include "meanvar.h"
#include "poisson.h"

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
\t-M monte-carlo iterations\n\
\t-Z <use crude zero inflation compensation>\n\
\t-V <verbose mode>\n";

#define ARGS "f:l:L:VP:I:R:M:Z"


char Fname[255];
int Verbose;
int Sample_count;
int ARLags;
int LLags;
int Period;
double Rate;
int Iterations;
int MC;
int Zero_compensate;

double ACPLogLike(double *data, int fields, int f, int t, double *lam, double omega,
                              double *alphas, int arlags, 
			      double *betas, int llags, int period)
{
	double i;
	double lgnfac = 0.0;
	double ll;
	double y_t;
	double sum;

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

double TotLogLike(double *data, int fields, int f, int recs, double *lam, double omega,
                              double *alphas, int arlags, 
			      double *betas, int llags, int period)
{
	int i;
	double ll;
	double y_t;
	double al;

	for(i=(arlags+llags); i < recs; i++) {
	}
	ll = 0;
	for(i=(arlags+llags); i < recs; i++) {
		al = ACPLogLike(data,fields,f,i,lam,omega,alphas,arlags,betas,llags,period);
		ll += al;

	}

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
 *
 * format: [omega, alpha1, alpha2...alpha(p1), alpha(p2)...beta1, * beta1...betap1, betap2]
 */
double *Partialmu_t(double *data, int fields, int f, int t, 
			double omega, double *alphas, int acount, double *betas,
			int bcount, double *mu_history, double **partial_hist, int period)
{
	double *part;
	double *prev_part;
	int p;
	int i;
	double sum;
	int k;
	int j;

	if(period > 0) {
		part = (double *)malloc(2*(acount+bcount+1)* sizeof(double));
	} else {
		part = (double *)malloc((acount+bcount+1)* sizeof(double));
	}
	if(part == NULL) {
		exit(1);
	}

	/*
	 * compute the sum of previous partials
	 *
	 * first element in partial is omega
	 */
	sum = 0;
	if(period > 0) {
		for(p = 0; p < bcount; p++) {
			prev_part = partial_hist[p];
			sum += (prev_part[0] * betas[p]);
		}
		k = 0;
		for(p = bcount; p < 2*bcount; p++) {
			prev_part = partial_hist[period-k-1];
			sum += (prev_part[0] * betas[p]);
			k++;
		}
	} else {
		for(p = 0; p < bcount; p++) {
			prev_part = partial_hist[p];
			sum += (prev_part[0] * betas[p]);
		}
	}

	/*
	 * 0 is the omega term
	 */
	sum += 1;
	part[0] = sum;

	/*
	 * now do the alphas (after omgea in partials array)
	 */
	if(period > 0) {
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
			sum += data[(t-i)*fields+f];
			part[i+1] = sum;
		}
		for(i=acount; i < 2*acount; i++) {
			sum = 0;
			/*
			 * previous partials
			 */
			k = 0;
			for(p = bcount; p < 2*bcount; p++) {
				prev_part = partial_hist[period-k-1];
				sum += (prev_part[i+1] * betas[p]);
				k++;
			}
			/*
			 * current partial
			 */
			sum += data[(t-period-i)*fields+f];
			part[i+1] = sum;
		}
	} else {
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
			sum += data[(t-i)*fields+f];
			part[i+1] = sum;
		}
	}

	/*
	 * now do the betas
	 */
	if(period > 0) {
		for(i=0; i < bcount; i++) {
			sum = 0;
			for(p = 0; p < bcount; p++) {
				prev_part = partial_hist[p];
				//sum += (prev_part[i+1+acount+1] * betas[p]);
				sum += (prev_part[2*acount+1+i] * betas[p]);
			}
			sum += mu_history[i]; /* first element is mu_(t-1) */
			//part[i+1+acount+1] = sum;
			part[2*acount+1+i] = sum;
		}
		k = 0;
		j = 0;
		for(i=bcount; i < 2*bcount; i++) {
			sum = 0;
			for(p = bcount; p < 2*bcount; p++) {
				prev_part = partial_hist[period-k-1];
				//sum += (prev_part[i+1+acount+1] * betas[p]);
				sum += (prev_part[2*acount+1+i] * betas[p]);
				k++;
			}
			sum += mu_history[period-j-1]; /* first element is mu_(t-1) */
			//part[i+1+acount+1] = sum;
			part[2*acount+1+i] = sum;
			j++;
		}
	} else {
		for(i=0; i < bcount; i++) {
			sum = 0;
			for(p = 0; p < bcount; p++) {
				prev_part = partial_hist[p];
				sum += (prev_part[i+1+acount] * betas[p]);
			}
			sum += mu_history[i]; /* first element is mu_(t-1) */
			part[i+1+acount] = sum;
		}
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
		double *lam_history, double **part_hist, int period)
{
	double *new;
	int i;
	double *part;

	/*
	 * get partials for u_t
	 */
	part = Partialmu_t(data,fields,f,t,
			   omega,alphas,acount,betas,bcount,
			   lam_history,part_hist,period);
	if(part == NULL) {
		exit(1);
	}

	if(period > 0) {
		new = (double *)malloc(2*(acount+bcount+1) * sizeof(double));
	} else {
		new = (double *)malloc((acount+bcount+1) * sizeof(double));
	}
	if(new == NULL) {
		exit(1);
	}

	if(period > 0) {
/*
		for(i=0; i < (acount+bcount+1); i++) {
			new[i] = ((data[t*fields+f] - lam[t]) / lam[t]) * part[i];
		}
		for(i=acount+bcount+1; i < 2*(acount+bcount)+1; i++) {
			new[i] = ((data[(t-period)*fields+f] - lam[t-period]) / lam[t-period]) * part[i];
		}
*/
		/*
		 * omega in [0] and the alphas for time t
		 */
		for(i=0; i < (acount+1); i++) {
			new[i] = ((data[t*fields+f] - lam[t]) / lam[t]) * part[i];
		}
		/*
		 * do the alphas for the period
		 */
		for(i=(acount+1); i < 2*acount+1; i++) {
			new[i] = ((data[(t-period)*fields+f] - lam[t-period]) / lam[t-period]) * part[i];
		}
		/*
		 * do the betas for t
		 */
		for(i=2*acount+1; i < (2*acount+1+bcount); i++) {
			new[i] = ((data[t*fields+f] - lam[t]) / lam[t]) * part[i];
		}
		/*
		 * to the betas for the period
		 */
		for(i=(2*acount+1+bcount); i < (2*(acount+bcount)+1); i++) {
			new[i] = ((data[t*fields+f] - lam[t]) / lam[t]) * part[i];
		}
	} else {
		for(i=0; i < (acount+bcount+1); i++) {
			new[i] = ((data[t*fields+f] - lam[t]) / lam[t]) * part[i];
		}
	}

	free(part);

	return(new);
}

int ChooseTheta(double *alphas, int acount, double *betas, int bcount, int period)
{
	int i;
	double sum;
	int done = 0;
	double rest;

	while(!done) {
		rest = 1.0;
		sum = 0;
		if(period > 0) {
			for(i=0; i < 2*acount; i++) {
				alphas[i] = drand48() * drand48() * rest;
				sum += (alphas[i] * alphas[i]);
				rest = rest - sum;
				
			}
			for(i=0; i < 2*bcount; i++) {
				betas[i] = drand48() * drand48() * rest;
				sum += (betas[i] * betas[i]);
				rest = rest - sum;
			}
		} else {
			for(i=0; i < acount; i++) {
				alphas[i] = drand48() * drand48() * rest;
				sum += (alphas[i] * alphas[i]);
				rest = rest - sum;
				
			}
			for(i=0; i < bcount; i++) {
				betas[i] = drand48() * drand48() * rest;
				sum += (betas[i] * betas[i]);
				rest = rest - sum;
			}
		}
		if(sum < 1) {
			done = 1;
		}
	}

	return(0);
}
	

	
int main(int argc, char *argv[])
{
	int c;
	int err;
	int i;
	int j;
	int t;
	int f;
	int k;
	int m;
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
	double totll;
	double old_ll;
	double new_ll;
	double *max_a;
	double *max_b;
	double max_o;
	double max_ll;
	int y_t;
	double *history;
	double phi;
	int nan_done;

	memset(Fname,0,sizeof(Fname));
	ARLags = 0;
	LLags = 0;
	Iterations = 1;
	Rate = 0.01;
	MC = 1;
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
			case 'M':
				MC = atoi(optarg);
				break;
			case 'Z':
				Zero_compensate = 1;
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

	if(Period != 0) {
		if((Period < ARLags) || (Period < LLags)) {
			fprintf(stderr,
		"period specified must be bigger than AR and cond mu lags\n"
);
			fprintf(stderr,"%s",Usage);
			exit(1);
		}
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

	if(Period > 0) {
		alphas = (double *)malloc(2*ARLags * sizeof(double));
	} else {
		alphas = (double *)malloc(ARLags * sizeof(double));
	}
		
	if(alphas == NULL) {
		exit(1);
	}

	if(Period > 0) {
		betas = (double *)malloc(2*LLags * sizeof(double));
	} else {
		betas = (double *)malloc(LLags * sizeof(double));
	}
	if(betas == NULL) {
		exit(1);
	}

	lam = (double *)malloc(recs * sizeof(double));
	if(lam == NULL) {
		exit(1);
	}

	/*
	 * initialize first set of lambdas with mu
	 */
	if(Period > 0) {
		start = Period;
	} else if(ARLags > LLags) {
		start = ARLags;
	} else {
		start = LLags;
	}
	for(i=0; i < start; i++) {
		lam[i] = mu;
	}

	if(Period > 0) {
		totgrad = (double *)malloc(2*(ARLags+LLags+1) * sizeof(double));
		lam_history = (double *)malloc(Period * sizeof(double));
		part_history = (double **)malloc(Period * sizeof(double *));
	} else {
		totgrad = (double *)malloc((ARLags+LLags+1) * sizeof(double));
		lam_history = (double *)malloc(LLags * sizeof(double));
		part_history = (double **)malloc(LLags * sizeof(double *));
	}
	if(totgrad == NULL) {
		exit(1);
	}
	if(lam_history == NULL) {
		exit(1);
	}
	if(part_history == NULL) {
		exit(1);
	}
	if(Period > 0) {
		for(i=0; i < Period; i++) {
			part_history[i] = (double *)malloc(2*(ARLags+LLags+1)*sizeof(double));
			if(part_history[i] == NULL) {
				exit(1);
			}
			memset(part_history[i],0,2*(ARLags+LLags+1)*sizeof(double));
		}
	} else {
		for(i=0; i < LLags; i++) {
			part_history[i] = (double *)malloc((ARLags+LLags+1)*sizeof(double));
			if(part_history == NULL) {
				exit(1);
			}
			memset(part_history[i],0,(ARLags+LLags+1)*sizeof(double));
		}
	}

	err = MeanVar(data_mio,f,&mu,&var);
	if(err < 0) {
		exit(1);
	}

	if(Period > 0) {
		max_a = (double *)malloc(2*ARLags * sizeof(double));
		max_b = (double *)malloc(2*LLags * sizeof(double));
	} else {
		max_a = (double *)malloc(ARLags * sizeof(double));
		max_b = (double *)malloc(LLags * sizeof(double));
	}
	if(max_a == NULL) {
		exit(1);
	}
	if(max_b == NULL) {
		exit(1);
	}

	max_ll = -9999999999999.9;

	for(m=0; m < MC; m++) {
		if(Period > 0) {
			memset(lam_history,0,Period*sizeof(double));
			for(i=0; i < Period; i++) {
				memset(part_history[i],0,2*(ARLags+LLags+1)*sizeof(double));
			}
		} else {
			memset(lam_history,0,LLags*sizeof(double));
			for(i=0; i < LLags; i++) {
				memset(part_history[i],0,(ARLags+LLags+1)*sizeof(double));
			}
		}
		err = ChooseTheta(alphas,ARLags,betas,LLags,Period);
		if(err < 0) {
			exit(1);
		}
		sum = 0;
		if(Period > 0) {
			for(i=0; i < 2*ARLags; i++) {
				sum += alphas[i];
			}
			for(i=0; i < 2*LLags; i++) {
				sum += betas[i];
			}
		} else {
			for(i=0; i < ARLags; i++) {
				sum += alphas[i];
			}
			for(i=0; i < LLags; i++) {
				sum += betas[i];
			}
		}
		/*
		 * unconditional expectation mu = omega / (1 - sum(alphas and betas)
		 */
		omega = mu * (1 - sum);
		
		nan_done = 0;
		for(j=0; j < Iterations; j++) {
			ll = 0;
			totll = 0;
			totalgrad = 0;
			if(Period > 0) {
				memset(totgrad,0,2*(ARLags+LLags+1)*sizeof(double));
				for(i=0; i < Period; i++) {
					memset(part_history[i],0,2*(ARLags+LLags+1)*sizeof(double));
				}
				memset(lam_history,0,Period * sizeof(double));
				for(i=0; i < Period; i++) {
					lam_history[i] = mu;
				}
			} else {
				memset(totgrad,0,(ARLags+LLags+1)*sizeof(double));
				for(i=0; i < LLags; i++) {
					memset(part_history[i],0,(ARLags+LLags+1)*sizeof(double));
				}
				memset(lam_history,0,LLags * sizeof(double));
				for(i=0; i < LLags; i++) {
					lam_history[i] = mu;
				}
			}
			
			for(t=0; t < (start+ARLags+LLags); t++) {
				lam[t] = mu;
			}
			for(t=start+ARLags+LLags; t < recs; t++) {
	//			if(data[t*data_fields+f] == 0) {
	//				continue;
	//			}
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
				if(Period > 0) {
					for(i=0; i < ARLags; i++) {
						lam[t] += (alphas[i+ARLags] * data[(t-Period-i)*data_fields+f]);
					}
					for(i=0; i < LLags; i++) {
						lam[t] += (betas[i+LLags] * lam_history[Period-i-1]);
					}
				}
				grad = ACPGrad(data,data_fields,f,t,lam,
						omega, alphas,ARLags, betas, LLags,
						lam_history,part_history,Period);
				if(Period > 0) {
					for(i=0; i < 2*(ARLags+LLags+1); i++) {
						totgrad[i] += grad[i];
					}
				} else {
					for(i=0; i < (ARLags+LLags+1); i++) {
						totgrad[i] += grad[i];
					}
				}
				totalgrad++;
				ll += ACPLogLike(data,data_fields,f,t,lam,
						omega,alphas,ARLags,betas,LLags,Period); 
//printf("ll[%d]: %f\n",t,ll);
				totll++;
				plam = Partialmu_t(data,data_fields,f,t,
						   omega,alphas,ARLags,betas,LLags,
						   lam_history,part_history,Period);
				if(plam == NULL) {
					exit(1);
				}
				free(grad);
				/*
				 * age the histories
				 * most recent values are in [0] location
				 * so we age off the end
				 */
				if(Period > 0) {
					free(part_history[Period-1]);
					for(i=Period-2; i >= 0; i--) {
						lam_history[i+1] = lam_history[i];
						part_history[i+1] = part_history[i];
					}
				} else {
					free(part_history[LLags-1]);
					for(i=LLags-2; i >= 0; i--) {
						lam_history[i+1] = lam_history[i];
						part_history[i+1] = part_history[i];
					}
				}
				lam_history[0] = lam[t];
				part_history[0] = plam;
			}


			/*
			 * update the parameters using LL and average gradient
			 */
			ll = TotLogLike(data,data_fields,f,t,lam,
					omega,alphas,ARLags,betas,LLags,Period); 

			omega = omega + (Rate * (totgrad[0] / totalgrad));

			if(isnan(omega)) {
				nan_done = 1;
			}

			if(Period > 0) {
				for(i=0; i < 2*ARLags; i++) {
					alphas[i] = alphas[i] + (Rate*(totgrad[i+1]/totalgrad));
					if(isnan(alphas[i])) {
						nan_done = 1;
					}
				}
				for(i=0; i < 2*LLags; i++) {
					betas[i] = betas[i] + (Rate*(totgrad[i+(2*ARLags)+1]/totalgrad));
					if(isnan(betas[i])) {
						nan_done = 1;
					}
				}
			} else {
				for(i=0; i < ARLags; i++) {
					alphas[i] = alphas[i] + (Rate*(totgrad[i+1]/totalgrad));
					if(isnan(alphas[i])) {
						nan_done = 1;
					}
				}
				for(i=0; i < LLags; i++) {
					betas[i] = betas[i] + (Rate*(totgrad[i+ARLags+1]/totalgrad));
					if(isnan(betas[i])) {
						nan_done = 1;
					}
				}
			}

			if(nan_done == 1) {
				break;
			}

			
			ll = TotLogLike(data,data_fields,f,t,lam,
					omega,alphas,ARLags,betas,LLags,Period); 

			if(isnan(ll)) {
				break;
			}

			if(ll > max_ll) {
				max_ll = ll;
				max_o = omega;
				if(Period > 0) {
					for(i=0; i < 2*ARLags; i++) {
						max_a[i] = alphas[i];
					}
					for(i=0; i < 2*LLags; i++) {
						max_b[i] = betas[i];
					}
				} else {
					for(i=0; i < ARLags; i++) {
						max_a[i] = alphas[i];
					}
					for(i=0; i < LLags; i++) {
						max_b[i] = betas[i];
					}
				}
			}

			if(Verbose == 1) {
				printf("[%d] iteration: %d, avg ll: %f\n",m,j,ll);
				printf("\tomega: %f avggrad: %f\n",omega,totgrad[0] / totalgrad);
				if(Period > 0) {
					for(i=0; i < 2*ARLags; i++) {
						printf("\ta[%d]: %f avggrad: %f\n",i,alphas[i], (totgrad[i+1] / totalgrad));
					}
					for(i=0; i < 2*LLags; i++) {
						printf("\tb[%d]: %f avggrad: %f\n",i,betas[i], totgrad[i+(2*ARLags)+1]/totalgrad);
					}
				} else {
					for(i=0; i < ARLags; i++) {
						printf("\ta[%d]: %f avggrad: %f\n",i,alphas[i], (totgrad[i+1] / totalgrad));
					}
					for(i=0; i < LLags; i++) {
						printf("\tb[%d]: %f avggrad: %f\n",i,betas[i], totgrad[i+ARLags+1]/totalgrad);
					}
				}
			}

		}
	}


	if(Verbose == 1) {
		printf("MLE theta: %f\n",max_ll);
		printf("\tomega: %f\n",max_o);
		if(Period > 0) {
			for(i=0; i < 2*ARLags; i++) {
				printf("\ta[%d]: %f\n",i,max_a[i]);
			}
			for(i=0; i < 2*LLags; i++) {
				printf("\tb[%d]: %f\n",i,max_b[i]);
			}
		} else {
			for(i=0; i < ARLags; i++) {
				printf("\ta[%d]: %f\n",i,max_a[i]);
			}
			for(i=0; i < LLags; i++) {
				printf("\tb[%d]: %f\n",i,max_b[i]);
			}
		}
		exit(1);
	}

	/*
	 * generate an artificial series using max_o, max_a and max_b
	 */

	history = (double *)malloc(Period * sizeof(double));
	if(history == NULL) {
		exit(1);
	}

	if(Zero_compensate == 1) {
		/*
		 * compute probability that any value is a zero
		 */
		sum = 0;
		for(i=start; i < recs; i++) {
			if(data[i*data_fields+f] == 0) {
				sum += 1;
			}
		}
		phi = sum / (double)(recs - start);
	}

	for(t = 0; t < start; t++) {
		history[t] = data[t*data_fields+f];
	}

	memset(lam_history,0,Period*sizeof(double));

	for(i=0; i < start; i++) {
		lam_history[i] = mu;
	}

	for(t=start; t < recs; t++) {
		sum = max_o;
		for(i=0; i < ARLags; i++) {
			sum += max_a[i] * history[i];
		}
		j = 0;
		for(i=ARLags; i < 2*ARLags; i++) {
			sum += max_a[i] * history[Period-j-1];
			j++;
		}
		for(i=0; i < LLags; i++) {
			sum += max_b[i] * lam_history[i];
		}
		j = 0;
		for(i=LLags; i < 2*LLags; i++) {
			sum += max_b[i] * lam_history[Period-j-1];
			j++;
		}
		if(Zero_compensate == 0) {
			y_t = InvertedPoissonCDF(sum);
		} else {
			if(drand48() < (1.0 - phi)) {
				y_t = InvertedPoissonCDF(sum);
			} else {
				y_t = 0;
			}
		}
		printf("%10.0f %d\n",data[t*data_fields+0],y_t);

		for(i = Period-2; i >= 0; i--) {
			history[i+1] = history[i];
			lam_history[i+1] = lam_history[i];
		}
		history[0] = (double)y_t;
		lam_history[0] = sum;
	}
			 
	free(history);

	if(Period > 0) {
		for(i=0; i < Period; i++) {
			free(part_history[i]);
		}
	} else {
		for(i=0; i < LLags; i++) {
			free(part_history[i]);
		}
	}
	free(lam_history);

	free(part_history);
			
	MIOClose(data_mio);
	free(alphas);
	free(betas);
	free(lam);
	free(max_a);
	free(max_b);

	return(0);
}

