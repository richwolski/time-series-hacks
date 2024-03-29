#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "simple_input.h"
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

char *Usage = "markov-seasonal-series -f filename\n\
\t-B busy_epoch_count\n\
\t-c sample_count\n\
\t-i time interval (for synthetic trace)\n\
\t-I number-of-evdnts-to-add\n\
\t-P cyclic_period (measured in lags)\n\
\t-3 <use 3D method>\n\
\t-V <verbose mode>\n";

#define ARGS "f:Vc:i:P:3I:B:"

char Fname[255];
int Verbose;
int Sample_count;
double Interval;
int Period;
int Use3D;
int DoIncrement;
int BusyCount;

/*
 * this function adds arrivals to the epoch according to the conditional
 * probabilities for the epoch
 *
 * it returns the total number of arrivals it added to the epoch
 * assumes that otcounts and ottables are initialized as out parameters
 */
int IncrementEpoch(int epoch, double **tcounts, double **ttables, int scount,
                   double **otcounts, double **ottables)
{
	double *trans;
	double *counts;
	double *ocounts;
	double *otrans;
	int i;
	int j;
	int k;
	double p;
	double q;
	double r;
	double psum;
	double *cps; /* probabilities of each conditional pair */
	double cptot;
	int incr;
	double scheck;
	double m;

	trans = ttables[epoch]; /* transition counts of next state based on pair */
	counts = tcounts[epoch]; /* count of each conditional pair */
	otrans = ottables[epoch]; /* output version */
	ocounts = otcounts[epoch]; /* output version */

	cps = (double *)malloc(scount*scount*sizeof(double));
	if(cps == NULL) {
		exit(1);
	}
	/*
	 * count up total number of conditional pairs
	 */
	cptot = 0.0;
	for(i=0; i < scount; i++) {
		for(j=0; j < scount; j++) {
			/*
			 * skip cases where all probability is concentrated at
			 * zero
			 */
			p = trans[i*scount*scount+j*scount+0] / counts[i*scount+j];
			if(p == 1.0) {
				continue;
			}
			cptot += counts[i*scount+j];
		}
	}
	/*
	 * sanity check
	 */
	if(cptot == 0.0) {
		free(cps);
		return(-1);
	}

	/*
	 * compute probs of each pair
	 */
	for(i=0; i < scount; i++) {
		for(j=0; j < scount; j++) {
			cps[i*scount+j] = 0.0;
			p = trans[i*scount*scount+j*scount+0] / counts[i*scount+j];
			if(p == 1.0) {
				continue;
			}
			cps[i*scount+j] = counts[i*scount+j]/cptot;
		}
	}


	/*
	 * now for the tricky part
	 *
	 * choose a pair and a count based on the joint prob of a pair and the
	 * count in the transition array
	 * skip the zeros since this is an increment function
	 */
	psum = 0.0;
	r = drand48();
	for(i=0; i < scount; i++) {
		for(j=0; j < scount; j++) {
			/*
			 * skip case where we haven't seen a pair
			 */
			if(counts[i*scount+j] == 0) {
//printf("no counts (%d,%d)\n",i,j);
				continue;
			}
			/*
			 * skip the case where all of the condition prob is
			 * concentrated at zero
			 */
			p = trans[i*scount*scount+j*scount+0] / counts[i*scount+j];
			if(p == 1.0) {
//printf("all zero\n");
				continue;
			}
			scheck = 0.0;
			for(k=0; k < scount; k++) {
				if(k == 0) {
					q = 1-(trans[i*scount*scount+j*scount+k]/counts[i*scount+j]);
					if(q == 0) {
						q = 0.000001;
					}
				} else {
					/* condition on dest not zero */
					p = trans[i*scount*scount+j*scount+k] / counts[i*scount+j];
					m = (cps[i*scount+j]*(p / q));
//if(m > 0) {
//printf("(%d,%d) -- %d p: %f q: %f cps: %f counts: %d t: %d\n",
//i,j,k,p,q,cps[i*scount+j],(int)counts[i*scount+j],(int)trans[i*scount*scount+j*scount+k]);
//}
					psum += m;
//printf("psum: %f\n",psum);
//scheck += (p/q);
//printf("scheck: %f\n",scheck);
					if(r <= psum) {
						/*
						 * found it
						 */
						incr = k;
						/*
						 * update output version
						 */
						otrans[i*scount*scount+j*scount+k] += 1;
						ocounts[i*scount+j] += 1;
						if(Verbose == 1) {
  	printf("incr (%d,%d) -- %d p: %f q: %f psum: %f r: %f\n",
							i,j,k,
							p,q,psum,r);
						}
						free(cps);
						return(incr);
					} else {
//printf("cps: %f p: %f q: %f psum: %f\n",cps[i*scount+j],p,q,psum);
					}
				}
			}
		}
	}

	printf("error: couldn't increment epoch %d\n",epoch);
	scheck = 0.0;
	for(i=0; i < scount; i++) {
		for(j=0; j < scount; j++) {
			printf("(%d,%d) count: %d prob: %f\n",
				i,j,(int)counts[i*scount+j],cps[i*scount+j]);
		}
		printf("r: %f, q: %f, psum: %f\n",r,q,psum);
	}
exit(1);

	free(cps);
	return(-1);
}

int *BusyEpochs(int ecount, int epochs, double **tcounts, double **ttables, int scount)
{
	int *busy_epochs;
	double *trans;
	double *counts;
	double cptot;
	double *cps;
	int i;
	int j;
	int e;
	RB *list;
	RB *rb;
	double p;
	double q;
	Hval hv;

	busy_epochs = (int *)malloc(ecount * sizeof(int));
	if(busy_epochs == NULL) {
		exit(1);
	}
	memset(busy_epochs,0,ecount*sizeof(int));

	cps = (double *)malloc(scount*scount*sizeof(double));
	if(cps == NULL) {
		exit(1);
	}

	list = RBInitD();
	if(list == NULL) {
		exit(1);
	}

	for(e=0; e < epochs; e++) {
		trans = ttables[e]; /* transition counts of next state based on pair */
		counts = tcounts[e]; /* count of each conditional pair */

		/*
	 	 * count up total number of conditional pairs
	 	 */
		cptot = 0.0;
		for(i=0; i < scount; i++) {
			for(j=0; j < scount; j++) {
				cptot += counts[i*scount+j];
			}
		}
		/*
	 	 * sanity check
	 	 */
		if(cptot == 0.0) {
			free(cps);
			free(busy_epochs);
			RBDestroyD(list);
			return(NULL);
		}

		/*
	 	 * compute probs of each pair
	 	 */
		for(i=0; i < scount; i++) {
			for(j=0; j < scount; j++) {
				cps[i*scount+j] = counts[i*scount+j]/cptot;
			}
		}

		p = 0;
		/*
		 * compute total prob of zero in this epoch
		 */
		for(i=0; i < scount; i++) {
			for(j=0; j < scount; j++) {
				if(counts[i*scount+j] == 0) {
					continue;
				}
				p += (cps[i*scount+j] * 
					(trans[i*scount*scount+j*scount+0]/counts[i*scount+j]));
			}
		}
		if(p == 0) {
			continue;
		}
		q = 1 - p; /* prob of non zero */
		hv.i = e;
		RBInsertD(list,q,hv);
	}

	/*
	 * take off the top ecount epochs
	 */
	e = 0;
	RB_BACKWARD(list,rb) {
		busy_epochs[e] = rb->value.i;
		e++;
		if(e >= ecount) {
			break;
		}
	}

	RBDestroyD(list);
	free(cps);

	return(busy_epochs);

}

double **Copyttables(double **ttables, int scount, int epochs)
{
	double **ottables;
	int i;
	int j;
	int k;
	int e;
	double *trans;
	double *otrans;

	ottables = (double **)malloc(epochs * sizeof(double *));
	if(ottables == NULL) {
		fprintf(stderr,"no space for ttables copy\n");
		exit(1);
	}
	for(e=0; e < epochs; e++) {
		otrans = (double *)malloc(scount*scount*scount*sizeof(double));
		if(otrans == NULL) {
			fprintf(stderr,"no space for otrans at %d\n",e);
			exit(1);
		}
		trans = ttables[e];
		for(i=0; i < scount; i++) {
			for(j=0; j < scount; j++) {
				for(k=0; k < scount; k++) {
					otrans[i*scount*scount+j*scount+k] =
					  trans[i*scount*scount+j*scount+k];
				}
			}
		}
		ottables[e] = otrans;
	}

	return(ottables);
}

double **Copytcounts(double **tcounts, int scount, int epochs)
{
	int i;
	int j;
	int e;
	double **otcounts;
	double *counts;
	double *ocounts; 

	otcounts = (double **)malloc(epochs*sizeof(double *));
	if(otcounts == NULL) {
		fprintf(stderr,"no space for otcounts\n");
		exit(1);
	}

	for(e=0; e < epochs; e++) {
		ocounts = (double *)malloc(scount*scount*sizeof(double));
		if(ocounts == NULL) {
			fprintf(stderr,"no space for ocounts at %d\n",e);
			exit(1);
		}
		counts = tcounts[e];
		for(i=0; i < scount; i++) {
			for(j=0; j < scount; j++) {
				ocounts[i*scount+j] = counts[i*scount+j];
			}
		}
		otcounts[e] = ocounts;
	}

	return(otcounts);
}
		

	
						
int main(int argc, char *argv[])
{
	int c;
	int err;
	void *input_data;
	int field_count;
	double *values;
	double *transitions;
	RB *states;
	RB *rb;
	double count;
	Hval hv;
	int state_count;
	int i;
	int j;
	int k;
	int e;
	int src;
	int psrc;
	int dest;
	double max;
	double min;
	double *counts;
	double *history_period;
	double now;
	double r;
	double p;
	double start_ts;
	double **ttables;
	double **tcounts;
	double **ottables;
	double **otcounts;
	int tp;
	int theperiod;
	int lastperiod;
	int incr;
	int added;
	int *busy_epochs;

	memset(Fname,0,sizeof(Fname));
	Interval = 1;
	Period = 0;
	Use3D = 0;
	DoIncrement = 0;
	BusyCount = 1;
	while((c = getopt(argc,argv,ARGS)) != EOF)
	{
		switch(c)
		{
			case 'f':
				strncpy(Fname,optarg,sizeof(Fname));
				break;
			case 'B':
				BusyCount = atoi(optarg);
				break;
			case 'V':
				Verbose = 1;
				break;
			case 'c':
				Sample_count = atoi(optarg);
				break;
			case 'i':
				Interval = atof(optarg);
				break;
			case 'I':
				DoIncrement = atoi(optarg);
				break;
			case 'P':
				Period = atoi(optarg);
				break;
			case '3':
				Use3D = 1;
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

	if(Period == 0) {
		fprintf(stderr,"musty specify period (in lags)\n");
		fprintf(stderr,"usage: %s",Usage);
		exit(1);
	}

	field_count = GetFieldCount(Fname);

	if(field_count <= 0)
	{
		fprintf(stderr,"bad field count in %s\n",Fname);
		exit(1);
	}

	values = (double *)malloc(field_count*sizeof(double));
	if(values == NULL)
	{
		fprintf(stderr,"no space for values\n");
		exit(1);
	}

	err = InitDataSet(&input_data,field_count);
	if(err == 0)
	{
		fprintf(stderr,"couldn't init data for %s\n",Fname);
		fflush(stderr);
		free(values);
		exit(1);
	}

	err = LoadDataSet(Fname,input_data);
	if(err == 0)
	{
		fprintf(stderr,"couldn't load data from %s\n",Fname);
		fflush(stderr);
		free(values);
		exit(1);
	}

	states = RBInitD();
	if(states == NULL) {
		exit(1);
	}


	max = -10000000.0;
	min = 10000000.0;
	state_count = 0;
	start_ts = -1;
	while(ReadData(input_data,field_count,values) != 0) {
		/*
		 * get starting ts for synth trace
		 */
		if(start_ts == -1) {
			start_ts = values[0];
		}
		/*
		 * assume value is in last field
		 */
		rb = RBFindD(states,values[field_count-1]);
		if(rb != NULL) {
			count = rb->value.d + 1.0; 
			hv.d = count;
			RBInsertD(states,values[field_count-1],hv);
		} else {
			if(values[field_count-1] > max) {
				max = values[field_count-1];
			}
			if(values[field_count-1] < min) {
				min = values[field_count-1];
			}
	
			hv.d = 1.0;
			RBInsertD(states,values[field_count-1],hv);
		}
	}

	if(Use3D == 1) {
		ttables = (double **)malloc(Period * sizeof(double *));
		tcounts = (double **)malloc(Period * sizeof(double *));
	} else {
		ttables = (double **)malloc(sizeof(double *));
		tcounts = (double **)malloc(sizeof(double *));
	}

	if(ttables == NULL) {
		exit(1);
	}

	if(tcounts == NULL) {
		exit(1);
	}


	/*
	 * cover all possible states
	 */
	state_count = (int)(max - min) + 1;

	if(Use3D == 1) {
		for(i=0; i < Period; i++) {
			ttables[i] = (double *)malloc(state_count*state_count*state_count*sizeof(double));
			if(ttables[i] == NULL) {
				fprintf(stderr,"no space for ttable %d\n",i);
				exit(1);
			}
			memset(ttables[i],0,state_count*state_count*state_count*sizeof(double));
			tcounts[i] = (double *)malloc(state_count*state_count*sizeof(double));
			if(tcounts[i] == NULL) {
				fprintf(stderr,"no space for tcount %d\n",i);
				exit(1);
			}
			memset(tcounts[i],0,state_count*state_count*sizeof(double));
		}
	} else {

		/*
		 * condition each transition on previous value and period
		 * value
		 */
		ttables[0] = (double *)
			malloc(state_count*state_count*state_count*sizeof(double));
		if(ttables[0] == NULL) {
			exit(1);
		}
		memset(ttables[0],0,state_count*state_count*state_count*sizeof(double));

		tcounts[0] = (double *)malloc(state_count*state_count*sizeof(double));
		if(tcounts[0] == 0) {
			exit(1);
		}
		memset(tcounts[0],0,state_count*state_count*sizeof(double));
	}

	history_period = (double *)malloc(Period * sizeof(double));
	if(history_period == NULL) {
		exit(1);
	}
	/*
	 * note that history_period is initialized to state (0 - min)
	 */
	memset(history_period,0-min,Period*sizeof(double));

	Rewind(input_data);
	/*
	 * compute the transition counts and the count of each state
	 * as a src
	 */
	src = -10000000;
	theperiod = 0;
	while(ReadData(input_data,field_count,values) != 0) {
		psrc = history_period[0]; // get oldest period value
		if(src == -10000000) {
			src = (int)(values[field_count-1] - min);
			history_period[Period-1] = src;
			continue;
		}
		dest = (int)(values[field_count-1]-min);
		transitions = ttables[theperiod];
		transitions[(src*state_count*state_count)+(psrc*state_count)+dest] += 1.0;
		counts = tcounts[theperiod];
		counts[(src*state_count)+psrc] += 1.0;
		/*
		 * add src to history period
		 */
		for(i=0; i < (Period-1); i++) {
			history_period[i] = history_period[i+1]; // age the period 
		}
		history_period[Period - 1] = dest;
		src = dest;
		if(Use3D == 1) {
			theperiod = (theperiod + 1) % Period;
		}
	}

#if 0
	if(Verbose) {
		if(Use3D == 0) {
			lastperiod = 1;
		} else {
			lastperiod = Period;
		}

		for(tp=0; tp < lastperiod; tp++) {
			transitions = ttables[tp];
			counts = tcounts[tp];
			printf("               ");
			for(i=0; i < state_count; i++) {
				printf("%7.7d ",i+(int)min);
			}
			printf("\n");
			for(i=0; i < state_count; i++) {
				for(j=0; j < state_count; j++) {
					printf("[%5.5d %5.5d]: ",i+(int)min,j+(int)min);
					for(k=0; k < state_count; k++) {
						if(counts[state_count*i+j] == 0) {
							printf("%5.5f ",0.0);
						} else {
							p =
							transitions[(state_count*state_count*i)+
								    (state_count*j)+
								    k] / counts[state_count*i+j];
							printf("%5.5f ",p);
							if(p > 1.0) {
								printf("%f %f\n",
									transitions[(state_count*i)+ (state_count*j)+ k],
									counts[state_count*i+j]);
								exit(1);
							}
						}
					}
					printf("\n");
				}
			}
		}
		exit(1);
	}
#endif

	/*
	 * if we are incrementing, do it here
	 */
	if(DoIncrement > 0) {
		if(Use3D == 1) {
			busy_epochs = BusyEpochs(BusyCount,Period,tcounts,ttables,state_count);
		} else {
			busy_epochs = BusyEpochs(BusyCount,1,tcounts,ttables,state_count);
		}
		if(busy_epochs == NULL) {
			printf("failed to get busy epochs\n");
			exit(1);
		}
		if(Verbose) {
			for(i=0; i < BusyCount; i++) {
				printf("busy[%d]: %d\n",i,busy_epochs[i]);
			}
		}
		/*
		 * make a copy that gets updated so that early updates
		 * don't influence later conditional probabilities
		 */
		if(Use3D == 0) {
			ottables = Copyttables(ttables,state_count,1);
			otcounts = Copytcounts(tcounts,state_count,1);
		} else {
			ottables = Copyttables(ttables,state_count,Period);
			otcounts = Copytcounts(tcounts,state_count,Period);
		}
		added = 0;
		j = 0;
		for(i=0; i < DoIncrement; i++) {
			if(Use3D == 0) {
				incr = IncrementEpoch(0,tcounts,ttables,state_count,otcounts,ottables);
			} else {
				incr = IncrementEpoch(busy_epochs[j],tcounts,ttables,state_count,otcounts,ottables);
			}
			if(incr < 0) {
				printf("increment failed\n");
				exit(1);
			}
			added += incr;
			if(added > DoIncrement) {
				break;
			}
			if(Use3D == 1) {
				j = (j + 1) % BusyCount;
			}
		}
		free(busy_epochs);
		/*
		 * replace the original ttables and tcounts with updated
		 * verion
		 */
		if(Use3D == 1) {
			for(e=0; e < Period; e++) {
				free(ttables[e]);
				free(tcounts[e]);
				ttables[e] = ottables[e];
				tcounts[e] = otcounts[e];
			}
		} else {
			free(ttables[0]);
			free(tcounts[0]);
			ttables[0] = ottables[0];
			tcounts[0] = otcounts[0];
		}
		free(ottables);
		free(otcounts);
		if(Verbose == 1) {
			exit(0);
		}
	}
				

	if(Sample_count == 0) {
		Sample_count = SizeOf(input_data);
	}

	/*
	 * generate a synthetic series using the conditional transition
	 * probabilities
	 */


	/*
	 * start with zero state
 	 */
	theperiod = 0;
	src=0-min;
	memset(history_period,0-min,Period*sizeof(double));
	now = start_ts;
	for(i=0; i < Sample_count; i++) {
		r = drand48();
		/*
		 * find the dest state
		 */
		dest = -1;
		psrc = history_period[0];
		p = 0;
		counts = tcounts[theperiod];
		transitions = ttables[theperiod];
		for(j=0; j < state_count; j++) {
			if(counts[(src*state_count)+psrc] > 0.0) {
				p = (transitions[(src*state_count*state_count)+(psrc*state_count)+j] / 
					counts[(src*state_count)+psrc]) + p;
			}
			if(r <= p) {
				dest = j;
				break;
			}
		}
		/*
		 * we could generate a state we have never actually seen
		 *
		 * choose zero state in this case
		 */
	
		if(dest == -1) {
			if(counts[(src*state_count)+psrc] == 0.0) {
				dest = 0-(int)min;
			} else {
				fprintf(stderr,
				"error: i: %d, src: %d, psrc: %d dest -1\n",i,src,psrc);
				fprintf(stderr,"counts: %f\n",counts[(src*state_count)+psrc]);
				exit(1);
			}
		}
		printf("%10.0f %d\n",now,dest+(int)min);
		for(j=0; j < (Period-1); j++) {
			history_period[j] = history_period[j+1];
		}
		history_period[Period-1] = dest;
		src=dest;
		now += Interval;
		if(Use3D == 1) {
			theperiod = (theperiod + 1) % Period;
		}
	}

	if(Use3D == 1) {
		for(i=0; i < Period; i++) {
			free(ttables[i]);
			free(tcounts[i]);
		}
	} else {
		free(ttables[0]);
		free(tcounts[0]);
	}

		
	FreeDataSet(input_data);
	RBDestroyD(states);
	free(values);
	free(ttables);
	free(tcounts);
	return(0);
}

