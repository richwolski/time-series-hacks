#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "simple_input.h"
#include "red_black.h"

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
\t-c sample_count\n\
\t-i time interval (for synthetic trace)\n\
\t-P cyclic_period (measured in lags)\n\
\t-3 <use 3D method>\n\
\t-V <verbose mode>\n";

#define ARGS "f:Vc:i:P:3"

char Fname[255];
int Verbose;
int Sample_count;
double Interval;
int Period;
int Use3D;


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
	int tp;
	int theperiod;
	int lastperiod;

	memset(Fname,0,sizeof(Fname));
	Interval = 1;
	Period = 0;
	Use3D = 0;
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
			case 'i':
				Interval = atof(optarg);
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

	states = RBTreeInit();
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
		rb = RBFind(states,values[field_count-1]);
		if(rb != NULL) {
			count = rb->value.d + 1.0; 
			hv.d = count;
			RBInsert(states,values[field_count-1],hv);
		} else {
			if(values[field_count-1] > max) {
				max = values[field_count-1];
			}
			if(values[field_count-1] < min) {
				min = values[field_count-1];
			}
	
			hv.d = 1.0;
			RBInsert(states,values[field_count-1],hv);
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
	RBDeleteTree(states);
	free(values);
	free(ttables);
	free(tcounts);
	return(0);
}

