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

char *Usage = "markov-series -f filename\n\
\t-c sample_count\n\
\t-i time interval (for synthetic trace)\n\
\t-V <verbose mode>\n";

#define ARGS "f:Vc:i:"

char Fname[255];
int Verbose;
int Sample_count;
double Interval;


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
	int src;
	int dest;
	double max;
	double min;
	double *counts;
	double now;
	double r;
	double p;
	double prev_p;
	double start_ts;

	memset(Fname,0,sizeof(Fname));
	Interval = 1;
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


	/*
	 * cover all possible states
	 */
	state_count = (int)(max - min) + 1;

	transitions = (double *)malloc(state_count*state_count*sizeof(double));
	if(transitions == NULL) {
		exit(1);
	}
	memset(transitions,0,state_count*state_count*sizeof(double));

	counts = (double *)malloc(state_count*sizeof(double));
	memset(counts,0,state_count*sizeof(double));

	Rewind(input_data);
	/*
	 * compute the transition counts and the count of each state
	 * as a src
	 */
	src = -10000000;
	while(ReadData(input_data,field_count,values) != 0) {
		if(src == -10000000) {
			src = (int)(values[field_count-1] - min);
			continue;
		}
		dest = (int)(values[field_count-1]-min);
		transitions[src*state_count+dest] += 1.0;
		counts[src] += 1.0;
		src = dest;
	}

	if(Verbose == 1) {
		count = SizeOf(input_data);
		printf("       ");
		for(i=0; i < state_count; i++) {
			printf("%6.6d ",i+(int)min);
		}
		printf("\n");
		for(i=0; i < state_count; i++) {
			printf("--------");
		}
		printf("\n");
		for(i=0; i < state_count; i++) {
			printf("%4.4d | ",i+(int)min);
			for(j=0; j < state_count; j++) {
				printf("%4.4f ",transitions[i*state_count+j]/counts[i]);
	//			printf("%4.4f ",transitions[i*state_count+j]);
	//			printf("%4.4f ",transitions[i*state_count+j]/(double)count);
			}
			printf("\n");
		}
	}

	if(Sample_count == 0) {
		Sample_count = SizeOf(input_data);
	}
	/*
	 * generate a synthetic series using the conditional transition
	 * probabilities
	 */
	src=0;
	now = start_ts;
	for(i=0; i < Sample_count; i++) {
		r = drand48();
		/*
		 * find the dest state
		 */
		dest = -1;
		p = 0;
		for(j=0; j < state_count; j++) {
			prev_p = p;
			p = (transitions[src*state_count+j] / counts[src]) + p;
			if(r <= p) {
				dest = j;
				break;
			}
/*
			if(j != 0) {
				prev_p = transitions[src*state_count+(j-1)] / counts[src];
			}
			if((j == 0) && (r <= p)) {
				dest = j;
				break;
			}
*/
/*
			} else if((j != 0) && (r > prev_p) && (r <= p)) {
				dest = j;
				break;
			} else if((j == (state_count-1)) && (r > p)) {
				dest = j;
				break;
			}
*/
		}
		// sanity check
		if(dest == -1) {
			fprintf(stderr,
			"error: i: %d, src: %d, dest -1\n",i,src);
			exit(1);
		}
		printf("%10.0f %d\n",now,dest+(int)min);
		src=dest;
		now += Interval;
	}

		
	FreeDataSet(input_data);
	RBDeleteTree(states);
	free(values);
	free(transitions);
	return(0);
}

