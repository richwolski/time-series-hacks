#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

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

char *Usage = "markov-series -f filename\n";

#define ARGS "f:"

char Fname[255];


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
	double *counts;

	memset(Fname,0,sizeof(Fname));
	while((c = getopt(argc,argv,ARGS)) != EOF)
	{
		switch(c)
		{
			case 'f':
				strncpy(Fname,optarg,sizeof(Fname));
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

	max = -1;
	state_count = 0;
	while(ReadData(input_data,field_count,values) != 0) {
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
			hv.d = 1.0;
			RBInsert(states,values[field_count-1],hv);
		}
	}


	/*
	 * cover all possible states
	 */
	state_count = (int)max + 1;

	transitions = (double *)malloc(state_count*state_count*sizeof(double));
	if(transitions == NULL) {
		exit(1);
	}
	memset(transitions,0,state_count*state_count*sizeof(double));

	counts = (double *)malloc(state_count*sizeof(double));
	memset(counts,0,state_count*sizeof(double));

	Rewind(input_data);
	src = -1;
	while(ReadData(input_data,field_count,values) != 0) {
		if(src == -1) {
			src = (int)values[field_count-1];
			continue;
		}
		dest = (int)values[field_count-1];
		transitions[src*state_count+dest] += 1.0;
		counts[src] += 1.0;
		src = dest;
	}

	count = SizeOf(input_data);

	printf("       ");
	for(i=0; i < state_count; i++) {
		printf("%6.6d ",i);
	}
	printf("\n");
	for(i=0; i < state_count; i++) {
		printf("--------");
	}
	printf("\n");
	for(i=0; i < state_count; i++) {
		printf("%4.4d | ",i);
		for(j=0; j < state_count; j++) {
//			printf("%4.4f ",transitions[i*state_count+j]/counts[i]);
//			printf("%4.4f ",transitions[i*state_count+j]);
			printf("%4.4f ",transitions[i*state_count+j]/(double)count);
		}
		printf("\n");
	}

		
	FreeDataSet(input_data);
	RBDeleteTree(states);
	free(values);
	free(transitions);
	return(0);
}

