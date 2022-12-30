#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

#include "simple_input.h"

/*
 * for WTB interarrival times where cameras take sequences
 * when ever they are triggered
 *
 * interval is max interval between images in a run
 */

char *Usage = "interarrival-runs -f filename\n\
\t -i interval\n";

#define ARGS "f:i:"

double Interval;
char Fname[255];


int main(int argc, char *argv[])
{
	int c;
	int err;
	void *input_data;
	int field_count;
	double *values;
	double last_arrival;
	double run_start;
	double interarrival;
	int run_length;
	double run_duration;

	memset(Fname,0,sizeof(Fname));
	while((c = getopt(argc,argv,ARGS)) != EOF)
	{
		switch(c)
		{
			case 'f':
				strncpy(Fname,optarg,sizeof(Fname));
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

	last_arrival = -1;
	while(ReadData(input_data,field_count,values) != 0) {
		/*
		 * assume timestamp is in field 0
		 * and that the data is in time sorted order
		 */
		if(last_arrival == -1) {
			last_arrival = values[0];
			run_duration = 0.0;
			run_length = 0;
			run_start = last_arrival;
			continue;
		}
		interarrival = values[0] - last_arrival;
		if(interarrival <= Interval) {
			run_duration += interarrival;
			run_length++;
			last_arrival = values[0];
			continue;
		}

		/*
		 * run has ended here
		 */
		printf("%10f %f %d %f\n",
			run_start,
			run_duration,
			run_length,
			run_duration / (float)run_length);
		last_arrival = run_start = values[0];
		run_duration = 0.0;
		run_length = 0;
	}

		
	FreeDataSet(input_data);
	free(values);
	return(0);
}

