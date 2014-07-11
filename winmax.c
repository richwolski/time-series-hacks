#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

#include "simple_input.h"

char *Usage = "winmax -f filename -W timewindow\n";

#define ARGS "f:W:"

int Window;
char Fname[255];

	
#ifdef TEST
int main(int argc, char *argv[])
{
	int c;
	int err;
	void *input_data;
	void *xform_data;
	int fields;
	double *values;
	int i;
	double last_ts;
	double ts;
	double max_value;
	double value;

	memset(Fname,0,sizeof(Fname));
	while((c = getopt(argc,argv,ARGS)) != EOF)
	{
		switch(c)
		{
			case 'f':
				strncpy(Fname,optarg,sizeof(Fname));
				break;
			case 'W':
				Window = atoi(optarg);
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

	fields = GetFieldCount(Fname);

	if(fields <= 0)
	{
		fprintf(stderr,"bad field count in %s\n",Fname);
		exit(1);
	}

	values = (double *)malloc(fields*sizeof(double));
	if(values == NULL)
	{
		fprintf(stderr,"no space for values\n");
		exit(1);
	}

	err = InitDataSet(&input_data,fields);
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

	if(fields <= 1) {
		fprintf(stderr,
	"file must have at least two fields and first col is time stamps\n");
		exit(1);
	}

	last_ts = -1;
	max_value = 0;
	while(ReadData(input_data,fields,values) != 0)
	{
		/*
		 * assume ts is in fields 0
		 */
		ts = values[0];
		/*
		 * assume value of interest in is last field
		 */
		value = values[fields-1];
		if(last_ts == -1) {
			last_ts = ts;
			max_value = value;
		} else {
			if ((ts - Window) < last_ts) {
				if(value > max_value) {
					 max_value = value;
				}
			} else {
				printf("%10.0f %f\n",last_ts,max_value);
				last_ts = ts;
				max_value = value;
			}
		}
	}

	FreeDataSet(input_data);
	free(values);
	return(0);
}

#endif
