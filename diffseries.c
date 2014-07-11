#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

#include "simple_input.h"

char *Usage = "diffseries -f filename\n\
\t -d difference\n";

#define ARGS "f:d:"

int Diff;
char Fname[255];

void *DiffSeries(void *series, int fields)
{
	int err;
	void *diff;
	double *data;
	double *last_data;
	double old_datum;
	int i;

	data = (double *)malloc(fields * sizeof(double));
	if(data == NULL)
	{
		return(NULL);
	}

	last_data = (double *)malloc(fields * sizeof(double));
	if(last_data == NULL)
	{
		free(data);
		return(NULL);
	}
	
	err = InitDataSet(&diff,fields);
	if(err <= 0)
	{
		free(data);
		free(last_data);
		return(NULL);
	}

	Rewind(series);

	err = ReadData(series,fields,last_data);
	if(err <= 0)
	{
		free(data);
		free(last_data);
		FreeDataSet(diff);
		return(NULL);
	}

	while(ReadData(series,fields,data))
	{
		old_datum = data[fields-1];
		data[fields-1] = data[fields-1] - last_data[fields-1];
		WriteData(diff,fields,data);
		data[fields-1] = old_datum;
		for(i=0; i < fields; i++)
		{
			last_data[i] = data[i];
		}
	}

	free(data);
	free(last_data);

	Rewind(diff);
	return(diff);
}
	

#ifdef TEST
int main(int argc, char *argv[])
{
	int c;
	int err;
	void *input_data;
	void *xform_data;
	int field_count;
	double *values;
	int i;
	int diff;

	memset(Fname,0,sizeof(Fname));
	while((c = getopt(argc,argv,ARGS)) != EOF)
	{
		switch(c)
		{
			case 'f':
				strncpy(Fname,optarg,sizeof(Fname));
				break;
			case 'd':
				Diff = atoi(optarg);
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

	for(i=0; i < Diff; i++)
	{
		xform_data = DiffSeries(input_data,field_count);
		if(xform_data == NULL)
		{
			fprintf(stderr,"couldn't comput %d Diff\n",i);
			free(values);
			exit(1);
		}
		FreeDataSet(input_data);
		input_data = xform_data;
	}

	Rewind(input_data);
	while(ReadData(input_data,field_count,values) != 0)
	{
		for(i=0; i < field_count; i++)
		{
			fprintf(stdout,"%f ", values[i]);
		}
		fprintf(stdout,"\n");
		fflush(stdout);
	}

	FreeDataSet(input_data);
	free(values);
	return(0);
}

#endif
