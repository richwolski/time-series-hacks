#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mio.h"
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

char *Usage = "inar-seasonal-series -f filename\n\
\t-c count\n\
\t-P cyclic_period (measured in lags)\n\
\t-V <verbose mode>\n";

#define ARGS "f:P:V"

char Fname[255];
int Verbose;
int Sample_count;
int Period;

/*
 * taken from A Poisson INAR(1) process with a seasonal structure
 */
double EstimatePhi(double *data, unsigned long size, int f, int period)
{
	int i;
	int j;
	double y_t;
	double y_t_minus_s;
	double sum1 = 0.0;
	double sum2 = 0.0;
	double sum3 = 0.0;
	double sum4 = 0.0;
	double phi;

	

	for(i=period; i < size; i++) {
		y_t = data[i*f+(f-1)];
		y_t_minus_s = data[(i-period)*f+(f-1)];

		sum2 += y_t;
		sum3 += y_t_minus_s;
		sum4 += (y_t_minus_s * y_t_minus_s);
	}

	for(i=period; i < size; i++) {
		y_t = data[i*f+(f-1)];
		y_t_minus_s = data[(i-period)*f+(f-1)];

		sum1 += (y_t * y_t_minus_s);
	}


	phi = (((size - period)*sum1) - (sum2*sum3)) /
		(((size - period)*sum4) - (sum3*sum3));

	return(phi);
}
	


int main(int argc, char *argv[])
{
	int c;
	int err;
	int i;
	int j;
	int k;
	MIO *data_mio;
	MIO *raw_mio;
	double *data;
	int data_fields;
	unsigned long size;
	unsigned long recs;
	double phi;

	memset(Fname,0,sizeof(Fname));
	Period = 0;
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
			case 'P':
				Period = atoi(optarg);
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
	recs = MIOTextRecords(raw_mio);

	MIOClose(raw_mio);

	data = MIOAddr(data_mio);

	phi = EstimatePhi(data, recs, data_fields, Period);

	if(Verbose == 1) {
		printf("CLS phi: %f\n",phi);
	}

	MIOClose(data_mio);

	return(0);
}

