#include <stdlib.h>
#include <stdio.h>
#include "io.h"

#define LINESIZE 128

char *inputfilename = DEFAULTINPUT;

inline static void strtodata(char *str, delays_mt *data, long i);



void parse(int argc, char *argv[])
{
	if(argc>1)
		inputfilename = argv[1];
	return;
}





delays_mt getinput(void)
{
	FILE *file = fopen(inputfilename, "r");
	if(file == NULL)
	{
		fprintf(stderr, "Cannot open input file\n");
		exit(EXIT_FAILURE);
	}

	delays_mt retval;

	char *str = NULL;
	str = malloc(LINESIZE);
	if(str == NULL)
	{
		fprintf(stderr, "Cannot malloc\n");
		exit(EXIT_FAILURE);
	}

	/*find number of lines (number of measurements)*/
	long long nblines = 0;
	while(fgets(str, LINESIZE, file) != NULL)
	{
		nblines++;
	}
	rewind(file);
	retval.length = nblines;

	retval.time = malloc(nblines*sizeof(time_mt));
	retval.delay = malloc(nblines*sizeof(delay_mt));

	
	for(long i = 0; fgets(str, LINESIZE, file); i++)
	{
			strtodata(str, &retval, i);
	}
	return retval;
}

/*The implementation changes if the delay_mt or time_mt change*/
inline static void strtodata(char *str, delays_mt *data, long i)
{
	char *fin;
	data->time[i] = strtod(str, &fin);
	data->delay[i] = strtod(fin, &str);
}









void plot1D(float_mt *data, long length)
{
	char tempname[L_tmpnam];
	tmpnam(tempname);
	FILE *file = fopen(tempname, "w");
	for(int i = 0; i<length; i++)
	{
		fprintf(file, "%f\n", data[i]);
	}
	fclose(file);
	
	int L = 100;
	char str[L];
	snprintf(str, L, "./plot.py %s\n", tempname);
	system(str);
	remove(tempname);
}





