#include <stdlib.h>
#include <stdio.h>
#include "io.h"
#include <stdarg.h>
#define LINESIZE 128

char *inputfilename = DEFAULTINPUT;

inline static void strtodata(char *str, delays_mt *data, long long i);



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
inline static void strtodata(char *str, delays_mt *data, long long i)
{
	char *fin;
	data->time[i] = strtod(str, &fin);
	data->delay[i] = strtod(fin, &str);
}









void plot1D(float_mt *data, long long length)
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




#define NUMSPACE 2
void dispstr(int level, char *str, ...)
{
#ifdef VERBOSE
	va_list params;
	va_start(params, str);
	for(int i = 0; i<level*NUMSPACE; i++)
		printf(" ");
	vprintf(str, params);
	va_end(params);
#endif
}




const char *space = "    ";
void write_results(struct params *p)
{
	FILE *file = fopen(OUTPUTPARAMS, "w");
	fprintf(file, "Those are the parameters of a HSMM trained with the data of the file ...\n\n");
	fprintf(file, "A:     %s", space);
	for(int i=0; i<p->N; i++)
	{
		if(i>0)
			fprintf(file, "       %s", space);
		for(int j=0; j<p->N; j++)
		{
			fprintf(file, "%2.4f%s", p->A[i*p->N+j], space);
		}
		fprintf(file, "\n");
	}
	fprintf(file, "shape: %s", space);
	for(int i = 0; i<p->N; i++)
		fprintf(file,"%2.4f%s", p->shape[i], space);
	fprintf(file,"\n");
	fprintf(file, "scale: %s", space);
	for(int i = 0; i<p->N; i++)
		fprintf(file,"%2.4f%s", p->scale[i], space);
	fprintf(file,"\n");
	fprintf(file, "dshape:%s", space);
	for(int i = 0; i<p->N; i++)
		fprintf(file,"%2.4f%s", p->dshape[i], space);
	fprintf(file,"\n");
	fprintf(file, "dscale:%s", space);
	for(int i = 0; i<p->N; i++)
		fprintf(file,"%2.4f%s", p->dscale[i], space);
	fprintf(file,"\n");

	fclose(file);
}





void write_alpha_beta(float_mt *alpha, float_mt *beta, int N, long long L)
{
	FILE *file = fopen("alpha_beta.txt","w");

	for(long t=0; t<L; t++)
	{
		for(int i=0; i<N; i++)
		{
			fprintf(file, "%2.4f%s",alpha[t*N+i],space);
		}
		fprintf(file, "%s%s", space, space);
		for(int i=0; i<N; i++)
		{
			fprintf(file, "%2.4f%s",beta[t*N+i],space);
		}
		fprintf(file, "\n");

	}

	fclose(file);
}














