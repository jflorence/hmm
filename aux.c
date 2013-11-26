#include <stdlib.h>
#include "aux.h"
#include <stdio.h>
static inline long partition(delays_mt *A, long p, long r)
{
	time_mt x = A->time[r];
	int i = p-1;
	int k;
	time_mt tempt;
	delay_mt tempd;
	for(int j=p; j<r;j++)
	{
		//printf("%d\n", j);
		if(A->time[j]<=x)
		{
			i++;

			tempt = A->time[i];
			tempd = A->delay[i];
			A->time[i] = A->time[j];
			A->delay[i] = A->delay[j];
			A->time[j] = tempt;
			A->delay[j] = tempd;
		}
	}
	k = i+1;
	tempt = A->time[k];
	tempd = A->delay[k];
	A->time[k] = A->time[r];
	A->delay[k] = A->delay[r];
	A->time[r] = tempt;
	A->delay[r] = tempd;
	return k;
}


void myquicksort(delays_mt *A, long p, long r)
{
	if(p>=r)
		return;
	long q = partition(A, p, r);

	myquicksort(A, p, q-1);
	myquicksort(A, q+1,r);
}



void mysort(delays_mt *data)
{
	myquicksort(data, 0, data->length-1);	
	return;
}










