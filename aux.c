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


static void myquicksort(delays_mt *A, long p, long r)
{
	if(p>=r)
		return;
	long q = partition(A, p, r);

	myquicksort(A, p, q-1);
	myquicksort(A, q+1,r);
}


static void insertion(delays_mt *A, long L)
{
	time_mt key;
	delay_mt d;
	long i;
	for(long j=1; j<L; j++)
	{
		key = A->time[j];
		d = A->delay[j];
		i = j-1;
		while(i>=0 && A->time[i]>key)
		{
			A->time[i+1] = A->time[i];
			A->delay[i+1] = A->delay[i];
			i--;
		}
		A->time[i+1] = key;
		A->delay[i+1] = d;
	}
}


void mysort(delays_mt *data)
{
	//myquicksort(data, 0, data->length-1);	
	insertion(data, data->length);
	return;
}


float_mt prod(float_mt *data, int n)
{
	float_mt res = 1.0;
	if(n<1)
		return 1.0;
	for(int i =0; i<n; i++)
	{
		res *= data[i];
	}
	return res;
}


inline float_mt sum(float_mt *data, int n)
{
	float_mt res = 0.0;
	if(n<1)
		return 0.0;
	for(int i = 0; i<n; i++)
	{
		res += data[i];
	}
	return res;
}













