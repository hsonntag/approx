#include <pthread.h>
#include <stdio.h>
#include <cspl/cspl.h>
#define NUM_THREADS  256
#define SIGNAL_LENGTH 65536
void *xcorr(void *threadid)
{
	long tid;
	tid = (long)threadid;
	int i;
	double c_xy[SIGNAL_LENGTH];
	double x[SIGNAL_LENGTH];
	double y[SIGNAL_LENGTH];
	for (i = 0; i < SIGNAL_LENGTH; i++)
	{
		if (i < 100)
		{
			x[i] = (double) i;
			y[i] = (double) i;
		}
		else if (i < 200)
		{
			x[i] = (double) (100 - i);
			y[i] = (double) (100 - i);
		}
		else if (i < 300)
		{
			x[i] = 0.0;
			y[i] = 0.0;
		}
		else if (i < 400)
		{
			x[i] = (double) i;
			y[i] = 0.0;
		}
		else
		{
			x[i] = (double) (100 - i);
			y[i] = 0.0;
		}
	}
	cspl_radix2_xcorr (c_xy, x, y, SIGNAL_LENGTH); 
		pthread_exit(NULL);
}

int main (int argc, char *argv[])
{
	pthread_t threads[NUM_THREADS];
	int rc;
	long t;
	for(t=0; t<NUM_THREADS; t++){
		printf("In main: creating thread %ld\n", t);
		rc = pthread_create(&threads[t], NULL, xcorr, (void *)t);
		if (rc){
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
	}
	pthread_exit(NULL);
}
