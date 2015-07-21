#include <sys/resource.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int sys_tracker_(){
	int i=0;
	long mb=0;
	float s_sec = 0.0;
	float u_sec = 0.0;
	struct rusage r_usage;

	getrusage(RUSAGE_SELF,&r_usage);
	mb = r_usage.ru_maxrss/1000;
	s_sec = (float)r_usage.ru_stime.tv_sec + 1e-6*(float)r_usage.ru_stime.tv_usec;
	u_sec = (float)r_usage.ru_utime.tv_sec + 1e-6*(float)r_usage.ru_utime.tv_usec;
	printf("	Maximum Memory usage = %ld mb\n",mb);
	printf("	System time elapsed  = %5.2f s\n",s_sec);
	printf("	User time elapsed    = %5.2f s\n",u_sec);
	return 0;

}