#include <sys/resource.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int sys_tracker_(){
	int i=0;
	long mb=0;
	struct rusage r_usage;
	getrusage(RUSAGE_SELF,&r_usage);
	mb = r_usage.ru_maxrss;
	printf("Maximum Memory usage = %ld mb\n",mb/1000);
	return 0;

}