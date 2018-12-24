#include <iostream>
#include <sys/sysinfo.h>
#include <cstdio>
using namespace std;


void printMem()
{
	struct sysinfo myinfo;
	sysinfo(&myinfo); // there is a struct called sysinfo and a function called sysinfo as well
	double total_mem = myinfo.mem_unit * myinfo.totalram;
	total_mem /= (1024.*1024.);
	double free_mem = myinfo.mem_unit * myinfo.freeram;
	free_mem /= (1024.*1024.);

	printf("printMem node d: Memory: total: %.2f MB, avail: %.2f MB, used %.2f MB\n", total_mem, free_mem, total_mem-free_mem);
}


int main(){
	printMem();
	return 0;
}
