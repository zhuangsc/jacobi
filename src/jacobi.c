
#include "jacobi_setup.h"

void *Ahb;
void *Ahbh;
double *v_b;
double *v_x;
int format;
int bs;

int main(int argc, char* argv[]){

	if(argc != 5){
		printf("Usage %s [HB file(A)] [block size] [csc(0)/csr(1)] [Vector(b)]\n", argv[0]);
		exit(1);
	}

	jacobi_setup(argc, argv);

	jacobi_shutdown();

	return 0;
}
