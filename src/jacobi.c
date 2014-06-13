
#include "jacobi_setup.h"

void *Ahb;
void *Ahbh;
double *v_b;
double *v_x;
double *v_x0;
int format;
int bs;

int main(int argc, char* argv[]){

	if(argc != 4){
		printf("Usage %s [HB file(A)] [block size] [csc(0)/csr(1)]\n", argv[0]);
		exit(1);
	}

	if (jacobi_setup(argc, argv)){
		printf("Unable to proceed, exit\n");
		exit(1);
	}

	jacobi_shutdown();

	return 0;
}
