#include <math.h>
#include "jacobi_setup.h"
#include "jacobi_main.h"

void *Ahb;
void *Ahbh;
double *v_b;
double *v_x;
double *v_x0;
int format;
int bs;
int dim;
int dim0;
long iter;

int main(int argc, char* argv[]){

	if(argc != 4){
		printf("Usage %s [HB file(A)] [block size] [csc(0)/csr(1)]\n", argv[0]);
		exit(1);
	}

	if (jacobi_setup(argc, argv)){
		printf("Unable to proceed, exit\n");
		exit(1);
	}

	jacobi_main_csr(Ahbh, v_x0, v_b, bs);

	printf("Converge at iteration: %ld\n", iter);
	double x = 0;
	double x0 = 0;
//	for(int i = 0; i < dim0; ++i)
//		printf("%lf, ", v_x[i]);
//	puts("");
//	for(int i = 0; i < dim0; ++i)
//		printf("%lf, ", v_x0[i]);
//	puts("");

	for(int i = 0; i < dim0; ++i) 
		x += v_x[i] * v_x[i];
	x = sqrt(x);
	for(int i = 0; i < dim0; ++i) 
		x0 += v_x0[i] * v_x0[i];
	x0 = sqrt(x0);
	printf("2-norm x: %lf, 2-norm x0: %lf\n", x, x0);
	printf("2-norm : %lf\n", (x-x0)/x);

//	jacobi_shutdown();

	return 0;
}
