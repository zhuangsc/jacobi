#include "jacobi_setup.h"

extern void *Ahb;
extern void *Ahbh;
extern double *v_b;
extern double *v_x;
extern double *v_x0;
extern int format;
extern int bs;
extern int dim;
extern int dim0;

int jacobi_setup (int argc, char* argv[]) {

	char* matrix_file = argv[1];
	bs = atoi(argv[2]);
	format = atoi(argv[3]);

	Ahb = (hbmat_t*) malloc(sizeof(hbmat_t));
	hbmat_t *A = Ahb;

	readHB_newmat_double(matrix_file, &(A->m), &(A->n), &(A->elemc), &(A->vptr), &(A->vpos), (double **)&(A->vval));

	one2zero(A);
//	int *work = malloc(bs * sizeof(int));
//	hbmat_t *tmp = hbh2hb(ompss_csr_dchol_ll(A->n, A, work));
//	hb_print_CSC2("L111.dat", tmp);

	dim0 = A->n;
	Ahbh = hb2hbh(A, bs, format);
	hbmat_t *A0 = Ahbh;

	dim = A0->n * bs;
	v_b = malloc(dim * sizeof(double));
	v_x = malloc(dim * sizeof(double));
	v_x0 = calloc(dim, sizeof(double));
	result_gen(v_x, dim);
	mkl_cspblas_dcsrgemv("N", &dim0, A->vval, A->vptr, A->vpos, v_x, v_b);

#if 0
	for (int i = 0; i < dim0; ++i)
		printf("%lf ", v_x[i]);
	puts("\n");
	for (int i = 0; i < dim0; ++i)
		printf("%lf ", v_b[i]);
	puts("\n\n");
#endif


	return 0;
}

void jacobi_shutdown () {
	hb_free(Ahb);
	hbh_free(Ahbh);
	free(v_b); free(v_x); free(v_x0);
}

void result_gen (double *vector, int length) {
	srand48(time(0));
	for(int i = 0; i < length; ++i){
		vector[i] = drand48() * 9;
	}
}
