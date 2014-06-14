#include "jacobi_setup.h"

extern void *Ahb;
extern void *Ahbh;
extern double *v_b;
extern double *v_x;
extern double *v_x0;
extern int format;
extern int bs;

int jacobi_setup (int argc, char* argv[]) {

	char* matrix_file = argv[1];
	int bs = atoi(argv[2]);
	int format = atoi(argv[3]);

	Ahb = (hbmat_t*) malloc(sizeof(hbmat_t));
	hbmat_t *A = Ahb;

	readHB_newmat_double(matrix_file, &(A->m), &(A->n), &(A->elemc), &(A->vptr), &(A->vpos), (double **)&(A->vval));

	one2zero(A);
	hb_print_CSC2("A000.dat", Ahb);

	int *work = malloc(A->m * sizeof(int));
	hbmat_t *ttt;
	ttt = ompss_csr_dchol_ll(bs, Ahb, work);
	print_matrix(ttt,1,"ttt");

	int dim = A->m;
	v_b = malloc(dim * sizeof(double));
	v_x = malloc(dim * sizeof(double));
	v_x0 = calloc(dim, sizeof(double));
	result_gen(v_x, dim);
	mkl_cspblas_dcsrgemv("N", &dim, A->vval, A->vptr, A->vpos, v_x, v_b);

#if 0
	for (int i = 0; i < dim; ++i)
		printf("%lf \n", v_x[i]);
	puts("");
	for (int i = 0; i < dim; ++i)
		printf("%lf \n", v_b[i]);
	puts("");
	Ahbh = hb2hbh(A, bs, format);
	hbmat_t *A1 = Ahbh;
	for (int i = 0; i < A1->n; ++i)
		printf("%d\n", (A1->vdiag)[i]);
#endif

	Ahbh = hb2hbh(A, bs, format);
//	hbmat_t *A0;
//	A0 = hbh2hb(Ahbh);
//	hb_print_CSC2("L000.dat", A0);
//	hb_free(A0);

	return 0;
}

void jacobi_shutdown () {
	hb_free(Ahb);
	hbh_free(Ahbh);
}

void result_gen (double *vector, int length) {
	srand48(time(0));
	for(int i = 0; i < length; ++i){
		vector[i] = drand48() * 50;
	}
}
