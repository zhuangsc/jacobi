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
#endif

	one2zero(A);
	hb_print_CSC2("A000.dat", Ahb);

	Ahbh = hb2hbh(A, bs, format);

	hbmat_t *A0;
	A0 = hbh2hb(Ahbh);
	hb_print_CSC2("L000.dat", A0);
	hb_free(A0);

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

void print_matrix (const hbmat_t* matrix_info, int h, char* name) {
	double* value = (double*) matrix_info->vval;
	hbmat_t** address = matrix_info->vval;
	printf("------------------------------------\n");
	printf("\t Displaying : %s\n", name);
	printf("Rows = %d,Columns = %d,Non-zeros = %d\n",matrix_info->m,matrix_info->n,matrix_info->elemc);
	printf("Virtual Address = %p\n", matrix_info);
	for (int i = 0; i <= matrix_info->n; i++)
		printf("%d ", matrix_info->vptr[i]);
	printf("\n");
	for (int j = 0; j < matrix_info->elemc; j++)
		printf("%d ", matrix_info->vpos[j]);
	printf("\n");
	for (int j = 0; j < matrix_info->elemc; j++){
		if(h == 0)
			printf("%f ", value[j]);
		else{
			printf("vval[%d]: %p\n", j, address[j]);
			print_matrix(address[j],0,"");
		}
	}
	printf("\n");
}
