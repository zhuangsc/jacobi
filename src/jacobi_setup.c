#include "jacobi_setup.h"

int main(int argc, char** argv){

	if(argc != 4){
		printf("Usage %s [filename] [block size] [csc(0)/csr(1)]\n", argv[0]);
		exit(1);
	}
	char* matrix_file = argv[1];
	int bs = atoi(argv[2]);
	int format = atoi(argv[3]);

	hbmat_t *L_csr, *L_csc, *L, *A, *A_csr, *A0;
	A = (hbmat_t*) malloc(sizeof(hbmat_t));

	unsigned long elapsed = 0;
	readHB_newmat_double(matrix_file, &(A->m), &(A->n), &(A->elemc), &(A->vptr), &(A->vpos), (double **)&(A->vval));
	struct timeval start, stop;
	gettimeofday(&start, NULL);

	one2zero(A);
	hb_print_CSC2("A000.dat", A);

	L = hb2hbh(A,bs);

	A0 = hbh2hb(L);
	hb_print_CSC2("L000.dat", A0);
	hbh_free(L);
	hb_free(A0);
	hb_free(A);
	return 0;
}

void print_matrix(const hbmat_t* matrix_info, int h, char* name){
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
