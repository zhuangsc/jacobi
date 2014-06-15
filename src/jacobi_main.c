
#include "jacobi_main.h"

extern long iter;
extern int dim;
extern int dim0;

int jacobi_main_csr(hbmat_t *Acsr, double *x, double *b, int bs){

	int M = Acsr->m; int N = Acsr->n;
	int *vdiag = Acsr->vdiag;
	int *vptr = Acsr->vptr; int *vpos = Acsr->vpos; hbmat_t **vval = Acsr->vval;
	int converge = 0;

	hbmat_t *Adiag, *Acurrent;
	hbmat_t **diagL = malloc(N * sizeof(hbmat_t*) );
	double *vtmp = malloc(bs * N * sizeof(int));
	int *work = malloc(bs * sizeof(int));
	for ( int j = 0; j < N; ++j ) {
		hbmat_t** Mtmp;
		Adiag = vval[vdiag[j]];
		Mtmp = ompss_csr_dchol_ll(bs, Adiag, work);
#pragma omp taskwait
		diagL[j] = hbh2hb(Mtmp);
//		print_matrix(diagL[j],0,"diagL");
		hbh_free(Mtmp);
	}
	free(work);
	puts("pass cholesky");

#pragma omp taskwait

	double *x0 = malloc(bs * sizeof(double));
	double *x1 = malloc(bs * sizeof(double));
	iter = 0;
	while(1) {
		++iter;
		for(int k = 0; k < dim; ++k)
			vtmp[k] = 0.0;

		for ( int I = 0; I < N; ++I ) {
			for ( int J = vptr[I]; J < vptr[I+1]; ++J ) {
				if ( vpos[J] == I )
					continue;
				Acurrent = vval[J];
				jacobi_dgemv_csr( Acurrent, &(x[vpos[J] * bs]), &(vtmp[I * bs]) );
//				puts("------------------------");
//				for ( int k = 0; k < dim; ++k)
//					printf("%lf ",vtmp[k]);
//				puts("");

			}
			jacobi_dsubvv( b, vtmp, I, bs );
//			puts("++++++++++++++++++++++++");
//			for ( int k = 0; k < dim; ++k)
//				printf("%lf ",vtmp[k]);
//			puts("");

			/*
			 *	A[i,i] * x[i] = vtmp[i]
			 */
			Adiag = diagL[I];
			jacobi_dtrsm_csr( Adiag, &(vtmp[I * bs]), x0 );
			jacobi_dtrsmt_csr( Adiag, x0, x1 );
			for ( int k = 0; k < bs; ++k ) {
				vtmp[I * bs + k] = x1[k];
			}
		}
		for ( int k = 0; k < dim0; ++k) {
			if ( x[k] - vtmp[k] ) {
				converge = 0;
				break;
			}
			converge = 1;
		}

//		printf("Iter: %d\n", iter);
//		for ( int k = 0; k < dim; ++k)
//			printf("%lf ",vtmp[k]);
//		puts("");

		for ( int k = 0; k < dim; ++k) {
			x[k] = vtmp[k];
		}

		if (converge || iter == 1000) 
			break;
	}

	return 0;
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
