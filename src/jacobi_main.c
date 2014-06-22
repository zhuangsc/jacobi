
#include "jacobi_main.h"
#define MAX_ITER 1000

extern long iter;
extern int dim;
extern int dim0;

int jacobi_main_csr(hbmat_t *Acsr, double *x, double *b, int bs){

	int M = Acsr->m; int N = Acsr->n;
	int *vdiag = Acsr->vdiag;
	int *vptr = Acsr->vptr; int *vpos = Acsr->vpos; hbmat_t **vval = Acsr->vval;
	int converge = 0;

	hbmat_t *Adiag, *Adiag0, *Acurrent;
	hbmat_t **diagL = malloc(N * sizeof(hbmat_t*) );
	double *vtmp = malloc(dim * sizeof(double));

	double **x0 = malloc(N * sizeof(double*));
	for(int k = 0; k < N; ++k){
		x0[k] = malloc(bs * sizeof(double));
	}
		
	double *x1 = malloc(dim * sizeof(double));
	
////////////////////1st Iteration//////////////////

	hbmat_t **A0 = malloc(dim *sizeof(hbmat_t*));
	int **etree0 = malloc(dim * sizeof(int*));
	for(int k = 0; k < dim; ++k){
		vtmp[k] = 0.0;
	}
	int **work = malloc(N * sizeof(int*));
	for(int k = 0; k < N; ++k)
		work[k] = malloc(bs*sizeof(int));
	for ( int I = 0; I < N; ++I ) {
		etree0[I] = etree(vval[vdiag[I]]);
		A0[I] = hb2hbh_sym_etree_csr_p(vval[vdiag[I]], bs, etree0[I]);
		diagL[I] = ((hbmat_t**)A0[I]->vval)[0];
		potrf_sparse_csr(diagL[I]);
		for ( int J = vptr[I]; J < vptr[I+1]; ++J ) {
			if ( vpos[J] != I ){
				Acurrent = vval[J];
				jacobi_dgemv_csr( Acurrent, &(x[vpos[J] * bs]), &(vtmp[I * bs]) );
			}
		}
		jacobi_dsubvv( b, vtmp, I, bs );

  	/*
  	 *	A[i,i] * x[i] = vtmp[i]
  	 */
		jacobi_dtrsm_csr( diagL[I], &(vtmp[I * bs]), x0[I]);
		jacobi_dtrsmt_csr( diagL[I], x0[I], &(x1[I*bs]) );
	}
#pragma omp taskwait

//	for ( int k = 0; k < dim; ++k) {
//		vtmp[k] = x1[k];
//	}

	for ( int k = 0; k < dim; ++k) {
		x[k] = x1[k];
	}

///////////////////////////////////////////
	iter = 0;
	while(!converge && iter < MAX_ITER) {
		++iter;
//		printf("Iter: %d\n", iter);
		for(int k = 0; k < dim; ++k){
			vtmp[k] = 0.0;
		}
		for ( int I = 0; I < N; ++I ) {
			for ( int J = vptr[I]; J < vptr[I+1]; ++J ) {
				if ( vpos[J] != I ){
					Acurrent = vval[J];
					jacobi_dgemv_csr( Acurrent, &(x[vpos[J] * bs]), &(vtmp[I * bs]) );
				}
			}
			jacobi_dsubvv( b, vtmp, I, bs );

			/*
			 *	A[i,i] * x[i] = vtmp[i]
			 */
			jacobi_dtrsm_csr( diagL[I], &(vtmp[I * bs]), x0[I] );
			jacobi_dtrsmt_csr( diagL[I], x0[I], &(x1[I*bs]) );
		}

#pragma omp taskwait

		for ( int k = 0; k < dim; ++k) {
			vtmp[k] = x1[k];
		}

		converge = 1;
		for ( int k = 0; k < dim0 && converge; ++k) {
			if ( x[k] - vtmp[k] ) 
				converge = 0;
		}

		for ( int k = 0; k < dim; ++k) {
			x[k] = vtmp[k];
		}
		converge = 0;
	}

	puts("through");
	free(x1);
	for(int k = 0; k < N; ++k)
		free(x0[k]);
	free(x0);

	//free(work); free(vtmp);
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
