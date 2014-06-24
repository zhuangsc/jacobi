
#include "jacobi_main.h"
#define MAX_ITER 5

extern long iter;
extern int dim;
extern int dim0;


int jacobi_main_csr(hbmat_t *Acsr, double *x, double *b, int bs){

	int M = Acsr->m; int N = Acsr->n;
	int *vdiag = Acsr->vdiag;
	int *vptr = Acsr->vptr; int *vpos = Acsr->vpos; hbmat_t **vval = Acsr->vval;
	int converge = 0;

	hbmat_t *Acurrent;
	hbmat_t **diagL = malloc(N * sizeof(hbmat_t*) );

	double *vtmp = calloc(dim , sizeof(double));
	double *ltmp = calloc( N * N * bs , sizeof(double));

	double *x0 = malloc(N * N * bs * sizeof(double));
	double *x1 = malloc(N * N * bs * sizeof(double));
	
	double *b0 = malloc(N * bs * sizeof(double));
	double *b1 = malloc(N * bs * sizeof(double));
	
////////////////////1st Iteration//////////////////

	hbmat_t **A0 = malloc(dim *sizeof(hbmat_t*));
	int **etree0 = malloc(dim * sizeof(int*));

	for ( int I = 0; I < N; ++I ) {
		etree0[I] = etree(vval[vdiag[I]]);
		A0[I] = hb2hbh_sym_etree_csr_p(vval[vdiag[I]], bs, etree0[I]);
		diagL[I] = ((hbmat_t**)A0[I]->vval)[0];
		potrf_sparse_csr(diagL[I]);
		int lc = 0;
		for ( int J = vptr[I]; J < vptr[I+1]; ++J ) {
			if ( vpos[J] != I ){
				Acurrent = vval[J];
				jacobi_dgemv_csr( Acurrent, &(x[vpos[J] * bs]), &(ltmp[I * N * bs + lc * bs]));
				lc++;
			}
		}

  	/*
  	 *	A[i,i] * x[i] = vtmp[i]
  	 */
		jacobi_dtrsm_csr( diagL[I], &(b[I*bs]), &(b0[I*bs]));
		jacobi_dtrsmt_csr( diagL[I], &(b0[I*bs]), &(b1[I*bs]));
		for ( int k = 0; k < lc; ++k) {
			jacobi_dtrsm_csr( diagL[I], &(ltmp[I * N * bs + k *bs]), &(x0[I * N * bs + k * bs]));
			jacobi_dtrsmt_csr( diagL[I], &(x0[I*N*bs+k*bs]), &(x1[I*N*bs+k*bs]) );
			jacobi_dsubvv( b1, vtmp, I, bs, &(x1[I * N * bs+k*bs]));
		}
	}
#pragma omp taskwait

	for ( int k = 0; k < dim; ++k) {
		x[k] = vtmp[k];
	}

//////////////End of Iter 1///////////////

	iter = 0;
	while(!converge && iter < MAX_ITER) {
		++iter;
//		printf("Iter: %d\n", iter);
		for(int k = 0; k < dim; ++k){
			vtmp[k] = 0.0;
		}
		for (int k = 0; k < N * N * bs; ++k)
			ltmp[k] = 0.0;

		for ( int I = 0; I < N; ++I ) {
			int lc = 0;
			for ( int J = vptr[I]; J < vptr[I+1]; ++J ) {
				if ( vpos[J] != I ){
					Acurrent = vval[J];
					jacobi_dgemv_csr( Acurrent, &(x[vpos[J] * bs]), &(ltmp[I * N * bs + lc * bs]));
					lc++;
				}
			}

			/*
			 *	A[i,i] * x[i] = vtmp[i]
			 */
			jacobi_dtrsm_csr( diagL[I], &(b[I*bs]), &(b0[I*bs]));
			jacobi_dtrsmt_csr( diagL[I], &(b0[I*bs]), &(b1[I*bs]));
			for ( int k = 0; k < lc; ++k) {
				jacobi_dtrsm_csr( diagL[I], &(ltmp[I * N * bs + k *bs]), &(x0[I * N * bs + k * bs]));
				jacobi_dtrsmt_csr( diagL[I], &(x0[I*N*bs+k*bs]), &(x1[I*N*bs+k*bs]) );
				jacobi_dsubvv( b1, vtmp, I, bs, &(x1[I * N * bs+k*bs]));
			}
		}

#pragma omp taskwait

//		converge = 1;
//		for ( int k = 0; k < dim0 && converge; ++k) {
//			if ( x[k] - vtmp[k] ) 
//				converge = 0;
//		}

		for ( int k = 0; k < dim; ++k) {
			x[k] = vtmp[k];
		}
		converge = 0;
	}

	puts("through");
	free(x0); free(x1);
	free(b0); free(b1);
	free(vtmp); free(ltmp);
	for(int k = 0; k < N; ++k){
		free(etree0[k]);
	}

	return 0;
}

#if 0
int jacobi_main_csr(hbmat_t *Acsr, double *x, double *b, int bs){

	int M = Acsr->m; int N = Acsr->n;
	int *vdiag = Acsr->vdiag;
	int *vptr = Acsr->vptr; int *vpos = Acsr->vpos; hbmat_t **vval = Acsr->vval;
	int converge = 0;

	hbmat_t *Adiag, *Adiag0, *Acurrent;
	hbmat_t **diagL = malloc(N * sizeof(hbmat_t*) );

	double *vtmp = calloc(N * bs , sizeof(double));
	double *ltmp = calloc( N * N * bs , sizeof(double));

	double *x0 = malloc(N * bs * sizeof(double));
	double *x1 = malloc(N * bs * sizeof(double));
	
////////////////////1st Iteration//////////////////

	hbmat_t **A0 = malloc(dim *sizeof(hbmat_t*));
	int **etree0 = malloc(dim * sizeof(int*));
	int **work = malloc(N * sizeof(int*));
	for(int k = 0; k < N; ++k)
		work[k] = malloc(bs*sizeof(int));


	for ( int I = 0; I < N; ++I ) {
		etree0[I] = etree(vval[vdiag[I]]);
		A0[I] = hb2hbh_sym_etree_csr_p(vval[vdiag[I]], bs, etree0[I]);
		diagL[I] = ((hbmat_t**)A0[I]->vval)[0];
		potrf_sparse_csr(diagL[I]);
		int lc = 0;
		for ( int J = vptr[I]; J < vptr[I+1]; ++J ) {
			if ( vpos[J] != I ){
				Acurrent = vval[J];
				jacobi_dgemv_csr( Acurrent, &(x[vpos[J] * bs]), &(ltmp[I*N*bs+lc*bs]));
				lc++;
			}
		}
		for(int k = 0; k < lc; ++k)
			jacobi_dsubvv(b,vtmp,I,bs,&(ltmp[I*N*bs+k*bs]));

		/*
		 *	A[i,i] * x[i] = vtmp[i]
		 */
		jacobi_dtrsm_csr( diagL[I], &(vtmp[I * bs]), &(x0[I*bs]));
		jacobi_dtrsmt_csr( diagL[I], &(x0[I*bs]), &(x1[I*bs]) );
	}
#pragma omp taskwait

	for ( int k = 0; k < dim; ++k) {
		x[k] = x1[k];
	}
	puts("Iter 1 through");

//////////////End of Iter 1///////////////

	iter = 0;
	while(!converge && iter < MAX_ITER) {
		++iter;
//		printf("Iter: %d\n", iter);
		for(int k = 0; k < dim; ++k){
			vtmp[k] = 0.0;
		}
		for (int k = 0; k < N * N * bs; ++k)
			ltmp[k] = 0.0;

		for ( int I = 0; I < N; ++I ) {
			int lc = 0;
			for ( int J = vptr[I]; J < vptr[I+1]; ++J ) {
				if ( vpos[J] != I ){
					Acurrent = vval[J];
					jacobi_dgemv_csr( Acurrent, &(x[vpos[J] * bs]), &(ltmp[I * N * bs + lc * bs]));
					lc++;
				}
			}
			for(int k = 0; k < lc; ++k)
				jacobi_dsubvv(b,vtmp,I,bs,&(ltmp[I*N*bs+k*bs]));

			/*
			 *	A[i,i] * x[i] = vtmp[i]
			 */
			jacobi_dtrsm_csr( diagL[I], &(vtmp[I * bs]), &(x0[I*bs]));
			jacobi_dtrsmt_csr( diagL[I], &(x0[I*bs]), &(x1[I*bs]) );
		}

#pragma omp taskwait

//		converge = 1;
//		for ( int k = 0; k < dim0 && converge; ++k) {
//			if ( x[k] - vtmp[k] ) 
//				converge = 0;
//		}

		for ( int k = 0; k < dim; ++k) {
			x[k] = x1[k];
		}
		converge = 0;
	}

	puts("through");
	free(x1); free(x0);

	return 0;
}
#endif
