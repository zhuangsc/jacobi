#ifndef __CHOLS_KERNELS_H__
#define __CHOLS_KERNELS_H__


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hb.h"
#include "vector.h"
#include "hbconvrt.h"

#include "mkl.h"


static inline void __attribute__((always_inline)) array_d2s(hbmat_t* A, double* peel, int col) {
	int* vptr = A->vptr;
	int* vpos = A->vpos;
	double* vval = A->vval;
	int i;
	for ( i = vptr[col]; i < vptr[col+1]; i++ ) {
		vval[i] = peel[vpos[i]];
	}
}

static inline void __attribute__((always_inline)) array_s2d(hbmat_t* A, double* peel, int col) {
	int* vptr = A->vptr;
	int* vpos = A->vpos;
	double* vval = A->vval;
	int i;
	for ( i = vptr[col]; i < vptr[col+1]; i++ ) {
		peel[vpos[i]] = vval[i];
	}
}

static inline void __attribute__((always_inline)) array_clear(double* peel, int elemc) {
	int i;
	for ( i=0; i<elemc; i++ ) {
		peel[i] = 0;
	}
}


void test_print_matrix(const hbmat_t*, int, char*);
void hbcopy(hbmat_t*, hbmat_t*);

#pragma omp task inout([1]A) priority(1)
void potrf_sparse(hbmat_t* A);

#pragma omp task in([1]A) inout([1]B) priority(2)
void dsyrk_sparse(hbmat_t* A, hbmat_t* B);

#pragma omp task in([1]A, [1]B) inout([1]C)
void dgemm_sparse(hbmat_t* A, hbmat_t* B, hbmat_t* C);

#pragma omp task in([1]A) inout([1]B)
void dtrsm_sparse(hbmat_t* A, hbmat_t* B);

#pragma omp task inout([1]A) priority(4)
void potrf_sparse_csr(hbmat_t* A);

#pragma omp task in([1]A) inout([1]B) priority(3)
void dsyrk_sparse_csr(hbmat_t* A, hbmat_t* B);

#pragma omp task in([1]A, [1]B) inout([1]C)
void dgemm_sparse_csr(hbmat_t* A, hbmat_t* B, hbmat_t* C);

#pragma omp task in([1]A) inout([1]B)
void dtrsm_sparse_csr(hbmat_t* A, hbmat_t* B);

void dgemv_sparse_csr(hbmat_t* A, double *x, double *b);

#endif
