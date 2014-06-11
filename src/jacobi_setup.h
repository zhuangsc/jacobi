#ifndef __JACOBI_SETUP_H__
#define __JACOBI_SETUP_H__

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include "ompss_cholesky.h"

extern int readHB_newmat_double(const char* filename, int* M, int* N, int* nonzeros,
		int** colptr, int** rowind, double** val);
extern void one2zero(hbmat_t*);
extern void hb_free(hbmat_t*);
extern void hbh_free(hbmat_t*);
extern void hbh_free2(hbmat_t*);
extern hbmat_t* hbh2hb(hbmat_t*);

int bs;

#endif 
