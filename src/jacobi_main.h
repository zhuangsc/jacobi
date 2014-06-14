#ifndef __JACOBI_MAIN_H__
#define __JACOBI_MAIN_H__

#include <stdio.h>
#include "jacobi_kernels.h"
#include "ompss_cholesky.h"

int jacobi_main_csr(hbmat_t *Acsr, double *x, double *b, int bs);

#endif
