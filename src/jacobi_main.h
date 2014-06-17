#ifndef __JACOBI_MAIN_H__
#define __JACOBI_MAIN_H__

#include <stdio.h>
#include "hb.h"
#include "hbext.h"
#include "jacobi_kernels.h"

int jacobi_main_csr(hbmat_t *Acsr, double *x, double *b, int bs);

#endif
