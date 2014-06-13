#ifndef __JACOBI_SETUP_H__
#define __JACOBI_SETUP_H__

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
//#include "ompss_cholesky.h"

#include "hb.h"
#include "hbconvrt.h"
#include "iohb.h"
#include "mkl.h"

int jacobi_setup(int argc, char* argv[]);
void jacobi_shutdown();

#endif 
