#ifndef __JACOBI_SETUP_H__
#define __JACOBI_SETUP_H__

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

#include "matutil/hb.h"
#include "matutil/hbconvrt.h"
#include "matutil/iohb.h"
#include "mkl.h"

int jacobi_setup(int argc, char* argv[]);
void jacobi_shutdown();

#endif 
