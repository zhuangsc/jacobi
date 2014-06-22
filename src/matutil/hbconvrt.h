#ifndef __HBCONVRT_H__
#define __HBCONVRT_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#include "hb.h"
#include "vector.h"


void hbh_free(hbmat_t *A);
void hb_free(hbmat_t *A);
void hbh_free2(hbmat_t *A);
hbmat_t* hb2csr(hbmat_t *A);
hbmat_t* hbh2hb(hbmat_t *A);
hbmat_t* hb_transpose(hbmat_t *A);
hbmat_t* hbh_hyper_transpose(hbmat_t *A);
hbmat_t* hb2hbh_sym_etree_u(hbmat_t *A, int b, int* etree);
hbmat_t* hb2hbh_sym_etree(hbmat_t *A, int b, int* etree);
hbmat_t *hb2hbb(hbmat_t *A, int b);
hbmat_t* hbb2csrb(hbmat_t *A);

//#pragma omp task in([1]A) out([1]entry, [1]block)
void hb2hbh_csr_task(int I, int J, hbmat_t *A, int b, int *entry, hbmat_t *block);

//#pragma omp task in([1]A) out([1]entry, [1]block)
void hb2hbh_csc_task(int I, int J, hbmat_t *A, int b, int *entry, hbmat_t *block);

hbmat_t* hb2hbh(hbmat_t *A, int b, int is_csr);

//hbmat_t* ll2b(hbmatm_t *A);
//hbmatm_t* b2ll(hbmat_t *A);

//hbmatm_t* hbbm2csrbm(hbmatm_t *A);
//void llsetdiag( hbmatm_t *A );
//void m_sync(hbmatm_t *Ahb, hbmatm_t *Acsr);

void ereach_csr_p(hbmat_t *A, int r, int *etree, vector_t* sub_row, vector_t* sub_val);

//#pragma omp task in([1]A, [1]etree) out([1]entry, [1]block)
void symbolic_csr_task(int I, int J, hbmat_t *A, int b, int *etree, int *entry, hbmat_t *block);

hbmat_t* hb2hbh_sym_etree_csr_p(hbmat_t *A, int b, int *etree);


//#pragma omp task in([1]etree) out([1]entry)
void hyper_sym_csr_task0(int I, int J, hbmat_t *A, int b, int *etree, int *entry);

void hyper_sym_csr_task1(hbmat_t *block);

#define FILLINS 1
pthread_mutex_t mutexhb;
int *vptr_pool, *vpos_pool;
double *vval_pool;
int vptr_unit, vpos_unit, vval_unit;
int vptr_pp, vpos_pp, vval_pp;
void hyper_sym_csr_task2(hbmat_t *block);

hbmat_t* hb2hbh_hyper_sym_csr(hbmat_t *A, int b, int *etree);

#endif // __HBCONVRT_H__
