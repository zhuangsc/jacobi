
#include "jacobi_main.h"

int jacobi_main_csr(hbmat_t *Acsr, double *x, double *b, double *vtmp, int bs){

	int m = Acsr->m; int n = Acsr->n;
	int *vdiag = Acsr->vdiag;
	int *vptr = Acsr->vptr; int *vpos = Acsr->vpos; hbmat_t **vval = Acsr->vval;

	hbmat_t *Adiag, *Acurrent;

	for ( int J = 0; J < n; ++J ) {
		Adiag = vval[vdiag[J]];
		for ( int I = vptr[J]; I < vptr[J+1]; ++I ) {
			if ( vpos[I] == J )
				continue;
			Acurrent = vval[I];
			jacobi_dgemv_csr( Acurrent, x, vtmp, vpos[I], bs );
		}
		jacobi_dsubvv( b, vtmp, J, bs );
	}

}
