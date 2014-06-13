#include "hbconvrt.h"

hbmat_t* hb_transpose(hbmat_t *A){
	hbmat_t *B = (hbmat_t*) malloc(sizeof(hbmat_t));
	int m = A->m;
	int n = A->n;
	int elemc = A->elemc;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	double *vval = A->vval;
	int acc = 0;

	B->m = n; B->n = m; B->elemc = elemc;
	vector_t *trans_vptr, *trans_vpos, *trans_vval;
	trans_vptr = vector_create_size(n);
	trans_vpos = vector_create_size(elemc);
	trans_vval = vector_create_size(elemc);
	vector_clear(trans_vptr); 
	vector_clear(trans_vpos);
	vector_clear(trans_vval);
	vel_t vptr_vel, vpos_vel, vval_vel;

	for(int j = 0; j < m; j++){
		vptr_vel.i = trans_vpos->elemc + 1;
		vector_insert(trans_vptr, vptr_vel);
		for(int i = 0; i < n; i++){
			for(int k = vptr[i]; k < vptr[i+1]; k++){
				if (j == vpos[k-1] - 1){
					acc++;
				 	vpos_vel.i = i + 1;
					vector_insert(trans_vpos, vpos_vel);
					vval_vel.d = vval[k-1];
					vector_insert(trans_vval, vval_vel);
					break;
				}
			}
		}
	}
	vptr_vel.i = trans_vpos->elemc + 1;
	vector_insert(trans_vptr, vptr_vel);

	B->vptr = vector2int(trans_vptr);
	B->vpos = vector2int(trans_vpos);
	B->vval = vector2double(trans_vval);
	//printf("hb_transpose\naddr_a %p \t addr_b %p\n", A, B);
	//free(A);
	return B;
}

hbmat_t* hbh_hyper_transpose(hbmat_t *A){
	hbmat_t *B = (hbmat_t*) malloc(sizeof(hbmat_t));
	int m = A->m;
	int n = A->n;
	int elemc = A->elemc;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	hbmat_t **vval = A->vval;
	int acc = 0;

	B->m = n; B->n = m; B->elemc = elemc;
	hbmat_t **b_vval = (hbmat_t**) malloc(elemc*sizeof(hbmat_t*));
	vector_t *trans_vptr, *trans_vpos;
	trans_vptr = vector_create_size(n);
	trans_vpos = vector_create_size(elemc);
	vector_clear(trans_vptr); 
	vector_clear(trans_vpos);
	vel_t vptr_vel, vpos_vel;

	for(int j = 0; j < m; j++){
		vptr_vel.i = trans_vpos->elemc + 1;
		vector_insert(trans_vptr, vptr_vel);
		for(int i = 0; i < n; i++){
			for(int k = vptr[i]; k < vptr[i+1]; k++){
				if (j == vpos[k-1] - 1){
				 	vpos_vel.i = i + 1;
					vector_insert(trans_vpos, vpos_vel);
					b_vval[acc++] = vval[k-1];
					break;
				}
			}
		}
	}
	vptr_vel.i = trans_vpos->elemc + 1;
	vector_insert(trans_vptr, vptr_vel);
	B->vptr = vector2int(trans_vptr);
	B->vpos = vector2int(trans_vpos);
	B->vval = b_vval;

	//printf("hyper-transpose\n");
	for(int i = 0; i < B->elemc; i++){
		//((hbmat_t**)B->vval)[i] = hb_transpose(((hbmat_t**)B->vval)[i]);
		//printf("%p\t%p\n", ((hbmat_t**)B->vval)[i], b_vval[i]);
		((hbmat_t**)B->vval)[i] = hb_transpose(b_vval[i]);
		//printf("%p\t%p\n", ((hbmat_t**)B->vval)[i], b_vval[i]);
		//printf("---------------------------------\n");
	}

	return B;
}

hbmat_t* hb2hbh_sym_etree_u(hbmat_t *A, int b, int* etree){

	hbmat_t *Ab = (hbmat_t*) malloc(sizeof(hbmat_t)); 

	int m = A->m;
	int n = A->n;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	double *vval = A->vval;

	int M = ( m + b - 1 ) / b;
	int N = ( n + b - 1 ) / b;
	int num = ((1 + M) * N) / 2;  //total number of blocks in the lower triangular matrix

	Ab->m = M; Ab->n = N; Ab->elemc = 0; Ab->vdiag = NULL;
	Ab->vptr = malloc((N+1)*sizeof(int));
	Ab->vpos = malloc(num*sizeof(int));
	Ab->vval = malloc(num*sizeof(hbmat_t*));

	hbmat_t* acchb = malloc(num*sizeof(hbmat_t));
	int acc = 0 ;
	int vpos_count = 0 ;

	if ( M==0 || N==0 ) {
		fprintf( stderr, "block size %i too large\n", b);
	}

	for ( int J = 0; J < N; ++J ) { 	
		
		int jstart = J * b;
		int jc = n - jstart;
		jc = b;
		Ab->vptr[J] = acc;

		for ( int I = 0; I < J+1; ++I ) { 	

			int base_col = J * b ;
			int base_row = I * b ;
			vel_t vptr_current;
			vector_t* sub_tree = vector_create();
			vector_t* sub_vptr = vector_create();
			vector_t* sub_vpos = vector_create();
			vector_t* sub_vval = vector_create();
			vector_clear(sub_vptr); vector_clear(sub_vpos); vector_clear(sub_vval);

			for ( int j = 0 ; j < jc; j++ ) { 	// Innermost loop: column by column

				int pos_col = base_col + j ;	// Absolute column position (0 based)
				int max_row = pos_col ;		// Maximum row position in this column (0 based)
				int bborder = base_row + b - 1;
				int current_row ;
				max_row = max_row <= bborder ? max_row : bborder ;

				int min_nz = vptr[pos_col] ;
				int	max_nz = vptr[pos_col+1] ;

				vptr_current.i = sub_vpos->elemc;
				vector_insert(sub_vptr, vptr_current);
				vector_clear(sub_tree);			

				//Padding
				if (pos_col >= n){
					if(base_row+b > n){
						vptr_current.i = j;
						vector_insert(sub_vpos, vptr_current);
						vptr_current.d = 0;
						vector_insert(sub_vval, vptr_current);
					}
					continue;
				}

				for ( int current = min_nz; current < max_nz; ++current ) {
					int status ;
					current_row = vpos[current];
					vel_t vel_current = {i : current_row} ;

					if (current_row >= base_row && current_row < max_row){
						status = vector_insert_t(sub_tree, vel_current) ;
						if (!status)
							continue;
					}
					else if (current_row == max_row){
						status = vector_insert_t(sub_tree, vel_current) ; 
						break ; //current row is the boarder, break the loop
					}
					else if (current_row > max_row){
						break ;	//current row is out of the boarder, break the loop
					}

					//Traverse the elimination tree
					for (int node = etree[current_row]; node != -1 && node <= max_row; node = etree[node]){
						if (node >= base_row){
							vel_t vel_node = {i : node} ;
							vector_insert_t(sub_tree, vel_node) ;
						}
					}
				}
				
				//Sort the vector
				vector_qsorti(sub_tree);
				for(int i = 0; i < sub_tree->elemc; i++){
					vel_t vval_current = {d : 0.0};
					vel_t vpos_c;
					int fill_in = 1;

					vpos_c.i = sub_tree->elem[i].i - base_row;
					vector_insert(sub_vpos, vpos_c);
					for (int current = min_nz; current < max_nz; current++){
						current_row = vpos[current];
						if (current_row > max_row) break;
						if (current_row == sub_tree->elem[i].i){
							vval_current.d = vval[current];
							vector_insert(sub_vval, vval_current);
							fill_in = 0;
							break;
						}
					}
					if (fill_in)
						vector_insert(sub_vval, vval_current);
				}
			}
			
			vector_free(sub_tree);
			if (sub_vpos->elemc != 0){
				acchb[acc].m = jc;
				acchb[acc].n = jc;
				acchb[acc].elemc = sub_vpos->elemc;
				acchb[acc].vdiag = NULL;
				vptr_current.i = sub_vpos->elemc;
				vector_insert(sub_vptr, vptr_current);
//				vector_printi(sub_vptr);
//				vector_printi(sub_vpos);
//				vector_printd(sub_vval);
				acchb[acc].vptr = vector2int(sub_vptr);
				acchb[acc].vpos = vector2int(sub_vpos);
				acchb[acc].vval = vector2double(sub_vval);
				Ab->vpos[vpos_count] = I;
				((hbmat_t**)Ab->vval)[vpos_count] = &(acchb[acc]);
				//((hbmat_t**)Ab->vval)[vpos_count] = acchb+acc;
				vpos_count++;
				Ab->elemc++;
				acc++;
			}
			else{
				vector_free(sub_vptr);
				vector_free(sub_vpos);
				vector_free(sub_vval);
			}
		}
	}
	Ab->vptr[N] = Ab->elemc ;
	
	//TODO Check the correctness
	Ab->vpos = (int*) realloc(Ab->vpos, acc * sizeof(int));
	Ab->vval = (hbmat_t**) realloc(Ab->vval, acc * sizeof(hbmat_t*));

	return Ab;
}


	void ereach_csr(hbmat_t *A, int r, int *etree, vector_t* sub_row){
	int m = A->m; int n = A->n; int elemc = A->elemc;
	int *vptr = A->vptr; int *vpos = A->vpos; double *vval = A->vval;

	vector_clear(sub_row);
	vel_t ab_vel;

	/*
	 * Finding the non-zero pattern
	 */
	for ( int k = vptr[r]; k < vptr[r+1]; ++k ) {
		ab_vel.i = vpos[k];
		vector_insert_t(sub_row, ab_vel);
		int inserted = 1;
		ab_vel.i = etree[ab_vel.i];
		while ( inserted && ab_vel.i > 0 && ab_vel.i <= r) {
			inserted = vector_insert_t(sub_row, ab_vel);
			ab_vel.i = etree[ab_vel.i];
		}
	}
	vector_qsorti(sub_row); // Sort vpos
}

hbmat_t* hb2hbh_sym_etree_csr(hbmat_t *A, int b, int* etree){

	int m = A->m; int n = A->n; int elemc = A->elemc;
	int *vptr = A->vptr; int *vpos = A->vpos; double* vval = A->vval;
	int M = ( m + b - 1 ) / b;
	int N = ( n + b - 1 ) / b;
	int num = ((1 + M) * N) / 2;

	hbmat_t* hyper = malloc(sizeof(hbmat_t));
	hyper->m = M; hyper->n = N; hyper->vdiag = NULL;
	hyper->vval = (hbmat_t**) malloc(num * sizeof(hbmat_t*));

	hbmat_t* acchb = malloc(num*sizeof(hbmat_t));
	int acc = 0;

	vector_t* ab_vptr = vector_create(); 
	vector_t* ab_vpos = vector_create();
	vector_clear(ab_vptr); vector_clear(ab_vpos);

	vel_t pos_val;

	if ( M==0 || N==0 ) {
		fprintf( stderr, "block size %i too large\n", b);
	}

	vector_t **vvptr = malloc(M * sizeof(vector_t*));
	vector_t **vvpos = malloc(M * sizeof(vector_t*));
	vector_t **vvval = malloc(M * sizeof(vector_t*));
	for( int i = 0; i < M; ++i){
		vvptr[i] = vector_create();
		vvpos[i] = vector_create();
		vvval[i] = vector_create();
	}

	vector_t *vpos_tmp = vector_create();
	double *vval_tmp = malloc(m * sizeof(double));
	for ( int I = 0; I < M; ++I ) {
		pos_val.i = ab_vpos->elemc;
		vector_insert(ab_vptr, pos_val);

		int istart = I*b; int iend = (I+1)*b;  //Absolute row position
		int blocks = I + 1; //column blocks
		for ( int J = 0; J < blocks; ++J ) {
			vector_clear(vvptr[J]);
			vector_clear(vvpos[J]);
			vector_clear(vvval[J]);
			pos_val.i = 0;
			vector_insert(vvptr[J], pos_val);
		}
		for ( int i = istart; i < iend; ++i ){  //row by row
			/*
			 * Padding
			 */
			if ( i >= m ) {
				pos_val.i = i - (blocks - 1) * b;
				vector_insert(vvpos[blocks-1], pos_val);
				pos_val.d = 1;
				vector_insert(vvval[blocks-1], pos_val);
				for ( int ll = 0; ll < blocks; ++ll) {
					pos_val.i = vvpos[ll]->elemc;
					vector_insert(vvptr[ll], pos_val);
				}
				continue;
			}

			ereach_csr(A, i, etree, vpos_tmp);
			for ( int tmp = 0; tmp < vpos_tmp->elemc; ++tmp){
				vval_tmp[tmp] = 0;
			}
			int ptr = 0;
			/*
			 * Setup vval properly
			 */
			for ( int j = vptr[i]; j < vptr[i+1]; ++j ){
				while(vpos_tmp->elem[ptr].i < vpos[j])
					++ptr;
				if (vpos_tmp->elem[ptr].i == vpos[j]){
					vval_tmp[ptr] = vval[j];
					++ptr;
				}
			}
			ptr = 0;
			/*
			 * Update row
			 */
			for ( int J = 0; J < blocks; ++J ) {
				int jstart = J * b; 
				int jend = (J + 1) * b;
				while ( vpos_tmp->elem[ptr].i < jend && ptr < vpos_tmp->elemc ) {
					pos_val.i = vpos_tmp->elem[ptr].i - jstart;
					vector_insert(vvpos[J], pos_val);
					pos_val.d = vval_tmp[ptr];
					vector_insert(vvval[J], pos_val);
					++ptr;
				}
				pos_val.i = vvpos[J]->elemc;
				vector_insert(vvptr[J], pos_val);
			}
		}
		for (int J = 0; J < blocks; ++J ) {
			if ( !vvpos[J]->elemc )
				continue;
			int vval_ptr;
			acchb[acc].m = b; acchb[acc].n = b; acchb[acc].elemc = vvpos[J]->elemc;
			acchb[acc].vdiag = NULL;
			acchb[acc].vptr = vector2int_nf(vvptr[J]);
			acchb[acc].vpos = vector2int_nf(vvpos[J]);
			acchb[acc].vval = vector2double_nf(vvval[J]);

			vval_ptr = ab_vpos->elemc;
			pos_val.i = J;
			vector_insert(ab_vpos, pos_val);
			((hbmat_t**)hyper->vval)[vval_ptr] = &(acchb[acc]);
			++acc;
		}
	}
	pos_val.i = ab_vpos->elemc;
	vector_insert(ab_vptr, pos_val);
	hyper->elemc = ab_vpos->elemc;
	hyper->vptr = vector2int(ab_vptr);
	hyper->vpos = vector2int(ab_vpos);

	vector_free(vpos_tmp); free(vval_tmp);
	for( int i = 0; i < M; ++i){
		vector_free(vvptr[i]);
		vector_free(vvpos[i]);
		vector_free(vvval[i]);
	}
	free(vvptr); free(vvpos); free(vvval);

	return hyper;
}




void ereach_csr_p(hbmat_t *A, int r, int *etree, vector_t *sub_row, vector_t *sub_val){
	int m = A->m; int n = A->n; int elemc = A->elemc;
	int *vptr = A->vptr; int *vpos = A->vpos; double *vval = A->vval;

	vector_clear(sub_row);
	vector_clear(sub_val);
	vel_t ab_vel;

	/*
	 * Finding the non-zero pattern
	 */
	for ( int k = vptr[r]; k < vptr[r+1]; ++k ) {
		ab_vel.i = vpos[k];
		vector_insert_t(sub_row, ab_vel);
		int inserted = 1;
		ab_vel.i = etree[ab_vel.i];
		while ( inserted && ab_vel.i > 0 && ab_vel.i <= r) {
			inserted = vector_insert_t(sub_row, ab_vel);
			ab_vel.i = etree[ab_vel.i];
		}
	}
	vector_qsorti(sub_row); // Sort vpos

	int ccol = vptr[r];
	for ( int i = 0; i < sub_row->elemc; ++i ) {
		while ( vpos[ccol] < sub_row->elem[i].i && ccol < vptr[r+1] )
			++ccol;
		if ( vpos[ccol] == sub_row->elem[i].i ){
			ab_vel.d = vval[ccol];
			vector_insert(sub_val, ab_vel);
			continue;
		}
		if (vpos[ccol] > sub_row->elem[i].i ) {
			ab_vel.d = 0;
			vector_insert(sub_val, ab_vel);
		}
	}
}

void symbolic_csr_task(int I, int J, hbmat_t *A, int b, int *etree, int *entry, hbmat_t *block){
	int m = A->m;
	int* vptr = A->vptr; int* vpos = A->vpos; 
	double* vval = A->vval;
	int brow = I*b; int erow = (I+1)*b;
	int bcol = J*b; int ecol = (J+1)*b;

	vector_t* ab_vptr = vector_create();
	vector_t* ab_vpos = vector_create();
	vector_t* ab_vval = vector_create();
	vector_clear(ab_vptr); vector_clear(ab_vpos); vector_clear(ab_vval);
	vel_t pos_val;

	*entry = 0;

	for ( int L = brow; L < erow; ++L ) {
		pos_val.i = ab_vpos->elemc; 
		vector_insert(ab_vptr, pos_val);
		
		//Padding
		if ( L >= m ) {
			if ( ecol < m )
				continue;
			pos_val.i = L - bcol;
			vector_insert(ab_vpos, pos_val);
			pos_val.d = 1;
			vector_insert(ab_vval, pos_val);
			continue;
		}

		int p_elemc = ab_vpos->elemc;
		for ( int k = vptr[L]; k < vptr[L+1]; ++k ) {

			if ( vpos[k] >= bcol && vpos[k] < ecol){
				pos_val.i = vpos[k] - bcol;
				vector_insert_t_partial(ab_vpos, pos_val, p_elemc);
			}
			int inserted = 1;
			pos_val.i = etree[vpos[k]];
			while ( inserted && pos_val.i < ecol && pos_val.i <= L) {
				if ( pos_val.i >= bcol ) {
					vel_t tmp;
					tmp.i = pos_val.i - bcol;
					inserted = vector_insert_t_partial(ab_vpos, tmp, p_elemc);
				}
				pos_val.i = etree[pos_val.i];
			}
		}

		//Sorting vpos
		vector_qsorti_partial(ab_vpos, p_elemc);

		int ccol = vptr[L];
		for ( int i = p_elemc; i < ab_vpos->elemc; ++i ) {
			while ( (vpos[ccol] - bcol) < ab_vpos->elem[i].i && ccol < vptr[L+1])
				++ccol;

			if ( (vpos[ccol] - bcol) == ab_vpos->elem[i].i ){
				pos_val.d = vval[ccol];
				vector_insert(ab_vval, pos_val);
				continue;
			}
			else {//if ((vpos[ccol] - bcol) > ab_vpos->elem[i].i ) {
				pos_val.d = 0;
				vector_insert(ab_vval, pos_val);
			}
		}
	}

	pos_val.i = ab_vpos->elemc;
	vector_insert(ab_vptr, pos_val);

	if ( ab_vpos->elemc ){
		block->m = b; block->n = b; block->elemc = ab_vpos->elemc;
		block->vdiag = NULL;
//		vector_printi(ab_vptr);
//		vector_printi(ab_vpos);
//		vector_printd(ab_vval);
		block->vptr = vector2int(ab_vptr); 
		block->vpos = vector2int(ab_vpos);
		block->vval = vector2double(ab_vval);
		*entry = 1;
	}
}

hbmat_t* hb2hbh_sym_etree_csr_p(hbmat_t *A, int b, int *etree){

	int m = A->m; int n = A->n; int elemc = A->elemc;
	int *vptr = A->vptr; int *vpos = A->vpos; double* vval = A->vval;
	int M = ( m + b - 1 ) / b;
	int N = ( n + b - 1 ) / b;
	int num = ((1 + M) * N) / 2;

	hbmat_t* hyper = malloc(sizeof(hbmat_t));
	hyper->m = M; hyper->n = N; hyper->vdiag = NULL;
	hyper->vval = malloc(num * sizeof(hbmat_t*));
	hbmat_t** hbmat_array = malloc(num * sizeof(hbmat_t*));
	int* hentry = malloc(num * sizeof(int));

	vector_t* ab_vptr = vector_create(); 
	vector_t* ab_vpos = vector_create();
	vector_clear(ab_vptr); vector_clear(ab_vpos);
	vel_t pos_val;

	if ( M==0 || N==0 ) {
		fprintf( stderr, "block size %i too large\n", b);
	}

	for(int i = 0; i < num; ++i)
		hbmat_array[i] = malloc(sizeof(hbmat_t));
	int acc = 0;
	int I, J;
	for ( I = 0; I < M; ++I ){
		for ( J = 0; J < I+1; ++J){
			symbolic_csr_task(I, J, A, b, etree, &(hentry[acc]), hbmat_array[acc]);
			++acc;
		}
	}

#pragma omp taskwait

	acc = 0;
	int acc0 = 0;
	for ( I = 0; I < M; ++I ) {
		pos_val.i = ab_vpos->elemc;
		vector_insert(ab_vptr, pos_val);
		for ( J = 0; J < I+1; ++J ) {
			if ( hentry[acc] ) {
				pos_val.i = J;
				vector_insert(ab_vpos, pos_val);
				((hbmat_t**)hyper->vval)[acc0] = hbmat_array[acc];
				++acc;
				++acc0;
			} else {
				free(hbmat_array[acc]);
				++acc;
			}
		}
	}

	pos_val.i = ab_vpos->elemc;
	vector_insert(ab_vptr, pos_val);
	hyper->elemc = ab_vpos->elemc;
	hyper->vptr = vector2int(ab_vptr);
	hyper->vpos = vector2int(ab_vpos);

	return hyper;
}

void hyper_sym_csr_task0(int I, int J, hbmat_t *A, int b, int *etree, int *entry){
	int m = A->m;
	int* vptr = A->vptr; int* vpos = A->vpos; 
	int brow = I*b; int erow = (I+1)*b;
	int bcol = J*b; int ecol = (J+1)*b;
	int ccol;

	*entry = 0;

	/*
	 * TODO Add path convergence check if necessary
	 */
	if (erow >= m)
		erow = m;
	for ( int L = brow; L < erow; ++L ) {
		for ( int k = vptr[L]; k < vptr[L+1]; ++k ) {
			ccol = vpos[k];
			while ( ccol < ecol && ccol <= L ) {
				if ( ccol >= bcol && ccol < ecol){
					*entry = 1;
					return;
				}
				ccol = etree[ccol];
			}
		}
	}
}

void hyper_sym_csr_task1(hbmat_t *block){

	/*
	 * Extract info from the block matrix
	 */
	int b = block->m;
	int I = block->orig_row; int J = block->orig_col;
	hbmat_t* A = block->orig;
	int* etree = block->e_tree;

	/*
	 * Setup variables
	 */
	int m = A->m;
	int* vptr = A->vptr; int* vpos = A->vpos; 
	double* vval = A->vval;
	int brow = I*b; int erow = (I+1)*b;
	int bcol = J*b; int ecol = (J+1)*b;

	vector_t* ab_vptr = vector_create();
	vector_t* ab_vpos = vector_create();
	vector_t* ab_vval = vector_create();
	vector_clear(ab_vptr); vector_clear(ab_vpos); vector_clear(ab_vval);
	vel_t pos_val;

	for ( int L = brow; L < erow; ++L ) {
		pos_val.i = ab_vpos->elemc; 
		vector_insert(ab_vptr, pos_val);
		
		/*
		 * Column padding
		 */
		if ( L >= m ) {
			if ( ecol < m )
				continue;
			pos_val.i = L - bcol;
			vector_insert(ab_vpos, pos_val);
			pos_val.d = 1;
			vector_insert(ab_vval, pos_val);
			continue;
		}

		int p_elemc = ab_vpos->elemc;
		for ( int k = vptr[L]; k < vptr[L+1]; ++k ) {

			if ( vpos[k] >= bcol && vpos[k] < ecol){
				pos_val.i = vpos[k] - bcol;
				vector_insert_t_partial(ab_vpos, pos_val, p_elemc);
			}
			int inserted = 1;
			pos_val.i = etree[vpos[k]];
			while ( inserted && pos_val.i < ecol && pos_val.i <= L) {
				if ( pos_val.i >= bcol ) {
					vel_t tmp;
					tmp.i = pos_val.i - bcol;
					inserted = vector_insert_t_partial(ab_vpos, tmp, p_elemc);
				}
				pos_val.i = etree[pos_val.i];
			}
		}

		//Sorting vpos
		vector_qsorti_partial(ab_vpos, p_elemc);

		int ccol = vptr[L];
		for ( int i = p_elemc; i < ab_vpos->elemc; ++i ) {
			while ( (vpos[ccol] - bcol) < ab_vpos->elem[i].i && ccol < vptr[L+1])
				++ccol;

			if ( (vpos[ccol] - bcol) == ab_vpos->elem[i].i ){
				pos_val.d = vval[ccol];
				vector_insert(ab_vval, pos_val);
				continue;
			}
			else {		//if ((vpos[ccol] - bcol) > ab_vpos->elem[i].i ) {
				pos_val.d = 0;
				vector_insert(ab_vval, pos_val);
			}
		}
	}

	pos_val.i = ab_vpos->elemc;
	vector_insert(ab_vptr, pos_val);

	if ( ab_vpos->elemc ){
		block->m = b; block->n = b; block->elemc = ab_vpos->elemc;
		block->vdiag = NULL;
//		vector_printi(ab_vptr);
//		vector_printi(ab_vpos);
//		vector_printd(ab_vval);
		block->vptr = vector2int(ab_vptr); 
		block->vpos = vector2int(ab_vpos);
		block->vval = vector2double(ab_vval);
	}else
		printf("Warning! task1 fail. I %d J %d\n", I, J);
}

hbmat_t* hb2hbh_hyper_sym_csr(hbmat_t *A, int b, int *etree){

	int m = A->m; int n = A->n; int elemc = A->elemc;
	int *vptr = A->vptr; int *vpos = A->vpos; double* vval = A->vval;
	int M = ( m + b - 1 ) / b;
	int N = ( n + b - 1 ) / b;
	int num = ((1 + M) * N) / 2;
	int helemc;

	hbmat_t* hyper = malloc(sizeof(hbmat_t));
	hyper->m = M; hyper->n = N; hyper->vdiag = NULL;
	hyper->vval = malloc(num*sizeof(hbmat_t*));
	hyper->e_tree = etree;
//	hbmat_t** hbmat_array = malloc(num*sizeof(hbmat_t*)); 
	int* hentry = malloc(num * sizeof(int));

	vector_t* ab_vptr = vector_create(); 
	vector_t* ab_vpos = vector_create();
	vector_clear(ab_vptr); vector_clear(ab_vpos);
	vel_t pos_val;

	if ( M==0 || N==0 ) {
		fprintf( stderr, "block size %i too large\n", b);
	}


	int acc = 0;
	int I, J;
	for ( I = 0; I < M; ++I ){
		for ( J = 0; J < I+1; ++J){
			hyper_sym_csr_task0(I, J, A, b, etree, &(hentry[acc]));
			++acc;
		}
	}

//#pragma omp taskwait

	acc = 0;
	int acc0 = 0;
	for ( I = 0; I < M; ++I ) {
		pos_val.i = ab_vpos->elemc;
		vector_insert(ab_vptr, pos_val);
		for ( J = 0; J < I+1; ++J ) {
//			printf("I = %d, J = %d, hentry[%d] = %d\n", I, J, acc, hentry[acc]);
			if ( hentry[acc] ) {
				pos_val.i = J;
				vector_insert(ab_vpos, pos_val);
				((hbmat_t**)hyper->vval)[acc0] = malloc(sizeof(hbmat_t));
				((hbmat_t**)hyper->vval)[acc0]->m = b;
				((hbmat_t**)hyper->vval)[acc0]->n = b;
				((hbmat_t**)hyper->vval)[acc0]->vptr = NULL;
				((hbmat_t**)hyper->vval)[acc0]->vpos = NULL;
				((hbmat_t**)hyper->vval)[acc0]->vval = NULL;
				((hbmat_t**)hyper->vval)[acc0]->orig = A;
				((hbmat_t**)hyper->vval)[acc0]->orig_row = I;
				((hbmat_t**)hyper->vval)[acc0]->orig_col = J;
				((hbmat_t**)hyper->vval)[acc0]->e_tree = etree;
				++acc0;
			}
			++acc;
		}
	}

	helemc = ab_vpos->elemc;

	vptr_unit = b + 1;
	vpos_unit = ceil(b * b * FILLINS);
	vval_unit = vpos_unit;
	int num_vptr = vptr_unit * helemc;
	int num_vpos = vpos_unit * helemc;
	int num_vval = num_vpos;
	vptr_pool = malloc(num_vptr * sizeof(int));
	vpos_pool = malloc(num_vpos * sizeof(int));
	vval_pool = malloc(num_vval * sizeof(double));
	vptr_pp = 0; vpos_pp = 0; vval_pp = 0;
	
	pos_val.i = helemc;
	vector_insert(ab_vptr, pos_val);
	hyper->elemc = helemc;
	hyper->vptr = vector2int(ab_vptr);
	hyper->vpos = vector2int(ab_vpos);

	free(hentry);

	return hyper;
}

void hyper_sym_csr_task2(hbmat_t *block){

	/*
	 * Extract info from the block matrix
	 */
	int b = block->m;
	int I = block->orig_row; int J = block->orig_col;
	hbmat_t* A = block->orig;
	int* etree = block->e_tree;

	/*
	 * Setup variables
	 */
	int m = A->m;
	int* vptr = A->vptr; int* vpos = A->vpos; 
	double* vval = A->vval;
	int brow = I*b; int erow = (I+1)*b;
	int bcol = J*b; int ecol = (J+1)*b;

	pthread_mutex_lock(&mutexhb);
	block->vptr = vptr_pool + vptr_pp * vptr_unit;
	block->vpos = vpos_pool + vpos_pp * vpos_unit;
	block->vval = vval_pool + vval_pp * vval_unit;
	++vptr_pp; ++vpos_pp; ++vval_pp;
	pthread_mutex_unlock(&mutexhb);

	vector_int* ab_vptr = vector_int_create(block->vptr, vptr_unit);
	vector_int* ab_vpos = vector_int_create(block->vpos, vpos_unit);
	vector_double* ab_vval = vector_double_create(block->vval, vval_unit);
	vector_int_clear(ab_vptr); vector_int_clear(ab_vpos); vector_double_clear(ab_vval);
	int val_int; double val_fp;

	for ( int L = brow; L < erow; ++L ) {
		val_int = ab_vpos->elemc; 
		vector_int_insert(ab_vptr, val_int);
		
		/*
		 * Column padding
		 */
		if ( L >= m ) {
			if ( ecol < m )
				continue;
			val_int = L - bcol;
			vector_int_insert(ab_vpos, val_int);
			val_fp = 1;
			vector_double_insert(ab_vval, val_fp);
			continue;
		}

		int p_elemc = ab_vpos->elemc;
		for ( int k = vptr[L]; k < vptr[L+1]; ++k ) {

			if ( vpos[k] >= bcol && vpos[k] < ecol){
				val_int = vpos[k] - bcol;
				vector_int_insert_t_partial(ab_vpos, val_int, p_elemc);
			}
			int inserted = 1;
			val_int = etree[vpos[k]];
			while ( inserted && val_int < ecol && val_int <= L) {
				if ( val_int >= bcol ) {
					int tmp;
					tmp = val_int - bcol;
					inserted = vector_int_insert_t_partial(ab_vpos, tmp, p_elemc);
				}
				val_int = etree[val_int];
			}
		}

		//Sorting vpos
		vector_int_qsorti_partial(ab_vpos, p_elemc);

		int ccol = vptr[L];
		for ( int i = p_elemc; i < ab_vpos->elemc; ++i ) {
			while ( (vpos[ccol] - bcol) < ab_vpos->elem[i] && ccol < vptr[L+1])
				++ccol;

			if ( (vpos[ccol] - bcol) == ab_vpos->elem[i] ){
				val_fp = vval[ccol];
				vector_double_insert(ab_vval, val_fp);
				continue;
			}
			else {
				val_fp = 0;
				vector_double_insert(ab_vval, val_fp);
			}
		}
	}

	val_int = ab_vpos->elemc;
	vector_int_insert(ab_vptr, val_int);

	if ( ab_vpos->elemc ){
		block->m = b; block->n = b; block->elemc = ab_vpos->elemc;
		block->vdiag = NULL;
	}else
		printf("Warning! task2 fail. I %d J %d\n", I, J);
}



void hb2hbh_csr_task(int I, int J, hbmat_t *A, int b, int *entry, hbmat_t *block){
	int m = A->m; int n = A->n;
	int* vptr = A->vptr; int* vpos = A->vpos; 
	double* vval = A->vval;
	int brow = I*b; int erow = (I+1)*b;
	int bcol = J*b; int ecol = (J+1)*b;

	vector_t* ab_vptr = vector_create();
	vector_t* ab_vpos = vector_create();
	vector_t* ab_vval = vector_create();
	vector_clear(ab_vptr); vector_clear(ab_vpos); vector_clear(ab_vval);
	vel_t pos_val;

	*entry = 0;

	for ( int L = brow; L < erow; ++L ) {
		pos_val.i = ab_vpos->elemc; 
		vector_insert(ab_vptr, pos_val);
		
		//Padding
		if ( L >= m ) {
			if ( ecol < n )
				continue;
			pos_val.i = L - bcol;
			vector_insert(ab_vpos, pos_val);
			pos_val.d = 1;
			vector_insert(ab_vval, pos_val);
			continue;
		}

		int p_elemc = ab_vpos->elemc;
		for ( int k = vptr[L]; k < vptr[L+1]; ++k ) {
			if ( vpos[k] >= bcol && vpos[k] < ecol){
				pos_val.i = vpos[k] - bcol;
				vector_insert(ab_vpos, pos_val);
				pos_val.d = vval[k];
				vector_insert(ab_vval, pos_val);
			}
		}
	}

	pos_val.i = ab_vpos->elemc;
	vector_insert(ab_vptr, pos_val);

	if ( ab_vpos->elemc ){
		block->m = b; block->n = b; block->elemc = ab_vpos->elemc;
		block->vdiag = NULL;
		block->vptr = vector2int(ab_vptr); 
		block->vpos = vector2int(ab_vpos);
		block->vval = vector2double(ab_vval);
		*entry = 1;
	}
}

void hb2hbh_csc_task(int I, int J, hbmat_t *A, int b, int *entry, hbmat_t *block){
	int m = A->m ;int n = A->n;
	int* vptr = A->vptr; int* vpos = A->vpos; 
	double* vval = A->vval;
	int brow = I*b; int erow = (I+1)*b;
	int bcol = J*b; int ecol = (J+1)*b;

	vector_t* ab_vptr = vector_create();
	vector_t* ab_vpos = vector_create();
	vector_t* ab_vval = vector_create();
	vector_clear(ab_vptr); vector_clear(ab_vpos); vector_clear(ab_vval);
	vel_t pos_val;

	*entry = 0;

	for ( int L = bcol; L < ecol; ++L ) {
		pos_val.i = ab_vpos->elemc; 
		vector_insert(ab_vptr, pos_val);
		
		//Padding
		if ( L >= n ) {
			if ( erow < m )
				continue;
			pos_val.i = L - brow;
			vector_insert(ab_vpos, pos_val);
			pos_val.d = 1;
			vector_insert(ab_vval, pos_val);
			continue;
		}

		int p_elemc = ab_vpos->elemc;
		for ( int k = vptr[L]; k < vptr[L+1]; ++k ) {
			if ( vpos[k] >= brow && vpos[k] < erow){
				pos_val.i = vpos[k] - brow;
				vector_insert(ab_vpos, pos_val);
				pos_val.d = vval[k];
				vector_insert(ab_vval, pos_val);
			}
		}
	}

	pos_val.i = ab_vpos->elemc;
	vector_insert(ab_vptr, pos_val);

	if ( ab_vpos->elemc ){
		block->m = b; block->n = b; block->elemc = ab_vpos->elemc;
		block->vdiag = NULL;
		block->vptr = vector2int(ab_vptr); 
		block->vpos = vector2int(ab_vpos);
		block->vval = vector2double(ab_vval);
		*entry = 1;
	}
}

hbmat_t* hb2hbh(hbmat_t *A, int b, int is_csr){

	int m = A->m; int n = A->n; int elemc = A->elemc;
	int *vptr = A->vptr; int *vpos = A->vpos; double* vval = A->vval;
	int M = ( m + b - 1 ) / b;
	int N = ( n + b - 1 ) / b;
	int num = M * N;

	hbmat_t* hyper = malloc(sizeof(hbmat_t));
	hyper->m = M; hyper->n = N; hyper->vdiag = NULL;
	hyper->vval = malloc(num * sizeof(hbmat_t*));
	hbmat_t** hbmat_array = malloc(num * sizeof(hbmat_t*));
	int* hentry = malloc(num * sizeof(int));

	vector_t* ab_vptr = vector_create(); 
	vector_t* ab_vpos = vector_create();
	vector_clear(ab_vptr); vector_clear(ab_vpos);
	vel_t pos_val;

	if ( M==0 || N==0 ) {
		fprintf( stderr, "block size %i too large\n", b);
	}

	for(int i = 0; i < num; ++i)
		hbmat_array[i] = malloc(sizeof(hbmat_t));

	int acc = 0;
	int I, J;
	if (is_csr){
		for ( I = 0; I < M; ++I ) {
			for ( J = 0; J < N; ++J ) {
				hb2hbh_csr_task(I, J, A, b, &(hentry[acc]), hbmat_array[acc]);
				++acc;
			}
		}
	}else{
		for ( J = 0; J < N; ++J ) {
			for ( I = 0; I < M; ++I ) {
				hb2hbh_csc_task(I, J, A, b, &(hentry[acc]), hbmat_array[acc]);
				++acc;
			}
		}
	}

#pragma omp taskwait

	acc = 0;
	int acc0 = 0;
	if ( is_csr ){
		for ( I = 0; I < M; ++I ) {
			pos_val.i = ab_vpos->elemc;
			vector_insert(ab_vptr, pos_val);
			for ( J = 0; J < N; ++J ) {
				if ( hentry[acc] ) {
					pos_val.i = J;
					vector_insert(ab_vpos, pos_val);
					((hbmat_t**)hyper->vval)[acc0] = hbmat_array[acc];
					++acc;
					++acc0;
				} else {
					free(hbmat_array[acc]);
					++acc;
				}
			}
		}
	} else {
		for ( J = 0; J < N; ++J ) {
			pos_val.i = ab_vpos->elemc;
			vector_insert(ab_vptr, pos_val);
			for ( I = 0; I < M; ++I ) {
				if ( hentry[acc] ) {
					pos_val.i = I;
					vector_insert(ab_vpos, pos_val);
					((hbmat_t**)hyper->vval)[acc0] = hbmat_array[acc];
					++acc;
					++acc0;
				} else {
					free(hbmat_array[acc]);
					++acc;
				}
			}
		}
	}

	pos_val.i = ab_vpos->elemc;
	vector_insert(ab_vptr, pos_val);
	hyper->elemc = ab_vpos->elemc;
//	vector_printi(ab_vptr);
//	vector_printi(ab_vpos);
	hyper->vptr = vector2int(ab_vptr);
	hyper->vpos = vector2int(ab_vpos);

	set_diag(hyper);

	return hyper;
}

void set_diag(hbmat_t *A){
	int m = A->m; int n = A->n;
	int *vdiag = calloc(n, sizeof(int));
	A->vdiag = vdiag;
	int *vptr = A->vptr; int *vpos = A->vpos;
	for (int j = 0; j < n; ++j) {
		for (int k = vptr[j]; k < vptr[j+1]; ++k) {
			if (vpos[k] == j)
				vdiag[j] = k;
		}
	}
}



hbmat_t* hb2hbh_sym_etree(hbmat_t *A, int b, int* etree){
	int m = A->m; int n = A->n; int elemc = A->elemc;
	int *vptr = A->vptr; int *vpos = A->vpos; double* vval = A->vval;
	int M = ( m + b - 1 ) / b;
	int N = ( n + b - 1 ) / b;

	hbmat_t* Ab = (hbmat_t*)malloc(sizeof(hbmat_t));
	Ab->m = M; Ab->n = N; Ab->elemc = 0;
	Ab->vdiag = NULL;

	int num = ((1 + M) * N) / 2;  //maximum total number of blocks in the lower triangular matrix
	hbmat_t* acchb = (hbmat_t*) malloc(num*sizeof(hbmat_t));
	int acc = 0 ;
	Ab->vval = (hbmat_t**) malloc(num*sizeof(hbmat_t*));
	int ab_count = 0 ;
	vector_t *ab_vptr, *ab_vpos;
	ab_vptr = vector_create(); ab_vpos = vector_create();
	vector_clear(ab_vptr); vector_clear(ab_vpos);

	vector_t** vec_col = (vector_t**) malloc(n*sizeof(vector_t*));
	vel_t ab_vel;

	if ( M==0 || N==0 ) {
		fprintf( stderr, "block size %i too large\n", b);
	}

	for(int i = 0; i < m; i++){
		vec_col[i] = vector_create();
		vector_clear(vec_col[i]);
	}

	for(int J = 0; J < N; J++){
		int jstart = J * b;
		int jc = n - jstart;
		jc = jc < b ? jc: b; 
		ab_vel.i = ab_vpos->elemc;
		vector_insert(ab_vptr, ab_vel);

		for(int j = jstart; j < jstart+jc; j++){
			for(int k = vptr[j]; k < vptr[j+1]; k++){
				ab_vel.i = vpos[k];
				vector_insert_t(vec_col[j], ab_vel);
			}
			vector_qsorti(vec_col[j]);
			for(int l = 2; l < vec_col[j]->elemc; l++){
				ab_vel.i = vec_col[j]->elem[l].i;
				vector_insert_t(vec_col[vec_col[j]->elem[1].i], ab_vel);
			}
		}

		for(int I = J; I < M; I++){
			int base_col = J * b ;
			int base_row = I * b ;
			int ic = m - base_row;
			int max_row = base_row+b;

			ic = ic < b ? ic : b;
			vector_t *sub_col, *sub_vptr, *sub_vpos, *sub_vval;
			vel_t acchb_vel;
			sub_col = vector_create_size(jc);
			sub_vptr = vector_create_size(jc);
			sub_vpos = vector_create_size(jc);
			sub_vval = vector_create_size(jc);
			vector_clear(sub_vptr);
			vector_clear(sub_vpos);
			vector_clear(sub_vval);

			for(int j = 0; j < b; j++){
				int current_col = base_col+j; //Absolute column position
				vector_clear(sub_col);
				ab_vel.i = sub_vpos->elemc;
				vector_insert(sub_vptr, ab_vel);

				//Column Padding
				if (current_col >= n){
					acchb_vel.i = current_col-base_row;
					vector_insert(sub_vpos, acchb_vel);
					acchb_vel.d = 1;
					vector_insert(sub_vval, acchb_vel);
					continue;
				}

				for(int k = 0; k < vec_col[current_col]->elemc; k++){
					if (vec_col[current_col]->elem[k].i < base_row)
						continue;
					if (vec_col[current_col]->elem[k].i >= max_row)
						break;
					acchb_vel.i = vec_col[current_col]->elem[k].i;
					vector_insert(sub_col, acchb_vel);
				}

				for(int k = 0; k < sub_col->elemc; k++){
					acchb_vel.i = sub_col->elem[k].i-base_row;
					vector_insert(sub_vpos, acchb_vel);
					for(int l = vptr[current_col]; l < vptr[current_col+1]; l++){
						acchb_vel.d = 0;
						if(sub_col->elem[k].i == vpos[l]){
							acchb_vel.d = vval[l];
							break;
						}
					}
					vector_insert(sub_vval, acchb_vel);
				}
			}
			if (sub_vpos->elemc != 0){
				ab_vel.i = sub_vpos->elemc;
				vector_insert(sub_vptr, ab_vel);
//				vector_printi(sub_vptr);
//				vector_printi(sub_vpos);
//				vector_printd(sub_vval);
				acchb[ab_vpos->elemc].m = b;
				acchb[ab_vpos->elemc].n = b;
				acchb[ab_vpos->elemc].elemc = sub_vpos->elemc;
				acchb[ab_vpos->elemc].vptr = vector2int(sub_vptr);
				acchb[ab_vpos->elemc].vpos = vector2int(sub_vpos);
				acchb[ab_vpos->elemc].vval = vector2double(sub_vval);
				acchb[ab_vpos->elemc].vdiag = NULL;
				((hbmat_t**)Ab->vval)[ab_vpos->elemc] = acchb + ab_vpos->elemc;
				ab_vel.i = I;
				vector_insert(ab_vpos, ab_vel);
			}
			else{
				vector_free(sub_vptr);
				vector_free(sub_vpos);
				vector_free(sub_vval);
			}
			vector_free(sub_col);
		}
	}

	for(int i = 0; i < n; i++)
		vector_free(vec_col[i]);
	free(vec_col);

	ab_vel.i = ab_vpos->elemc;
	vector_insert(ab_vptr, ab_vel);
	Ab->elemc = ab_vpos->elemc;
	Ab->vptr = vector2int(ab_vptr);
	Ab->vpos = vector2int(ab_vpos);

	return Ab;

}

hbmat_t *hbh2hb (hbmat_t *A){
	hbmat_t *B = (hbmat_t*) malloc(sizeof(hbmat_t));
	//Assuming the input matrix A is lower triangular
	int M = A->m; int N = A->n;
	int elemc = A->elemc;
	int* vptr = A->vptr;
	int* vpos = A->vpos;
	hbmat_t** vval = A->vval;
	
	vector_t *b_vptr, *b_vpos, *b_vval;
	b_vptr = vector_create(); b_vpos = vector_create(); b_vval = vector_create();
	vector_clear(b_vptr); vector_clear(b_vpos); vector_clear(b_vval);
	vel_t b_vptr_vel, b_vpos_vel, b_vval_vel;
	hbmat_t* sub_matrix;
	int bs = vval[0]->m; //Block size can be determined by the rows of the first sub-matrix
	int col_counter = 0;

	for(int J = 0; J < N; J++){
		sub_matrix = vval[vptr[J]]; //Fetch the first sub-matrix in this column
		int tot_col = sub_matrix->n;
		for(int j = 0; j < tot_col; j++){
			col_counter++;
			b_vptr_vel.i = b_vpos->elemc;
			vector_insert(b_vptr, b_vptr_vel);
			for(int I = vptr[J]; I < vptr[J+1]; I++){
				int c_row = vpos[I];
				int row_offset = c_row*bs;
				sub_matrix = vval[I];
				for(int jj = sub_matrix->vptr[j]; jj < sub_matrix->vptr[j+1]; jj++){
					if(1 || ((double*)sub_matrix->vval)[jj] != 0 ){
						b_vpos_vel.i = sub_matrix->vpos[jj] + row_offset;
						vector_insert(b_vpos, b_vpos_vel);
						b_vval_vel.d = ((double*)sub_matrix->vval)[jj];
						vector_insert(b_vval, b_vval_vel);
					}
				}
			}
		}
	}

	b_vptr_vel.i = b_vpos->elemc;
	vector_insert(b_vptr, b_vptr_vel);
	
//	vector_printi(b_vptr);vector_printi(b_vpos);vector_printd(b_vval);

	B->m = B->n = col_counter;
	B->elemc = b_vpos->elemc;
	B->vptr = vector2int(b_vptr);
	B->vpos = vector2int(b_vpos);
	B->vval = vector2double(b_vval);
	return B;
}

hbmat_t *hbh2hb_sym (hbmat_t *A){
	hbmat_t *B = (hbmat_t*) malloc(sizeof(hbmat_t));
	//Assuming the input matrix A is lower triangular
	int M = A->m; int N = A->n;
	int elemc = A->elemc;
	int* vptr = A->vptr;
	int* vpos = A->vpos;
	hbmat_t** vval = A->vval;
	
	vector_t *b_vptr, *b_vpos, *b_vval;
	b_vptr = vector_create(); b_vpos = vector_create(); b_vval = vector_create();
	vector_clear(b_vptr); vector_clear(b_vpos); vector_clear(b_vval);
	vel_t b_vptr_vel, b_vpos_vel, b_vval_vel;
	hbmat_t* sub_matrix;
	int bs = vval[0]->m; //Block size can be determined by the rows of the first sub-matrix
	int col_counter = 0;

	for(int J = 0; J < N; J++){
		sub_matrix = vval[vptr[J]]; //Fetch the first sub-matrix in this column
		int tot_col = sub_matrix->n;
		for(int j = 0; j < tot_col; j++){
			col_counter++;
			b_vptr_vel.i = b_vpos->elemc;
			vector_insert(b_vptr, b_vptr_vel);
			for(int I = vptr[J]; I < vptr[J+1]; I++){
				int c_row = vpos[I];
				int row_offset = c_row*bs;
				sub_matrix = vval[I];
				for(int jj = sub_matrix->vptr[j]; jj < sub_matrix->vptr[j+1]; jj++){
					if(1 || ((double*)sub_matrix->vval)[jj] != 0 ){
						b_vpos_vel.i = sub_matrix->vpos[jj] + row_offset;
						vector_insert(b_vpos, b_vpos_vel);
						b_vval_vel.d = ((double*)sub_matrix->vval)[jj];
						vector_insert(b_vval, b_vval_vel);
					}
				}
			}
		}
	}

	b_vptr_vel.i = b_vpos->elemc;
	vector_insert(b_vptr, b_vptr_vel);
	
//	vector_printi(b_vptr);vector_printi(b_vpos);vector_printd(b_vval);

	B->m = B->n = col_counter;
	B->elemc = b_vpos->elemc;
	B->vptr = vector2int(b_vptr);
	B->vpos = vector2int(b_vpos);
	B->vval = vector2double(b_vval);
	return B;
}

hbmat_t *hb2csr(hbmat_t *A){
	int m = A->m; int n = A->n; int elemc = A->elemc;
	int* vptr = A->vptr; int* vpos = A->vpos;
	double* vval = (double*)A->vval;
	hbmat_t* B = (hbmat_t*)malloc(sizeof(hbmat_t));
	vector_t *vptr_b, **vpos_b, **vval_b;
	vptr_b = vector_create();
	vpos_b = malloc(m*sizeof(vector_t*));
	vval_b = malloc(m*sizeof(vector_t*));
	vector_clear(vptr_b);
	vel_t b_vel;

	for(int i = 0; i < m; i++){
		vpos_b[i] = vector_create();
		vector_clear(vpos_b[i]);
		vval_b[i] = vector_create();
		vector_clear(vval_b[i]);
	}

	for (int j = 0; j < n; j++){
		for(int i = vptr[j]; i < vptr[j+1]; i++){
			b_vel.i = j;
			vector_insert(vpos_b[vpos[i]], b_vel);
			b_vel.d = vval[i];
			vector_insert(vval_b[vpos[i]], b_vel);
		}
	}

	//Update vptr_b
	b_vel.i = 0;
	vector_insert(vptr_b, b_vel);
	for(int i = 0; i < m; i++){
		b_vel = vector_get(vptr_b,i);
		b_vel.i += vpos_b[i]->elemc;
		vector_insert(vptr_b, b_vel);
	}

	for(int i = 1; i < m; i++){
		vpos_b[0] = vector_append(vpos_b[0],vpos_b[i]->elem, vpos_b[i]->elemc);
		vector_free(vpos_b[i]);
		vval_b[0] = vector_append(vval_b[0],vval_b[i]->elem, vval_b[i]->elemc);
		vector_free(vval_b[i]);
	}

	B->m = m; B->n = n; B->elemc = elemc;
	B->vptr = vector2int(vptr_b);
	B->vpos = vector2int(vpos_b[0]);
	B->vval = vector2double(vval_b[0]);
	free(vpos_b); free(vval_b);

	return B;
}

void hb_free(hbmat_t *A){
	free(A->vptr); free(A->vpos);
	free(A->vval);
	free(A);
}

void hbh_free(hbmat_t *A){
	int elemc = A->elemc;
	for(int i = 0; i < elemc; i++){
		free(((hbmat_t**)A->vval)[i]->vptr);
		free(((hbmat_t**)A->vval)[i]->vpos);
		free(((hbmat_t**)A->vval)[i]->vval);
	}
	free(((hbmat_t**)A->vval)[0]);
	hb_free(A);
}

void hbh_free2(hbmat_t *A){

	pthread_mutex_destroy(A->mtx);
	int elemc = A->elemc;
	free(A->e_tree);
	free(((hbmat_t**)A->vval)[0]);
	free(A->vptr_pool);
	free(A->vpos_pool);
	free(A->vval_pool);
	hb_free(A);

}

hbmat_t *hb2hbb(hbmat_t *A, int b) {
	//printf("hb2hbb %p %i\n", A, b);
	int m = A->m;
	int n = A->n;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	double *vval = A->vval;

	if ( vval == NULL ) {
		fprintf( stderr, "error: hb2hbb misses values\n");
		return NULL;
	}

	int M = ( m + b - 1 ) / b;
	int N = ( n + b - 1 ) / b;

	hbmat_t *Ab = (hbmat_t*) malloc( sizeof(hbmat_t) );
	Ab->m = M; Ab->n = N; Ab->b = b; Ab->elemc = A->elemc;
	Ab->vdiag = NULL;
	Ab->vptr = (int*) malloc( sizeof(int) * ( N + 1 ) );
	vector_t *accvpos = vector_create();
	double **accb = (double**) malloc( sizeof(double*) * A->elemc);

	vector_t *rfront = vector_create();
	vector_t *relemc = vector_create();
	int acc = 0;

	if ( M==0 || N==0 ) {
		fprintf( stderr, "block size %i too large\n", b);
	}


	int J;
	for ( J = 0; J < N; J++ ) {
		vector_clear(rfront);
		vector_clear(relemc);
		
		int jstart = J * b;
		int jc = n - jstart;
		jc = jc < b? jc: b;
		int j;
		for ( j = 0; j < jc; j++ ) {
			vel_t vel = { i: vptr[J*b+j] - 1};
			vector_insertat( rfront, vel, j );
			vel.i = 0;
			vector_insertat( relemc, vel, j );
		}

		Ab->vptr[J] = acc+1;
			
		int border = b;
		int I;	
		for ( I = 0; I < M; I++ ) {
			//printf("block (%i,%i) %i\n", I, J, relemc->elemc);
			double *cB = NULL;  

			int col = J * b;
			for ( j = 0 ; j < jc; j++ ) {
				int elemc = vector_get(relemc, j).i;
				int cnt = vptr[col+1] - vptr[col];
				int r = vpos[ vector_get( rfront, j ).i ];
				//printf("\tcol %i (%i) border %i initr %i\n", j, cnt, border, r);

				while ( r <= border && elemc < cnt) {
					if ( cB == NULL ) {
						cB = (double*) calloc( b * b, sizeof(double));
					}
			
					int frontier = vector_get( rfront, j ).i;
					cB[j*b + ((r - 1)%b)] = vval[frontier];
					//printf("(%i,%i) = %f front %i\n", (r-1)%b, j, vval[frontier], frontier);

					vel_t vel = { i: ++frontier};
					vector_insertat( rfront, vel, j);
					r = vpos[frontier];

					++elemc;
				}

				vel_t vel = { i: elemc };
				vector_insertat( relemc, vel, j);
				//vector_print( relemc );
				++col;	
			}

			if ( cB != NULL ) {
				vel_t vel = { i: I+1 };
				vector_insert( accvpos, vel);
				accb[acc++] = cB;
			}

			border += b;
		}
	}

	Ab->vptr[J] = acc+1;

	Ab->vval = (double**) realloc( accb, sizeof(double*) * acc );
	Ab->vpos = (int*) vector2int(accvpos);
	Ab->elemc = acc;

	vector_free( rfront );
	vector_free( relemc );

	return Ab;
}

	
hbmat_t* hbb2csrb(hbmat_t *A) {
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	double **vval = A->vval;
	int m = A->m;
	int n = A->n;
	int b = A->b;
	int elemc = A->elemc;

	if ( b==0 ) {
		fprintf(stderr, "warning: hbb2csrb: A has block size 0\n");
	}

	hbmat_t *Acsr = (hbmat_t*) malloc( sizeof(hbmat_t) );
	Acsr->vval = (void**) malloc( sizeof(double*) * elemc );
	Acsr->vpos = (int*) malloc( sizeof(int) * elemc );
	Acsr->vptr = (int*) malloc( sizeof(int) * (m + 1)) ;
	Acsr->elemc = elemc; Acsr->m = m; Acsr->n = n;
	Acsr->b = b;
	Acsr->vdiag = NULL;
	void** csr_vval = Acsr->vval;

	vector_t *cfront = vector_create();
	vector_t *celemc = vector_create();
	int j;
	for ( j = 0; j < n; j++ ) {
		vel_t vel = { i: vptr[j] - 1 };
		vector_insert( cfront, vel );
		vel.i = 0;
		vector_insert( celemc, vel );
	}
		

	int c = 0;
	int i;
	for ( i = 0; i < m; i++ ) {
		Acsr->vptr[i] = c + 1;
	
		int j;
		for ( j = 0; j < n ; j++ ) {
			int currc = vector_get( celemc, j ).i;
			int colc = vptr[j+1] - vptr[j];

			//printf("looking at col %i : %i %i\n", j, currc, colc);

			if ( currc < colc ) {
				int start = vector_get( cfront, j ).i;
				int row = vpos[start];

				if ( i == row - 1 ) {
// copies pointers to blocks, not blocks
					double *block = vval[start];

					Acsr->vpos[c] = j + 1;
					csr_vval[c++] = block;
	
					vel_t vel = { i: start+1 };
					vector_insertat ( cfront, vel, j );
					vel.i = currc+1;
					vector_insertat ( celemc, vel, j );
				}
			}
		}
	} 
	Acsr->vptr[i] = c + 1;

	vector_free(cfront);
	vector_free(celemc);

	return Acsr;
}

#if 0
hbmatm_t* b2ll(hbmat_t *A) {
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	double **vval = A->vval;
	int *vdiag = A->vdiag;
	int m = A->m;
	int n = A->n;
	int b = A->b;
	int elemc = A->elemc;

	if ( b==0 ) {
		fprintf(stderr, "warning: b2ll: A has block size 0\n");
	}

	hbmatm_t *Am = (hbmatm_t*) malloc( sizeof(hbmatm_t) );
	if ( Am == NULL ) {
		fprintf(stderr, "error: b2ll: allocation failed (1)\n");
		return NULL;
	}
	dll3_t *mvval = Am->vval = dll3_create();
	if ( mvval == NULL ) {
		free( Am );
		fprintf(stderr, "error: b2ll: allocation failed (2)\n");
		return NULL;
	}
	int *mvptr = Am->vptr = (int*) malloc( sizeof(int) * (m + 1) );
	if ( mvptr == NULL ) {
		free( Am );
		dll3_destroy(mvval);
		fprintf(stderr, "error: b2ll: allocation failed (3)\n");
		return NULL;
	}
	//dll_t *mvpos = Am->vpos = dll_create();
	int *mvdiag = Am->vdiag = (int*) malloc( sizeof(int) * m );
	if ( mvdiag == NULL ) {
		free( Am );
		dll3_destroy( mvval );
		free( mvptr );
		fprintf(stderr, "error: b2ll: allocation failed (4)\n");
		return NULL;
	}
	dll3_t **mvdiagl = Am->vdiagl = (dll3_t**) malloc( sizeof(dll3_t*) * m );
	dll3_t **mvptrv = Am->vptrv = (dll3_t**) malloc( sizeof(dll3_t*) * (m+1) );
	Am->elemc = elemc; Am->m = m; Am->n = n;
	Am->b = b;

	int c;
	for ( c = 0 ; c < elemc; c++ ) {
		//llel3_t e3; 
		//e3.A = vval[c];
		dll3_insert( mvval, vval[c], NULL, NULL, vpos[c], mvval);

#if 0
		dllel_t e;
		e.i = vpos[c];
		dll_insert( mvpos, e, mvpos);
#endif
	}


	for ( c = 0 ; c <= m; c++ ) {
		int idx = vptr[c];
		mvptr[c] = idx;
	}

	for ( c = 0; c < m; c++ ) {
		mvptrv[c] = dll3_get( mvval, mvptr[c] - 1 );
	}
	mvptrv[m] = mvval;

	if ( vdiag != NULL ) {
		int c;
		for ( c = 0 ; c < m; c++ ) {
			int d = vdiag[c];
			mvdiag[c] = d;
			mvdiagl[c] = dll3_get( mvval, d-1 );
		}
	}

	return Am;
}


hbmat_t* ll2b(hbmatm_t *A) {
	int *vptr = A->vptr;
	//dll_t *vpos = A->vpos;
	dll3_t *vval = A->vval;
	int *vdiag = A->vdiag;
	int m = A->m;
	int n = A->n;
	int b = A->b;
	int elemc = A->elemc;

	if ( b==0 ) {
		fprintf(stderr, "warning: ll2b: A has block size 0\n");
	}

	hbmat_t *Ab = (hbmat_t*) malloc( sizeof(hbmat_t) );
	int *bvptr = Ab->vptr = (int*) malloc( sizeof(int) * (m + 1) );
	int *bvpos = Ab->vpos = (int*) malloc( sizeof(int) * elemc );
	double **bvval = Ab->vval = (double**) malloc( sizeof(double*) * elemc );
	Ab->elemc = elemc; Ab->m = m; Ab->n = n;
	Ab->b = b;

	//dll3_t *pos = vpos->next;
	dll3_t *val = vval->next;
	int c;
	for ( c = 0 ; c < elemc; c++ ) {
		bvpos[c] = val->e.i;
		bvval[c] = val->e.A;
		
		//pos = pos->next;
		val = val->next;
	}

	for ( c = 0 ; c <= m; c++ ) {
		bvptr[c] = vptr[c]; 
	}


	if ( vdiag != NULL ) {
		int *bvdiag = Ab->vdiag = (int*) malloc( sizeof(int) * m );

		int c;
		for ( c = 0 ; c < m; c++ ) {
			bvdiag[c] = vdiag[c];
		}
	} else {
		Ab->vdiag = NULL;
	}

	return Ab;
}


hbmatm_t* hbbm2csrbm(hbmatm_t *A) {
	int *vptr = A->vptr;
	//dll_t *vpos = A->vpos;
	dll3_t *vval = A->vval;
	int M = A->m;
	int N = A->n;
	int b = A->b;
	int elemc = A->elemc;

	if ( b==0 ) {
		fprintf(stderr, "warning: hbbm2csrbm: A has block size 0\n");
	}

	hbmatm_t *Acsr = (hbmatm_t*) malloc( sizeof(hbmatm_t) );
	if ( Acsr == NULL ) {
		fprintf(stderr, "hbbm2csrbm: cannot allocate Acsr\n");
		return NULL;
	}
	dll3_t *csr_vval = Acsr->vval = dll3_create();
	if ( csr_vval == NULL ) {
		free( Acsr );
		fprintf(stderr, "hbbm2csrbm: cannot allocate csr_vval\n");
		return NULL;
	}
	//dll_t *csr_vpos = Acsr->vpos = dll_create();
	int *csr_vptr = Acsr->vptr = (int*) malloc( sizeof(int) * (M + 1)) ;
	if ( csr_vptr == NULL ) {
		free( Acsr );
		dll3_destroy( csr_vval );
		fprintf(stderr, "hbbm2csrbm: cannot allocate csr_vptr\n");
		return NULL;
	}
	Acsr->elemc = elemc; Acsr->m = M; Acsr->n = N;
	Acsr->b = b;
	Acsr->vdiag = NULL;
	Acsr->vdiagl = A->vdiagl;
	dll3_t **csr_vptrv = Acsr->vptrv = (dll3_t**) malloc( sizeof(dll3_t*) * (M+1));
	if ( csr_vptrv == NULL ) {
		free( Acsr );
		dll3_destroy( csr_vval );
		free(csr_vptr);
		fprintf(stderr, "hbbm2csrbm: cannot allocate csr_vptrv\n");
		return NULL;
	}
	Acsr->elemc = elemc; Acsr->m = M; Acsr->n = N;
	Acsr->b = b;
	

	vector_t *cfront = vector_create();
	vector_t *celemc = vector_create();
	int j;
	for ( j = 0; j < N; j++ ) {
		vel_t vel = { i: vptr[j] - 1 };
		vector_insert( cfront, vel );
		vel.i = 0;
		vector_insert( celemc, vel );
	}
		

	int c = 0;
	int i;
	for ( i = 0; i < M; i++ ) {
		csr_vptr[i] = c + 1;
		//printf("setting vptr[%i] to %i\n", i, c+1);

		dll3_t *first = NULL;
	
		int j;
		for ( j = 0; j < N ; j++ ) {
			int currc = vector_get( celemc, j ).i;
			int colc = vptr[j+1] - vptr[j];

			//printf("looking at col %i : %i %i\n", j, currc, colc);

			if ( currc < colc ) {
				int start = vector_get( cfront, j ).i;
				dll3_t *blocks = dll3_get( vval, start ); 
				int row = blocks->e.i;

				//printf("\tfound row %i\n", row );

				if ( i == row - 1 ) {
// copies pointers to blocks, not blocks
					//dll3_t *blocks = dll3_get( vval, start);
					
					//dllel_t llel = { i: j+1 };
					//dll_insert( csr_vpos, llel, csr_vpos );

					//llel3_t llel3;
					//printf("found %p\n", blocks->e.A);
					dll3_insert( csr_vval, blocks->e.A, blocks->e.S, blocks->e.P, j+1, csr_vval );
					++c;
	
					vel_t vel = { i: start+1 };
					vector_insertat ( cfront, vel, j );
					vel.i = currc+1;
					vector_insertat ( celemc, vel, j );

					if ( first == NULL ) {
						first = csr_vval->prev;
					}
				}
			}
		}

		csr_vptrv[i] = first;
	} 
	csr_vptr[i] = c + 1;
	csr_vptrv[i] = csr_vval;

	vector_free(cfront);
	vector_free(celemc);

	return Acsr;
}


void llsetdiag( hbmatm_t *A ){
	int m = A->m;
	printf("llsetdiag %p %i\n", A, m);

	if ( A->vdiagl != NULL ) {
		fprintf(stderr, "warning: llsetdiag: vdiagl not NULL\n");
		return; 
	}
	int *vdiag = A->vdiag;
	if ( vdiag == NULL ) {
		fprintf(stderr, "error: llsetdiag: vdiag NULL\n");
	}

	dll3_t **vdiagl = A->vdiagl = (dll3_t**) malloc( sizeof(dll3_t*) * m );
	dll3_t *vval = A->vval->next;

	int seen = 0;
	int c;
	for ( c=0; c<m; c++ ) {
		int idx = vdiag[c] - 1;
		int skip = idx - seen;
		seen = idx;
	
		vval = dll3_scan(vval, skip);
		vdiagl[c] = vval;
	}
}



void m_sync(hbmatm_t *Ahb, hbmatm_t *Acsr) {
	int *vptr = Ahb->vptr;
	//dll_t *vpos = Ahb->vpos;
	dll3_t *vval = Ahb->vval;
	int m = Ahb->m;
	int n = Ahb->n;
	int b = Ahb->b;

	if ( b==0 ) {
		fprintf(stderr, "warning: m_sync: A has block size 0\n");
	}

	int *cvptr = Acsr->vptr;
	//int *cvpos = Acsr->vpos;
	dll3_t *cvval = Acsr->vval;

	vector_t *cfront = vector_create();
	vector_t *celemc = vector_create();
	int j;
	for ( j = 0; j < n; j++ ) {
		vel_t vel = { i: vptr[j] - 1 };
		vector_insert( cfront, vel );
		vel.i = 0;
		vector_insert( celemc, vel );
	}
		

	int i;
	for ( i = 0; i < m; i++ ) {
		int csri = cvptr[i] - 1;
	
		int j;
		for ( j = 0; j < n ; j++ ) {
			int currc = vector_get( celemc, j ).i;
			int colc = vptr[j+1] - vptr[j];

			if ( currc < colc ) {
				int start = vector_get( cfront, j ).i;
				dll3_t *blocks = dll3_get( vval, start );
				int row = blocks->e.i;

				if ( i == row - 1 ) {
					// copies pointers to blocks, not blocks
					//dll3_t *blocks = dll3_get( vval, start );

					dll3_t *val = dll3_get( cvval, csri );
					while ( val->e.A !=  blocks->e.A ) {
						val = val->next;
					}

					val->e.S = blocks->e.S;
					val->e.P = blocks->e.P;

					if ( val->e.S == NULL || val->e.P == NULL ) {
						fprintf( stderr, "warning: m_sync: insert NULL\n");
					}

					vel_t vel = { i: start+1 };
					vector_insertat ( cfront, vel, j );
					vel.i = currc+1;
					vector_insertat ( celemc, vel, j );
				}
			}
		}
	} 

	vector_free(cfront);
	vector_free(celemc);
}
#endif
