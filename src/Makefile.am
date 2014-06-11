SUBDIRS = matutil

jacobi_sources = jacobi_setup.c
jacobi_headers = jacobi_setup.h

bin_PROGRAMS = jacobi

# Jacobi
jacobi_SOURCES = $(jacobi_sources) $(jacobi_headers)
jacobi_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -L/home/szhuang/mnt/Experiments/algebra/cholesky/build/lib -lompsschol
jacobi_CFLAGS = $(OMPSS_CFLAGS) -I$(srcdir)/matutil -I$(srcdir)/util -I/home/szhuang/mnt/Experiments/algebra/cholesky/build/include
jacobi_LDADD = $(builddir)/matutil/libmatutil.la
