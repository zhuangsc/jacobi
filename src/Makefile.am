SUBDIRS = matutil 

jacobi_sources = jacobi.c jacobi_setup.c jacobi_kernels.c jacobi_main.c
jacobi_headers = jacobi_setup.h jacobi_kernels.h jacobi_main.h

bin_PROGRAMS = jacobi

# Jacobi
jacobi_SOURCES = $(jacobi_sources) $(jacobi_headers)
jacobi_LDFLAGS = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) -lm
jacobi_CFLAGS = $(OMPSS_CFLAGS) -I$(srcdir)/matutil -I$(srcdir)/util
jacobi_LDADD = $(builddir)/matutil/libmatutil.la 
