
#include "jacobi_setup.h"

void *Ahb;
void *Ahbh;
int format;
int bs;

int main(int argc, char* argv[]){

	jacobi_setup(argc, argv);

	jacobi_shutdown();

	return 0;
}

