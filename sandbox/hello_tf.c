#include <stdio.h>
#include <tensorflow/c/c_api.h>
#include <random>

int main() {
	  printf("Hello from TensorFlow C library version %s\n", TF_Version());
	    return 0;
}
