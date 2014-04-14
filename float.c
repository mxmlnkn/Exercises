#include <math.h>
#include <stdio.h>

int main () {
	float a=1,b,c;
	for (int i=1; i<100; i++) {
		b = ldexp(1.0f, -i);
		c = a+b;
		if (a==c) {
			printf("Float Error at 1+2^(-%i)\n",i);
			break;
		}
	}
	
	double d=1,e,f;
	for (int i=1; i<100; i++) {
		e = ldexp(1.0d, -i);
		f = d+e;
		if (d==f) {
			printf("Double Error at 1+2^(-%i)\n",i);
			break;
		}
	}
	return 0;
}

//same result could be easily obtained with float.c makros: FLT_EPSILON, DBL_EPSILON.
