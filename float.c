#include <math.h>
#include <stdio.h>
#include <float.h>

int main () {
	float a=1,b,c;
	for (int i=1; i<100; i++) {
		b = ldexp(1.0f, -i);
		c = a+b;
		if (a==c) {
			printf(" Float Error at 1+2^(-%i) => epsilon:%e, FLT_EPSILON:%e\n",
                   i,b,FLT_EPSILON);
			break;
		}
	}
	
	double d=1,e,f;
	for (int i=1; i<100; i++) {
		e = ldexp(1.0d, -i);
		f = d+e;
		if (d==f) {
			printf("Double Error at 1+2^(-%i) => epsilon:%e, DBL_EPSILON:%e\n",
                   i,e,DBL_EPSILON);
			break;
		}
	}
    
    /* Because of this internally units ARE needed! Calculating with SI the
     * whole time uses up unneeded exponent "space" */
    printf("Max Float:%e, Float Min:%e\n", FLT_MAX, FLT_MIN);
    printf("Max Float:%i, Float Min:%i\n", FLT_MAX_EXP, FLT_MIN_EXP);
    /* why the fuck is this working? => "FLT_MAX_EXP: The maximal exponent
     * of a floating point value expressed in base FLT_RADIX; greater
     * exponents are principally possible (up to 16383), but not supported
     * in all math functions." */
    float fer = 1.00e-42;
    float rec = 1.0/fer;
    printf("1.0/%e = %e\n", fer, rec);
    /* printf shows the right number only if %e instead of %f is used.
     * This is because %f always shows non scientific representation to a
     * certain amount of numbers, which need to be rounded. Meaning the
     * following will show 0, although it isn't */
    printf("0? : %f\n", 1e-10);
    
	return 0;
}

