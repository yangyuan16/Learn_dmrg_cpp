#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_coupling.h>

// Convert float spin to 2j format (integers)

#define TWICE(x) ((int)(2 * (x) + 0.1))

int main(void){
    // Spins
    double jl = 0.5;
    double ml = 0.5;
    int two_jl = TWICE(jl);
    int two_ml = TWICE(ml);

    printf("two_jl=%d\n", two_jl);
    printf("two_ml=%d\n", two_ml);

    double jr = 0.5;
    double mr = 0.5;
    int two_jr = TWICE(jr);
    int two_mr = TWICE(mr);

    printf("two_jr=%d\n", two_jr);
    printf("two_mr=%d\n", two_mr);

    // Operator: vector operator T^{(1)}_0 corresponds to q=0, rank k = 1
    
    int k = 1;
    int q = 0;

    //Final state same: k, q

    int two_k = TWICE(k);
    int two_q = TWICE(q);

    printf("two_k=%d\n", two_k);
    printf("two_q=%d\n", two_q);


    // input jl, ml, jr, mr, k, q to "gsl_sf_coupling_3j" 
    int two_j1 = two_jl;
    int two_m1 = two_ml;

    int two_j2 = two_k;
    int two_m2 = two_q;

    int two_j3 = two_jr;
    int two_m3 = two_mr;
    
    double phase = pow(-1.0, ( two_j1 - two_m1)/2.0);
    double w3j = gsl_sf_coupling_3j(two_j1, two_j2, two_j3, -two_m1, two_m2, two_m3);
    printf("w3j coefficient: %.6f\n", w3j);
    double cg = w3j * phase;

    printf("cg coefficient: %.6f\n", cg);

    return 0;

}
