#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_coupling.h>

// Convert float spin to 2j format (integers)

#define TWICE(x) ((int)(2 * (x)))

int main(void){
    // Spins
    double j = 0.5;
    double m = 0.5;
    int two_j = TWICE(j);
    int two_m = TWICE(m);

    // Operator: vector operator T^{(1)}_0 corresponds to q=0, rank k = 1
    
    int k = 1;
    int q = 0;

    //Final state same: j', m'

    int two_jp = two_j;
    int two_mp = two_m;

    // Compute Clebsch-Gordan coefficient using 3j symbol:
    // C^{j', m'}_{j,m;k,q} = (-1)^{j'-j-q} * sqrt(2j'+1)*(j j' k; m -m' -q)

    double phase = pow(-1.0, (two_jp - two_j - 2*q)/2.0);
    double cg = gsl_sf_coupling_3j(two_j, TWICE(k), two_jp, two_m, TWICE(q), -two_mp);

    printf("3j coefficient: %.6f\n", cg);

    cg *= phase * sqrt(two_jp + 1.0);

    //Reduced matrix element for spin-1/2 operator

    double red = sqrt(3.0) / 2.0;

    // Apply Wigner-Eckart theorem

    double matrix_element = red * cg;

    printf("Clebsch-Gordan coefficient: %.6f\n", cg);
    printf("Reduced matrix element: %.6f\n", red);
    printf("Full matrix element: %.6f\n", matrix_element); // should be 0.5

    return 0;

}
