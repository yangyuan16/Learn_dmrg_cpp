
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_coupling.h>

// when I know the reduced matrix element <1/2 || T^1 || 1/2> = sqrt(3) / 2
// I can easily calculate the <1/2 1/2| T^1_{1} |1/2 1/2> (<j' m'| T^k_{q} |j m>)

int main(void) {
    // j = 1/2, m = -1/2; j' = 1/2, m' = +1/2; k = 1, q = +1
    int two_j = 1, two_m = -1;
    int two_jp = 1, two_mp = 1;
    int two_k = 2, two_q = 2;

    // Clebsch–Gordan coefficient via 3j:
    // C^{j', m'}_{j, m; k, q} = (-1)^{j'-j-q} * sqrt(2j'+1) * (j k j'; m q -m')
    double phase = pow(-1.0, (two_jp - two_j - two_q)/2.0);
    double three_j = gsl_sf_coupling_3j(two_j, two_k, two_jp,
                                        two_m, two_q, -two_mp);
    double cg = phase * sqrt(two_jp + 1.0) * three_j;

    printf("3j coefficient: %.6f\n", three_j);

    // Reduced matrix element <1/2 || T^1 || 1/2> = sqrt(3)/2
    double red = sqrt(3.0) / 2.0;
    double result = red * cg; //  <1/2 1/2| T^1_{1} |1/2 1/2> (<j' m'| T^k_{q} |j m>)

    printf("CG coefficient: %.6f\n", cg);
    printf("Reduced matrix element: %.6f\n", red);
    printf("Full matrix element: %.6f\n", result);

    return 0;
}
