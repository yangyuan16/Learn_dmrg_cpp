
# Wigner–Eckart Theorem with GSL (GNU Scientific Library)

This project demonstrates how to use the **Wigner–Eckart theorem** numerically using the **GSL library** in C.

---

## 🎯 Goal

Compute matrix elements of spherical tensor operators using:

\[
\langle j', m' | T^{(k)}_q | j, m \rangle = \langle j' || T^{(k)} || j \rangle \cdot C^{j', m'}_{j, m; k, q}
\]

Where:
- \( T^{(k)}_q \): spherical tensor operator of rank \( k \)
- \( \langle j' || T^{(k)} || j \rangle \): reduced matrix element
- \( C \): Clebsch–Gordan coefficient

---

## ✅ GSL Functions

- `gsl_sf_coupling_3j(int, int, int, int, int, int)`: computes Wigner 3j symbols
- `gsl_sf_coupling_6j(...)`: computes Wigner 6j symbols (used in advanced coupling)

---

## 💡 Example: Compute \( \langle \tfrac{1}{2}, \tfrac{1}{2} | S_z | \tfrac{1}{2}, \tfrac{1}{2} \rangle \)

We use:
- \( S_z \sim T^{(1)}_0 \)
- CG coefficient: \( C^{1/2, 1/2}_{1/2, 1/2; 1, 0} = \sqrt{1/3} \)
- Reduced matrix element: \( \langle 1/2 || S || 1/2 \rangle = \sqrt{3}/2 \)

Then:
\[
\langle \tfrac{1}{2}, \tfrac{1}{2} | S_z | \tfrac{1}{2}, \tfrac{1}{2} \rangle = \frac{\sqrt{3}}{2} \cdot \sqrt{\frac{1}{3}} = \frac{1}{2}
\]

---

## 📄 C Code: `wigner_eckart_gsl.c`

```c
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_coupling.h>

#define TWICE(x) ((int)(2 * (x)))

int main(void) {
    double j = 0.5, m = 0.5;
    int two_j = TWICE(j), two_m = TWICE(m);
    int k = 1, q = 0;
    int two_jp = two_j, two_mp = two_m;

    double phase = pow(-1.0, (two_jp - two_j - 2*q)/2.0);
    double cg = gsl_sf_coupling_3j(two_j, TWICE(k), two_jp,
                                    two_m, TWICE(q), -two_mp);
    cg *= phase * sqrt(two_jp + 1.0);

    double red = sqrt(3.0) / 2.0;
    double result = red * cg;

    printf("Clebsch–Gordan coefficient: %.6f\n", cg);
    printf("Reduced matrix element:     %.6f\n", red);
    printf("Full matrix element:        %.6f\n", result);  // should be 0.5

    return 0;
}
```

---

## 🧪 Compile and Run

```bash
gcc wigner_eckart_gsl.c -lgsl -lgslcblas -lm -o wigner
./wigner
```

### ✅ Output

```
Clebsch–Gordan coefficient: 0.577350
Reduced matrix element:     0.866025
Full matrix element:        0.500000
```

---

## 📘 Why `cg *= phase * sqrt(2j'+1)`?

GSL gives Wigner 3j symbols. To convert them to **Clebsch–Gordan coefficients**, use:

\[
C^{j', m'}_{j, m; k, q} = (-1)^{j' - j - q} \cdot \sqrt{2j' + 1} \cdot
\begin{pmatrix}
j & k & j' \\
m & q & -m'
\end{pmatrix}
\]

This formula is implemented in:

```c
cg *= phase * sqrt(two_jp + 1.0);
```

---

## ✅ Conclusion

- You can compute matrix elements using Wigner–Eckart with GSL.
- Use `gsl_sf_coupling_3j` to get 3j symbols and convert to CG coefficients.
- Multiply by your reduced matrix element to get the full answer.

