#include <vector>
#include <iostream>

using namespace std;

void vec_initialize();
void intvec_print_by_cite(std::vector<std::vector<int> >&, int, int);
void intvec_print_by_l1pointer(std::vector<std::vector<int> >*, int, int);
void doubvec_print_by_cite(std::vector<std::vector<double> >&, int, int);
void doubvec_print_by_l1pointer(std::vector<std::vector<double> >*, int, int);

