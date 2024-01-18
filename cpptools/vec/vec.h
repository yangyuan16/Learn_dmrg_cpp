#include <vector>
#include <string>
#include <iostream>
using namespace std;

class vec
{
private:
    /* data */
public:
    vec(/* args */);
    //------------ 1d vector  ----------------------------------------
    void initialize_1d(int );
    //------------ 2d vector ----------------------------------------
    void initialize_2d(int, int);
    // print 2d vector
    void intvec_print_by_cite(std::vector<std::vector<int> >&, int, int);
    //
    void intvec_print_by_l1pointer(std::vector<std::vector<int> >*, int, int);
    //
    void doubvec_print_by_cite(std::vector<std::vector<double> >&, int, int);
    //
    void doubvec_print_by_l1pointer(std::vector<std::vector<double> >*, int, int);
    //--------------------------------------------------------------------
    // convert2DTo1D vector
    template <typename T>
    void convert2DTo1D(const std::vector<std::vector<T>>& twoDArray, 
                        std::vector<T>& oneDArray);

    //
    ~ vec();
};


