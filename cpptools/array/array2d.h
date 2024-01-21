#ifndef _ARRAY2D_H_
#define _ARRAY2D_H_

#include<iostream>

class Array2d
{
private:
    /* data */
public:
    Array2d(/* args */);
    ~Array2d();
    void intarray2d_get_max_by_l2pointer(int**, int, int);
    void intarray2d_print_by_l2pointer(int**, int, int);

    void intarray2d_get_max_by_l1pointer(int*, int, int);
    void intarray2d_print_by_l1pointer(int*, int, int);

    //template <typename T>
    //void anyarray2d_get_max_by_l2pointer(T**, int, int);
};

#endif