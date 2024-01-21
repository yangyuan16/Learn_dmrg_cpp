#ifndef _ANYARRAY2D_H_
#define _ANYARRAY2D_H_

#include <vector>
class anyarray2d
{
private:
    /* data */
public:
    anyarray2d(/* args */);
    //
    //--------------------------------------------------------------------
    // 打印数组 和 取得数组中的最大值
    //--------------------------------------------------------------------
    // 传入2级指针
    template <typename T>
    void print_by_l2pointer(T**, int, int); // print array

    template <typename T>
    void get_max_by_l2pointer(T**, int, int); // query max elements
    //---------------------------------------------------------------------
    // 传入1级指针
    template <typename T>
    void print_by_l1pointer(T*, int, int); // print array

    template <typename T>
    void get_max_by_l1pointer(T*, int, int); // query max elements
    //
    //-------------- array convert to vector---------------------------    
    // 1. 2d array convert 2d vector
    template <typename T>
    std::vector<std::vector<T>> arr2dtovec2d(T*, int, int);
    // 2. 2d array convert 1d vector
    template <typename T>
    std::vector<T> arr2dtovec1d(T*, int rows, int cols);
    //
    //------------- 2d array conver to 1d array
    template <typename T>
    void arr2dtoarr1d(T** array2D, int rows, int cols, T* array1D);
    ~anyarray2d();

};

#endif
//std::vector<std::vector<int>> intarr2dtovec2d(int*, int, int);
//std::vector<std::vector<double>> doubarr2dtovec2d(double*, int, int);