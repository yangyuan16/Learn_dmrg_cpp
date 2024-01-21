#include<iostream>
using namespace std;
#include "anyarray2d.h"
#include <vector>
#include <cstring>
//
anyarray2d::anyarray2d(/* args */)
{
}

anyarray2d::~anyarray2d()
{
}
//-----------------打印数组 和 取得数组中的最大值------------------------------------------------------------
//  传入数组二级指针
// 1. print array
template <typename T>
void anyarray2d::print_by_l2pointer(T **array, int n, int m)
{
    for (int i = 0; i < n; i++)
    {
    for (int j = 0; j < m; j++)
    {
            cout << array[i][j] << " ";
    }
    cout << endl;
    }
}

// 模板函数必须在实现cpp文件中追加对应的特化方式，或者将模板函数的实例和声明全部写在 .h 文件中 
template void anyarray2d::print_by_l2pointer<double>(double**, int, int);
template void anyarray2d::print_by_l2pointer<int>(int**, int, int); 

// 2. query max elements in an 2d array
template <typename T>
void anyarray2d::get_max_by_l2pointer(T **array, int n, int m)
{
    T max = array[0][0];
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            if(max < array[i][j]) max = array[i][j];
    }
    }
    cout << "max: " << max << endl;
}
template void anyarray2d::get_max_by_l2pointer<double>(double**, int, int);
template void anyarray2d::get_max_by_l2pointer<int>(int**, int, int);  
//
//  传入数组一级指针
//  1. print array 
template <typename T>
void anyarray2d::print_by_l1pointer(T *array, int n, int m)
{
    for (int i = 0; i<n; i++)
    {
        for (int j = 0; j<m; j++)
        {
            cout << array[i * m + j] << " ";
        }
        cout << endl;
    }
}
template void anyarray2d::print_by_l1pointer<double>(double*, int, int);
template void anyarray2d::print_by_l1pointer<int>(int*, int, int);

//  2.query max elements in an 2d array
template <typename T>  
void anyarray2d::get_max_by_l1pointer(T *array, int n, int m)
{
    T max = array[0];
    for (int i = 0; i < n * m; i++){
        if(max < array[i]) max = array[i];
    }
    cout << "max: " << max << endl;
}
template void anyarray2d::get_max_by_l1pointer<double>(double*, int, int);
template void anyarray2d::get_max_by_l1pointer<int>(int*, int, int);

//-----------------2d array convert to 2d vector---------------------------------
template <typename T> 
std::vector<std::vector<T>>anyarray2d::arr2dtovec2d(T *arr2d, int nr, int nc)
{   
    vector<T> v; // 定义一维 vec
    vector<vector<T>> vec2d;  
    for (int i = 0; i<nr; i++)
    {   
        v.clear(); // 子数组返回时要清零
        for (int j = 0; j<nc; j++)
        {
            cout << arr2d[i * nc + j] << " ";
            v.push_back(arr2d[i * nc + j]);
        }
        cout << endl;
        vec2d.push_back(v);

    }
    return vec2d; 
}
template vector<vector<double>> anyarray2d::arr2dtovec2d<double>(double *arr2d, int nr, int nc);
template vector<vector<int>> anyarray2d::arr2dtovec2d<int>(int *arr2d, int nr, int nc);

//--------------------2d array convert to 1d vector----------------------
template <typename T>
std::vector<T> anyarray2d::arr2dtovec1d(T *arr2d, int rows, int cols)
{
    std::vector<T> vector1D;

    for (int i = 0; i<rows; i++){
        for(int j =0; j<cols; j++){
            vector1D.push_back(arr2d[i * cols + j]);
        }
    }
    return vector1D;
}
template vector<double> anyarray2d::arr2dtovec1d<double>(double *arr2d, int rows, int cols);
template vector<int> anyarray2d::arr2dtovec1d<int>(int *arr2d, int rows, int cols);

//---------------------2d array convert to 1d array---------------
template <typename T>
void anyarray2d::arr2dtoarr1d(T** array2D, int rows, int cols, T* array1D){
    for (int i = 0; i < rows; i++)
    {
        // Use pointer arithmetic to efficiently copy the row elements to 1D array
        std::memcpy(array1D + i * cols, array2D[i], cols * sizeof(T));
    }
}
template void anyarray2d::arr2dtoarr1d<double>(double** array2D, int rows, int cols, double* array1D);
template void anyarray2d::arr2dtoarr1d<int>(int** array2D, int rows, int cols, int* array1D);









//
//
//
//-------------------------------------------------------------------
// 输入一个二维数组，获取数组的行列信息
//-------------------------------------------------------------------
// 1  传入2级指针
/*
template <typename T>
void anyarray2d_info_size_by_name(T arr)
{
    cout << "二维数组占用的内存空间为: " << sizeof(arr) << endl << endl;
    cout << "二维数组第一行占用的内存空间为: " << sizeof(arr[0]) << endl << endl;
    cout << "二维数组第一个元素占用的内存空间为:" << sizeof(arr[0][0]) << endl << endl;
    cout << "二维数组的行数为: " << sizeof(arr) / sizeof(arr[0]) << endl << endl;
    cout << "二维数组的列数为: " << sizeof(arr[0]) / sizeof(arr[0][0]) << endl << endl;
    cout << "二维数组的首地址（十六进制）为: " << arr << endl << endl;
    cout << endl

}
template void anyarray2d_info_size_by_name<double>(double);
template void anyarray2d_info_size_by_name<int>(int);
*/

