#include<iostream>
using namespace std;
#include "anyarray2d.h"
//-----------------------------------------------------------------------------
// 打印数组 和 取得数组中的最大值
//-----------------------------------------------------------------------------
//  传入数组二级指针
// 1. print array
template <typename T>
void anyarray2d_print_by_l2pointer(T **array, int n, int m)
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
template void anyarray2d_print_by_l2pointer<double>(double**, int, int);
template void anyarray2d_print_by_l2pointer<int>(int**, int, int); 

// 2. query max elements in an 2d array
template <typename T>
void anyarray2d_get_max_by_l2pointer(T **array, int n, int m)
{
    T max = array[0][0];
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            if(max < array[i][j]) max = array[i][j];
    }
    }
    cout << "max: " << max << endl;
}
template void anyarray2d_get_max_by_l2pointer<double>(double**, int, int);
template void anyarray2d_get_max_by_l2pointer<int>(int**, int, int);  
//--------------------------------------------------------------------------------
//  传入数组一级指针
//  1. print array 
template <typename T>
void anyarray2d_print_by_l1pointer(T *array, int n, int m)
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
template void anyarray2d_print_by_l1pointer<double>(double*, int, int);
template void anyarray2d_print_by_l1pointer<int>(int*, int, int);

//  2.query max elements in an 2d array
template <typename T>  
void anyarray2d_get_max_by_l1pointer(T *array, int n, int m)
{
    T max = array[0];
    for (int i = 0; i < n * m; i++){
        if(max < array[i]) max = array[i];
    }
    cout << "max: " << max << endl;
}
template void anyarray2d_get_max_by_l1pointer<double>(double*, int, int);
template void anyarray2d_get_max_by_l1pointer<int>(int*, int, int);

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

