#include <iostream>
using namespace std;
#include "case1.h"
#include "./../array/anyarray2d.h"
#include "./../array/array2d.h"
//
//-------- example fun. for int 2d array -----------
void eg_intarr_1(int *pl1, int **pl2, int N, int M)
{
    // 采用 二级指针的方法获取 int型数组的最大值以及打印数组
    Array2d().intarray2d_get_max_by_l2pointer(pl2, N, M); //函数调用时传入指针 p
    Array2d().intarray2d_print_by_l2pointer(pl2, N, M);
    
    // 采用 一级指针的方法获取 int型数组的最大值以及打印数组
    Array2d().intarray2d_get_max_by_l1pointer(pl1, N, M);
    Array2d().intarray2d_print_by_l1pointer(pl1, N, M);
}
//---------example fun. for double 2d array ---------
void eg_doubarr_1(double *pl1_d, double **pl2_d, int N, int M)
{
    // 采用 二级指针的方法获取打印任意类型数组
    anyarray2d().print_by_l2pointer(pl2_d, N, M); 
    anyarray2d().get_max_by_l2pointer(pl2_d, N, M); // 采用 二级指针的方法获取任意型数组的最大值

    // 采用 一级指针的方法获取打印任意类型数组
    anyarray2d().print_by_l1pointer(pl1_d, N, M); 
    anyarray2d().get_max_by_l1pointer(pl1_d, N, M); // 采用 一级指针的方法获取任意型数组的最大值
}
//
void runcase1(){
    // --------- 声明和定义 int 型二级指针和 int 型一级指针
    const int N = 2; //要改变二维数组的行
    const int M = 3; //二维数组的列
    int a[N][M] = {{1, 2, 3}, {4, 5, 6}}; // int 型数组
    double b[N][M] = {{1.1,2.1,3.1}, {4.1,5.1,6.1}}; // double 型数组
    //
    int *pl2[N]; // 定义一个int型二级指针
    for (int i = 0; i<N; i++)
    {
        pl2[i] = &a[i][0]; // 为二级指针赋行地址值 // 也可以写为 p[i] = a[i];
    }
    int *pl1 = a[0]; // 定义一个一级指针
    // --------- 声明和定义 double 型 二级指针 和 double 
    double *pl2_d[N]; // 定义一个double 型二级指针 行指针
    pl2_d[0] = &b[0][0];
    pl2_d[1] = &b[1][0];
    double *pl1_d = b[0]; // 声明一个 double 型一级指针
    //
    cout <<"------print 2d array elements and query max elements-------------"<<endl;
    eg_intarr_1(pl1, pl2, N, M); // 传的是整型一级和二级指针
    eg_doubarr_1(pl1_d, pl2_d, N, M); // 传的是double 型一级和二级指针
    cout <<"-----------------------------------------------------------------"<<endl;
}