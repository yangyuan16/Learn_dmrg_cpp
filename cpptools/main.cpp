#include<iostream>
using namespace std;
//
//void fun(int** , int, int);

#include "array/array2d.h"
#include "array/anyarray2d.h"
int main()
{
    
    const int N = 2; //要改变二维数组的行和列
    const int M = 3; //只需改变这里即可
    int a[N][M] = {{1, 2, 3}, {4, 5, 6}}; // int 型数组
    double b[N][M] = {{1.1,2.1,3.1}, {4.1,5.1,6.1}}; // double 型数组

    int *pl2[N]; // 定义一个int型二级指针
    for (int i = 0; i<N; i++)
    {
        pl2[i] = &a[i][0]; // 为二级指针赋行地址值 // 也可以写为 p[i] = a[i];
    }

    int *pl1 = a[0]; // 定义一个一级指针
    
    // 采用 二级指针的方法获取 int型数组的最大值以及打印数组
    Array2d().intarray2d_get_max_by_l2pointer(pl2, N, M); //函数调用时传入指针 p
    Array2d().intarray2d_print_by_l2pointer(pl2, N, M);
    
    // 采用 一级指针的方法获取 int型数组的最大值以及打印数组
    Array2d().intarray2d_get_max_by_l1pointer(pl1, N, M);
    Array2d().intarray2d_print_by_l1pointer(pl1, N, M);

    
    double *pl2_d[N]; // 定义一个double 型二级指针 行指针
    pl2_d[0] = &b[0][0];
    pl2_d[1] = &b[1][0];

    double *pl1_d = b[0]; // 声明一个 double 型一级指针
    
    anyarray2d_print_by_l2pointer(pl2_d, N, M); // 采用 二级指针的方法获取打印任意类型数组
    anyarray2d_get_max_by_l2pointer(pl2_d, N, M); // 采用 二级指针的方法获取任意型数组的最大值

    anyarray2d_print_by_l1pointer(pl1_d, N, M); // 采用 一级指针的方法获取打印任意类型数组
    anyarray2d_get_max_by_l1pointer(pl1_d, N, M); // 采用 一级指针的方法获取任意型数组的最大值
    

    // 通过二级列指针打印数组
    double c[3][4] = {{1.1,2.1,3.1,4.1},{5.1,6.1,7.1,8.1},{9.1,10.1,11.1,12.1}};
    double *ppc[4];
    ppc[0] = &c[0][0];
    ppc[1] = &c[1][0];
    ppc[2] = &c[2][0];
    anyarray2d_print_by_l2pointer(ppc,3, 4);
    

    return 0;
}

