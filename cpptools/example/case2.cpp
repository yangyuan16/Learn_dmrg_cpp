#include <iostream>
using namespace std;
#include "case2.h"
#include "./../array/anyarray2d.h"
#include "./../array/array2d.h"

//---------example fun. 2d array transfer to 2d vec ------------------
void eg_arr2dtovec2d(int *pl1, double *pl1_d, int nr, int nc) // 采用一级指针方式传入数组
{
    vector<vector<int>> vec2d1; // 声明一个 int 2d vec
    vec2d1 = anyarray2d().arr2dtovec2d(pl1, nr, nc);
    vector<vector<double>> vec2d2; // 声明一个 double 2d vec
    vec2d2 = anyarray2d().arr2dtovec2d(pl1_d, nr, nc); 
}

void runcase2()
{
    const int N = 2; //要改变二维数组的行
    const int M = 3; //二维数组的列
    int a[N][M] = {{1, 2, 3}, {4, 5, 6}}; // int 型数组
    double b[N][M] = {{1.1,2.1,3.1}, {4.1,5.1,6.1}}; // double 型数组
    // --------- 声明和定义 int 型二级指针和 int 型一级指针
    int *pl2[N]; // 定义一个int型二级指针
    for (int i = 0; i<N; i++)
    {
        pl2[i] = &a[i][0]; // 为二级指针赋行地址值 // 也可以写为 p[i] = a[i];
    }
    int *pl1 = a[0]; // 定义一个一级指针
    // --------- 声明和定义 double 型 二级指针 和 double 新
    double *pl2_d[N]; // 定义一个double 型二级指针 行指针
    pl2_d[0] = &b[0][0];
    pl2_d[1] = &b[1][0];
    double *pl1_d = b[0]; // 声明一个 double 型一级指针
    //
    // Sec2: (1) 2d 指针数组转化为 2d Vec 容器
    cout << "---------------convert 2d array to 2d vector--------------" << endl;
    eg_arr2dtovec2d(pl1,pl1_d,N,M);
    cout << "----------------------------------------------------------" << endl;
}