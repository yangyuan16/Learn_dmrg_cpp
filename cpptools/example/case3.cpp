#include <iostream>
using namespace std;
#include "case3.h"
#include "./../array/anyarray2d.h"
#include "./../array/array2d.h"
#include "./../vec/vec.h"
#include <vector>

void runcase3()
{   
    const int N = 2; //要改变二维数组的行
    const int M = 3; //二维数组的列
    //int a[N][M] = {{1, 2, 3}, {4, 5, 6}}; // int 型数组
    double b[N][M] = {{1.1,2.1,3.1}, {4.1,5.1,6.1}}; // double 型数组
    double *pl1_d = b[0]; // 声明一个 double 型一级指针
    
   
    vector<vector<double>> vec2d2; // 声明一个 double 2d vec
    vec2d2 = anyarray2d().arr2dtovec2d(pl1_d, N, M);
    cout << "--- print elements in 2d vec " << endl;
    vec().doubvec_print_by_cite(vec2d2, N, M);

    //vec().doubvec_print_by_l1pointer(&vec2d2, N, M);
    //doubvec_print_by_l1pointer(&vec2d2, N, M); // 打印 2d vec, 采用指针的方式传参数
    //  传递 vector 最好采用地址的方式
}