#include<iostream>
using namespace std;
//
#include "array/array2d.h"
#include "array/anyarray2d.h"
#include<vector>
#include "array/arr2vec.h"
#include "vec/vec.h"

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
//---------example fun. 2d array transfer to 2d vec ------------------
void eg_arr2dtovec2d(int *pl1, double *pl1_d, int nr, int nc) // 采用一级指针方式传入数组
{
    vector<vector<int>> vec2d1; // 声明一个 int 2d vec
    vec2d1 = intarr2dtovec2d(pl1, nr, nc);
    vector<vector<double>> vec2d2; // 声明一个 double 2d vec
    vec2d2 = doubarr2dtovec2d(pl1_d, nr, nc); 
}
//-----------------------------------------------------------
//
//
int main()
{
    // Sec1: 2d 指针数组的打印和取其中max
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

    eg_intarr_1(pl1, pl2, N, M); // 传的是整型一级和二级指针
    eg_doubarr_1(pl1_d, pl2_d, N, M); // 传的是double 型一级和二级指针
    //--------------------------------------------------------------   
    // Sec2: (1) 2d 指针数组转化为 2d Vec 容器
    eg_arr2dtovec2d(pl1,pl1_d,N,M);
    // Sec2: (2) 2d 指针数组转化为 1d Vec 容器
    // Sec2: (3) 1d 指针数组转化为 1d Vec 容器
    //--------------------------------------------------------------
    // Sec3: (1) 2d vec 的打印
    // 得到一个 2d vec: 
    vector<vector<double>> vec2d2; // 声明一个 double 2d vec
    vec2d2 = doubarr2dtovec2d(pl1_d, N, M);
    // vec_initialize(); 
    vec().doubvec_print_by_cite(vec2d2, N, M);
    //doubvec_print_by_cite(vec2d2, N, M); // 打印2d vec， 采用引用的方式进行传参数     
    vec().doubvec_print_by_l1pointer(&vec2d2, N, M);
    //doubvec_print_by_l1pointer(&vec2d2, N, M); // 打印 2d vec, 采用指针的方式传参数
    //  传递 vector 最好采用地址的方式
    // 
    const int rows = 1000;
    const int cols = 1000;
    // Initialize a large 2D array
    std::vector<std::vector<int>> large2DArray(rows, std::vector<int>(cols, 0));
    // Initialize a 1D array with the same total size
    std::vector<int> oneDArray(rows * cols);
    // Convert 2D array to 1D array efficiently
    vec().convert2DTo1D(large2DArray, oneDArray);
    //
    // Access elements in the 1D array
    cout << "First element in 1D array: " << oneDArray[0] << endl;
    cout <<"-----------------------------------------------------"<< endl;
    std::vector<double> vec1d1_a(N * M); // 初始化一个 1d vector
    for (int i=0; i<N*M; i++)
    {
        cout<<vec1d1_a[i] << endl;
    }
    
    vec().convert2DTo1D(vec2d2, vec1d1_a); // covert 2d to 1d
    for (int i=0; i<N*M; i++)
    {
        cout<<vec1d1_a[i] << endl;
    }

    return 0;
}