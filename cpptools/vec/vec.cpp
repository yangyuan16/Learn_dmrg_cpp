#include <iostream>
#include "vec.h"
using namespace std;
//
//
vec::vec(/* args */){}
vec::~ vec(){}
//

//----------------- 1d vector ---------------------------------------------------
void vec::initialize_1d(int N){
    std::vector<int> oneDArray(N);
}
//----------------- 2d vector ---------------------------------------------------
void vec::initialize_2d(int rows, int cols){
    std::vector<std::vector<int>> large2DArray(rows, std::vector<int>(cols, 0));
}
// Sec1: (1) print the site and elements of int 2d vec by the way of cite
void vec::intvec_print_by_cite(std::vector<std::vector<int> >& vec2d, int nr, int nc)
{
    cout<<"---print the site and elements of int 2d vec---"<<endl;
    //打印 vec 的地址
    cout<<"site of 2dvec, &vec: "<<&vec2d<<endl;
    //打印 vec[i]的地址(即第一层 vector 的地址)
    cout<<"first layer of site of 2d vec, &vec[i]:"<<endl;
    for(int i=0; i<nr; i++)
        cout<<&vec2d[i]<<endl;
    // 打印vec的个元素地址
    cout<<"site of 2d vec elements, &vec[i][j]: "<<endl;
    for(int i=0; i<nr; i++)
    {
        for(int j=0; j<nc; j++)
            cout<<&vec2d[i][j]<<" ";
        cout<<endl;
    }
    cout<<"----print elements of int 2d vec----"<<endl;
    // 打印vec的各个元素值
    cout<<"elements of 2d vec, vec[i][j]:"<<endl;
    for(int i=0; i<nr; i++)
    {
        for(int j=0; j<nc; j++)
            cout<<vec2d[i][j]<<" ";
        cout<<endl;
    }
}
// Sec1: (2) print the site and elements of int 2d vec by the way of pointer
void vec::intvec_print_by_l1pointer(std::vector<std::vector<int> > *vec2d, int nr, int nc)
{
    cout<<"---print the site and elements of int 2d vec---"<<endl;
    //打印 vec 的地址
    cout<<"site of 2dvec, &vec: "<<&vec2d<<endl;
    //打印 vec[i]的地址(即第一层 vector 的地址)
    cout<<"first layer of site of 2d vec, &vec[i]:"<<endl;
    for(int i=0; i<nr; i++)
        cout<<&(*vec2d)[i]<<endl;
    // 打印vec的个元素地址
    cout<<"site of 2d vec elements, &vec[i][j]: "<<endl;
    for(int i=0; i<nr; i++)
    {
        for(int j=0; j<nc; j++)
            cout<<&(*vec2d)[i][j]<<" ";
        cout<<endl;
    }
    cout<<"-----print elements of int 2d vec--------"<<endl;
    // 打印vec的各个元素值
    cout<<"elements of 2d vec, vec[i][j]:"<<endl;
    for(int i=0; i<nr; i++)
    {
        for(int j=0; j<nc; j++)
            cout<<(*vec2d)[i][j]<<" ";
        cout<<endl;
    }
}
// 
// Sec1: (3) print the site and elements of double 2d vec by the way of cite
void vec::doubvec_print_by_cite(std::vector<std::vector<double> >& vec2d, int nr, int nc)
{
    cout<<"---print the site and elements of int 2d vec---"<<endl;
    //打印 vec 的地址
    cout<<"site of 2dvec, &vec: "<<&vec2d<<endl;
    //打印 vec[i]的地址(即第一层 vector 的地址)
    cout<<"first layer of site of 2d vec, &vec[i]:"<<endl;
    for(int i=0; i<nr; i++)
        cout<<&vec2d[i]<<endl;
    // 打印vec的个元素地址
    cout<<"site of 2d vec elements, &vec[i][j]: "<<endl;
    for(int i=0; i<nr; i++)
    {
        for(int j=0; j<nc; j++)
            cout<<&vec2d[i][j]<<" ";
        cout<<endl;
    }
    cout<<"----print elements of int 2d vec----"<<endl;
    // 打印vec的各个元素值
    cout<<"elements of 2d vec, vec[i][j]:"<<endl;
    for(int i=0; i<nr; i++)
    {
        for(int j=0; j<nc; j++)
            cout<<vec2d[i][j]<<" ";
        cout<<endl;
    }
}
// Sec1: (4) print the site and elements of double 2d vec by the way of pointer
void vec::doubvec_print_by_l1pointer(std::vector<std::vector<double> > *vec2d, int nr, int nc)
{
    cout<<"---print the site and elements of int 2d vec---"<<endl;
    //打印 vec 的地址
    cout<<"site of 2dvec, &vec: "<<&vec2d<<endl;
    //打印 vec[i]的地址(即第一层 vector 的地址)
    cout<<"first layer of site of 2d vec, &vec[i]:"<<endl;
    for(int i=0; i<nr; i++)
        cout<<&(*vec2d)[i]<<endl;
    // 打印vec的个元素地址
    cout<<"site of 2d vec elements, &vec[i][j]: "<<endl;
    for(int i=0; i<nr; i++)
    {
        for(int j=0; j<nc; j++)
            cout<<&(*vec2d)[i][j]<<" ";
        cout<<endl;
    }
    cout<<"-----print elements of int 2d vec--------"<<endl;
    // 打印vec的各个元素值
    cout<<"elements of 2d vec, vec[i][j]:"<<endl;
    for(int i=0; i<nr; i++)
    {
        for(int j=0; j<nc; j++)
            cout<<(*vec2d)[i][j]<<" ";
        cout<<endl;
    }
}
//
template <typename T>
void vec::convert2DTo1D(const std::vector<std::vector<T>>& twoDArray, std::vector<T>& oneDArray)
{
    const int rows = twoDArray.size();
    const int cols = twoDArray[0].size();

    for (int i = 0; i< rows; i++)
    {
        std::copy(twoDArray[i].begin(), twoDArray[i].end(),oneDArray.begin() + i * cols);
    }
        
}
template void vec::convert2DTo1D<double>(const std::vector<std::vector<double>>& twoDArray, 
                                    std::vector<double>& oneDArray);
template void vec::convert2DTo1D<int>(const std::vector<std::vector<int>>& twoDArray, 
                                    std::vector<int>& oneDArray);
// 
