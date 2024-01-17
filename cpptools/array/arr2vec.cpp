#include <iostream>
using namespace std;
#include <vector>
//----------------------------------------------------
// arr2d to vec2d; arr2d 采用 一级指针的方式传入; 
// 1. int 类型
vector<vector<int>> intarr2dtovec2d(int *arr2d, int nr, int nc)
{   
    vector<int> v; // 定义一维 vec
    vector<vector<int>> vec2d;  
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
// 
// 2. double 类型
vector<vector<double>> doubarr2dtovec2d(double *arr2d, int nr, int nc)
{   
    vector<double> v; // 定义一维 vec
    vector<vector<double>> vec2d;  
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
//