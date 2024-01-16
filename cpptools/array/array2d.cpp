#include<iostream>
using namespace std;
#include "array2d.h"

Array2d::Array2d(/* args */)
{
}

Array2d::~Array2d()
{
}

// ------------------------------------------------------------------------
//  传入的是二级指针
void Array2d::intarray2d_get_max_by_l2pointer(int **array, int n, int m)
{
    int max = array[0][0];
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            if(max < array[i][j]) max = array[i][j];
    }
    }
    cout << "max: " << max << endl;
}

void Array2d::intarray2d_print_by_l2pointer(int **array, int n, int m)
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
// ----------------------------------------------------------------------
//   传入的是一级指针
void Array2d::intarray2d_get_max_by_l1pointer(int *array, int n, int m)
{
    int max = array[0];
    for (int i = 0; i < n * m; i++){
        if(max < array[i]) max = array[i];
    }
    cout << "max: " << max << endl;
}

void Array2d::intarray2d_print_by_l1pointer(int *array, int n, int m)
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
// ------------------------------------------------------------------

