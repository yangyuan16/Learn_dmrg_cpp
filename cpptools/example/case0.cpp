#include <iostream>
using namespace std;
#include "case0.h"
#include "./../array/anyarray1d.h"

void runcase0()
{
    double arr[] = {1.1,1.2,1.3,1.4,1.5,1.6};
    int size = sizeof(arr) / sizeof(arr[0]);
    cout << "--------[print 1d array]--------" << endl;
    anyarray1d().printArray(arr, size);
}