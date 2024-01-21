#include<iostream>
using namespace std;
#include "example/case0.h"
#include "example/case0_2.h"
#include "example/case0_3.h"
#include "example/case0_4.h"
#include "example/case1.h"
#include "example/case2.h"
#include "example/case2_2.h"
#include "example/case2_3.h"
#include "example/case3.h"
#include "example/case3_2.h"
//
int main()
{   
    runcase0(); // Sec0: (1) 打印1d 数组，
    runcase0_2(); // Sec0: (2) 1d 数组转化为 2d 数组 in an 
    runcase0_3(); // Sec0: (3) 1d 指针数组转化为 1d Vec 容器
    runcase0_4(); // Sec0: (4) 1d 数组转化为 2d Vec
    // 
    runcase1(); // sec1: (1) 2d 指针数组的打印和取其中max
    runcase2(); // Sec2: (1) 2d 指针数组转化为 2d Vec 容器
    runcase2_2(); // Sec2: (2) 2d 指针数组转化为 1d Vec 容器
    runcase2_3(); // sec2: (3) 2d 指针数组转化为 1d 数组
    //
    //
    runcase3(); // Sec3: (1) 2d vec 的打印
    runcase3_2(); // Sec3: (2) 2d vec convert to 1d vec. 
    //    
    //
    return 0;
}