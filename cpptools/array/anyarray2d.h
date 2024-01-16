//-----------------------------------------------------------------------------
// 打印数组 和 取得数组中的最大值
//------------------------------------------------------------------------
// 传入2级指针
template <typename T>
void anyarray2d_print_by_l2pointer(T**, int, int); // print array

template <typename T>
void anyarray2d_get_max_by_l2pointer(T**, int, int); // query max elements
//------------------------------------------------------------------------
// 传入1级指针
template <typename T>
void anyarray2d_print_by_l1pointer(T*, int, int); // print array

template <typename T>
void anyarray2d_get_max_by_l1pointer(T*, int, int); // query max elements
//
//-------------------------------------------------------------------
// 输入一个二维数组，获取数组的行列信息
//-------------------------------------------------------------------
// 1. 传入 二级指针
/*
template <typename T>
void anyarray2d_info_size_by_name(T);
*/