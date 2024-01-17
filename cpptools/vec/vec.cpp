#include<iostream>
using namespace std;
#include<vector>

//----------------------------------------------------------------------------
// Sec0: initialize a 2d vector
void vec_initialize()
{
    
    // 创建 2*3 的vector 容器 v， 初始值均初始化为 0 1 2 1 2 3
    cout << "initialization of 2d vector" << endl;
    std::vector<std::vector<int> > v(2, std::vector<int>(3,0));
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<3; j++)
            v[i][j] = i+j;
    }
    // 打印v的地址
    cout<<"&v: "<<&v<<endl;
    // 打印v[i]的地址 （即第一层vector的地址）
    cout<<"&v[i]:"<<endl;
    for(int i=0; i<2; i++)
        cout<<&v[i]<<endl;
    // 打印 v 的各元素地址
    cout<<"&v[i][j]:"<<endl;
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<3;j++)
            cout<<&v[i][j]<<" ";
        cout<<endl;
    }
    cout<<"----------------------------------"<<endl;
    // 打印 v 的各元素值
    cout<<"v[i][j]:"<<endl;
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<3; j++)
            cout<<v[i][j]<<" ";
        cout<<endl;
    }
}
//
//----------------------------------------------------------------------------
// Sec1: (1) print the site and elements of int 2d vec by the way of cite
void intvec_print_by_cite(std::vector<std::vector<int> >& vec, int nr, int nc)
{
    cout<<"---print the site and elements of int 2d vec---"<<endl;
    //打印 vec 的地址
    cout<<"site of 2dvec, &vec: "<<&vec<<endl;
    //打印 vec[i]的地址(即第一层 vector 的地址)
    cout<<"first layer of site of 2d vec, &vec[i]:"<<endl;
    for(int i=0; i<nr; i++)
        cout<<&vec[i]<<endl;
    // 打印vec的个元素地址
    cout<<"site of 2d vec elements, &vec[i][j]: "<<endl;
    for(int i=0; i<nr; i++)
    {
        for(int j=0; j<nc; j++)
            cout<<&vec[i][j]<<" ";
        cout<<endl;
    }
    cout<<"----print elements of int 2d vec----"<<endl;
    // 打印vec的各个元素值
    cout<<"elements of 2d vec, vec[i][j]:"<<endl;
    for(int i=0; i<nr; i++)
    {
        for(int j=0; j<nc; j++)
            cout<<vec[i][j]<<" ";
        cout<<endl;
    }
}
// Sec1: (2) print the site and elements of int 2d vec by the way of pointer
void intvec_print_by_l1pointer(std::vector<std::vector<int> > *vec, int nr, int nc)
{
    cout<<"---print the site and elements of int 2d vec---"<<endl;
    //打印 vec 的地址
    cout<<"site of 2dvec, &vec: "<<&vec<<endl;
    //打印 vec[i]的地址(即第一层 vector 的地址)
    cout<<"first layer of site of 2d vec, &vec[i]:"<<endl;
    for(int i=0; i<nr; i++)
        cout<<&(*vec)[i]<<endl;
    // 打印vec的个元素地址
    cout<<"site of 2d vec elements, &vec[i][j]: "<<endl;
    for(int i=0; i<nr; i++)
    {
        for(int j=0; j<nc; j++)
            cout<<&(*vec)[i][j]<<" ";
        cout<<endl;
    }
    cout<<"-----print elements of int 2d vec--------"<<endl;
    // 打印vec的各个元素值
    cout<<"elements of 2d vec, vec[i][j]:"<<endl;
    for(int i=0; i<nr; i++)
    {
        for(int j=0; j<nc; j++)
            cout<<(*vec)[i][j]<<" ";
        cout<<endl;
    }
}
// 
// Sec1: (3) print the site and elements of double 2d vec by the way of cite
void doubvec_print_by_cite(std::vector<std::vector<double> >& vec, int nr, int nc)
{
    cout<<"---print the site and elements of int 2d vec---"<<endl;
    //打印 vec 的地址
    cout<<"site of 2dvec, &vec: "<<&vec<<endl;
    //打印 vec[i]的地址(即第一层 vector 的地址)
    cout<<"first layer of site of 2d vec, &vec[i]:"<<endl;
    for(int i=0; i<nr; i++)
        cout<<&vec[i]<<endl;
    // 打印vec的个元素地址
    cout<<"site of 2d vec elements, &vec[i][j]: "<<endl;
    for(int i=0; i<nr; i++)
    {
        for(int j=0; j<nc; j++)
            cout<<&vec[i][j]<<" ";
        cout<<endl;
    }
    cout<<"----print elements of int 2d vec----"<<endl;
    // 打印vec的各个元素值
    cout<<"elements of 2d vec, vec[i][j]:"<<endl;
    for(int i=0; i<nr; i++)
    {
        for(int j=0; j<nc; j++)
            cout<<vec[i][j]<<" ";
        cout<<endl;
    }
}
// Sec1: (4) print the site and elements of double 2d vec by the way of pointer
void doubvec_print_by_l1pointer(std::vector<std::vector<double> > *vec, int nr, int nc)
{
    cout<<"---print the site and elements of int 2d vec---"<<endl;
    //打印 vec 的地址
    cout<<"site of 2dvec, &vec: "<<&vec<<endl;
    //打印 vec[i]的地址(即第一层 vector 的地址)
    cout<<"first layer of site of 2d vec, &vec[i]:"<<endl;
    for(int i=0; i<nr; i++)
        cout<<&(*vec)[i]<<endl;
    // 打印vec的个元素地址
    cout<<"site of 2d vec elements, &vec[i][j]: "<<endl;
    for(int i=0; i<nr; i++)
    {
        for(int j=0; j<nc; j++)
            cout<<&(*vec)[i][j]<<" ";
        cout<<endl;
    }
    cout<<"-----print elements of int 2d vec--------"<<endl;
    // 打印vec的各个元素值
    cout<<"elements of 2d vec, vec[i][j]:"<<endl;
    for(int i=0; i<nr; i++)
    {
        for(int j=0; j<nc; j++)
            cout<<(*vec)[i][j]<<" ";
        cout<<endl;
    }
}
//-------------------------------------------------------------------------





















/*
int main()
{
    // 创建 2*3 的vector 容器 v， 初始值均初始化为 0 1 2 1 2 3
    std::vector<std::vector<int> > v(2, std::vector<int>(3,0));
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<3; j++)
            v[i][j] = i+j;
    }
    // 打印v的地址
    cout<<"&v: "<<&v<<endl;
    // 打印v[i]的地址 （即第一层vector的地址）
    cout<<"&v[i]:"<<endl;
    for(int i=0; i<2; i++)
        cout<<&v[i]<<endl;
    // 打印 v 的各元素地址
    cout<<"&v[i][j]:"<<endl;
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<3;j++)
            cout<<&v[i][j]<<" ";
        cout<<endl;
    }

    cout<<"----------------------------------"<<endl;
    // 打印 v 的各元素值
    cout<<"v[i][j]:"<<endl;
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<3; j++)
            cout<<v[i][j]<<" ";
        cout<<endl;
    }

    //function1(v);


    cout<<"------------------------------------"<<endl;
    // 打印 v 的各元素值
    cout<<"v[i][j]: "<<endl;
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<3; j++)
            cout<<v[i][j]<<" ";
        cout<<endl;
    }

    function2(v);

    cout<<"-------------------------------------"<<endl;
    // 打印v的个元素值
    cout<<"v[i][j]: "<<endl;
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<3; j++)
            cout<<v[i][j]<<" ";
        cout<<endl;
    }

    function3(&v);

    cout<<"-------------------------------------"<<endl;
    // 打印v的个元素值
    cout<<"v[i][j]: "<<endl;
    for(int i=0; i<2; i++)
    {
        for(int j=0; j<3; j++)
            cout<<v[i][j]<<" ";
        cout<<endl;
    }

    return 0;
}
*/