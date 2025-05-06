#include"parameter.h"
#include"subcommon.h"

class Sub:public SubCommon {

public:

        Sub( Parameter &para );           //Initialize blocks
        Sub( const int &i, const int &lsys);//Read in angular space from hard disk
	Sub(const int &i, const int &lsys, Sub &space);	//Read in operators

        Sub(const int &i, Parameter &para, Sub &old);   //New space for both sys and env

        Sub(char &sigh, Parameter &para, Sub &sys, Sub &sysnew_space);//Sys increase for Infinite and Finite sweep 
	Sub(Parameter &para, char &sigh, Sub &sys, Sub &sysnew_space);//Env increase for Infinite
	Sub(Sub &sys, Sub &sysnew_space, Parameter &para, char &sigh);//Env increase for Infinite

        Sub(Sub &space);    //Create Space 3 

        Sub(Parameter &para, Sub &density);     //Create space for truncation space in infinite sweep
	Sub(Parameter &para, Sub &density, const int &optimalstate);

        Sub(Parameter &para, char &sign, Sub &density, Sub &trun_space, Sub &old); //truncate operators


        Sub(Parameter &para, char *nooperator, Sub &old);//For the shrinking blocks in the finite sweep with only creating the new Hilbert space for diagonalization, but not to create the new operators.

	~Sub();                         //Delete the sub functions

private:

	int StartSite, mmin, SubDim;
        char trans_N, trans_T;
        double alpha, beta;

	inline void Find_Block_Num_hole(Parameter &para, Sub &old);
	inline void Find_Num_hole_block(Parameter &para, Sub &old);
	inline void Find_Num_J_hole_block(Parameter &para, Sub &old);
	inline void FindSpace_J_block(Parameter &para, Sub &old);

	inline void New_T(const int &block, Parameter &para, Sub &space, Sub &old, int *table_old, int *table_new);
        inline void NewS_Dia(const int &block, Parameter &para, Sub &space, Sub &old, int *table_old, int *table_new);
        inline void NewS_M_Dia(const int &block, Parameter &para, Sub &space, Sub &old, int *table_old, int *table_new);
	inline void New_NN(Parameter &para, Sub &space, Sub &old, int *table_old, int *table_new);
	inline void New_H_Sys(Parameter &para, Sub &space, Sub &old);
	inline void New_H_Env(Parameter &para, Sub &space, Sub &old);

        inline void Truncate(Parameter &para, Sub &old, const int &optimalstate);
        inline void Truncate_H(Sub &density, Sub &old);
        inline void Truncate_NN(Sub &density, Sub &old);
        inline void Truncate_T(Sub &density, Sub &old);
        inline void Truncate_S_Dia(Sub &density, Sub &old);
        inline void Truncate_S_M_Dia(Sub &density, Sub &old);

};
