#include "LSMatrixDiagonalization.h"

class Super:public LSMatrixDiagonalization {

public:
  //------Construct the space for wave function blocks
	Super(Parameter &para, Sub *sys_space1, Sub *env_space1, Sub *sysnew_space1, Sub *envnew_space1);
  
        Sub *sys_space;
        Sub *env_space;
        Sub *sysnew_space;
        Sub *envnew_space;

  //------Construct the space for diagonalization
        Super(char &sign, Parameter &para, Sub *sys1, Sub *env1, Sub *sysnew1, Sub *envnew1, Super &space);
  
        Sub *sys;
        Sub *env;
	Sub *sysnew;
	Sub *envnew;

	int QuantumNumber_hole, QuantumNumber_J;

	int Dim;
        int BlockNumber_for_TargetBlock;
        int BlockNumber_for_TargetBlock_config_3;

        int *H_sys, *H_env, *H_sysnew, *H_envnew, *J_sys, *J_env, *J_sysnew, *J_envnew, *Dim_block;
	int *H_sys_config_3, *H_env_config_3, *H_ns_config_3, *H_ne_config_3, *J_sys_config_3, *J_env_config_3, *J_sysnew_config_3, *J_envnew_config_3, *Dim_block_config_3;

        int *Table_1to2_Num;
        int *Table_1to2_Site;
        int **Table_2to1;

        double *WaveFunction;
        double **WaveFunction_block;

        double *six_j_1_sys_n; 
        double *six_j_1_e_env;
        double **nine_j_config_2;
	double **nine_j_config_2_inverse;

        double **six_j_T_config_3;
	double **six_j_J_config_3;
        double *six_j_J_config_3_diaele;
        double *six_j_2_config_3;
        double **nine_j_config_3;
        double **nine_j_config_3_inverse;

        void getMatrixDiagElement(double *f, const int &dim);
        void H_V(const double *f, double *g, const int &dim);
        void NormalizedCopy(const double *f1, double *f2);

	Super( Sub *sys_space1, Sub *env_space1, Sub *sysnew_space1, Sub *envnew_space1, Parameter &para );

        ~Super();       //Delete super functions

private:

	int N, StartSite, inc, index, SiteNum, site_s, site_e, old_J, dimension;
        int operator_number_T_sys, operator_number_T_env;
        int operator_number_J_sys, operator_number_J_env;
	int operator_number_N_sys, operator_number_N_env;
        int Dim_config_3;

        char destruct, side_L, side_R, uplo, trans_N, trans_T;

        double beta, alpha_p, beta_p, alpha;

        int *Table_T, *Table_T_sys, *Table_T_env;
        int *Table_J, *Table_J_sys, *Table_J_env;
        int *Table_N, *Table_N_sys, *Table_N_env;
        double *Interaction_T, *Interaction_J, *Interaction_N;

        double **f1;  double **g1;
        double **f2;  double **g2;
        double **f3;  double **g3;

	double *** T_sys;  double *** T_env;
        double *** S_Dia_sys;
        double *** S_M_Dia_sys;
        double *** S_Dia_env;
        double *** S_M_Dia_env;
	double *** NN_sys;
	double *** NN_env;

        inline void AllocateBlockNumber(Parameter &para);
        inline void AllocateBlockNumber_config_3(Parameter &para);

	inline void Allocate_6j_config_1_2(Super &space);
	inline void Allocate_9j_config_2(Parameter &para, Super &space);
	inline void Allocate_6j_config_3();
	inline void Allocate_9j_config_3(Parameter &para, Super &space);

        inline void H_Sys_and_Env();
        inline void H_Sys_Ns_T();
	inline void H_Sys_Ns_J();
	inline void H_Sys_Ns_N();
        inline void H_Env_Ne_T();
	inline void H_Env_Ne_J();
	inline void H_Env_Ne_N();
        inline void H_Sys_Ne_T();
	inline void H_Sys_Ne_J();
	inline void H_Sys_Ne_N();
        inline void H_Env_Ns_T();
        inline void H_Env_Ns_J();
	inline void H_Env_Ns_N();
        inline void H_Ns_Ne_T();
        inline void H_Ns_Ne_J();
	inline void H_Ns_Ne_N();
        inline void H_Sys_Env_T();
	inline void H_Sys_Env_J();
	inline void H_Sys_Env_N();

};
