#include"conjugate.h"

class SuperEnergy:public Conjugate {
public:
    //--Diagonalization of Hamiltonian by Matrix-Vector multiplication
        SuperEnergy(Parameter &para, Super &sup_space, Super &sup, const double &lan_precision);//infinite
        SuperEnergy(Parameter &para, Super &sup_space, Super &sup, char &sign, const double &lan_precision);//finite

        SuperEnergy(Parameter &para, Super &sup, Sub &old, Sub &trun, char &sign);

        SuperEnergy(Parameter &para, Super &sup_space, Super &sup, const int &n, const double &lan_precision);
        SuperEnergy(Super &sup, Sub &trun, const int &n);
        ~SuperEnergy();

        double eigenvalue;

private:

        inline void SupConjugate(Super &sup_space, Super &sup);

  //------Davidson method to diagonalization---------------------------------------------------------------------
        int n, eigennum, totalbase, minibase, maxiteration, input;
        double precision;
        inline void Davidson_Parameter(const int &sign, Super &sup, const double &lan_precision);
        inline void Davidson(Parameter &para, Super &sup_space, Super &sup);

  //------Jacobi-Davidson method to diagonalization:subroutine DPJDREVCOM of JADAMILU----------------------------
        int N, LX, NEIG, ISEARCH, NINIT, MADSPACE, INFO, ITER, IJOB, NDX1, NDX2, IPRINT;
        double SIGMA, TOL, SHIFT, DROPTOL, MEM, GAP;
        int *ICNTL, *JA, *IA;
        double *A, *X, *EIGS, *RES;
        inline void DPJDREVCOM_Parameter(Super &sup, const int &ninit);
        inline void DPJDREVCOM(Parameter &para, Super &sup_space, Super &sup);

  //------First step to truncate the wavefunction----------------------------------------------------------------
        char trans_N, trans_T;
        int index, SubDim, J_min, J_max, J_num, position_old, position_new;
        int truncated_block_number, oldH_sys, oldJ_sys, oldH_env, oldJ_env, untruncated_block_number;
        double alpha, beta;
	int *H_trun, *J_trun, *H_old, *J_old, *H_new, *J_new, *truncated_block_dim;
        int *H_sysnew_untrun, *J_sysnew_untrun, *H_envnew_untrun, *J_envnew_untrun, *H_sys_untrun, *J_sys_untrun, *H_env_untrun, *J_env_untrun, *untruncated_block_dim;

	double **truncated_wave_function;
        double **untruncated_wave_function;

	double ***truncated_density_eigenvector;

        inline void Truncate_sysnew_density_eigenvector(Parameter &para, Super &sup, Sub &old, Sub &trun);
        inline void Truncate_envnew_density_eigenvector(Parameter &para, Super &sup, Sub &old, Sub &trun);

        int truncated_Block_Num_hole;
        int *truncated_Num_hole_block, *truncated_Num_J_hole_block, *truncated_Old_hole;
        int **truncated_Value_J_block, **truncated_Dim_J_block, **truncated_density_dim, **truncated_Old_J;

        int new_Block_Num_hole;
        int *new_Num_hole_block, *new_Num_J_hole_block;
        int **new_Value_J_block, **new_Dim_J_block;
        int ***new_Hole_blockOld, ***new_J_blockOld, ***new_Start;

        double **six_j_basis_transformation;



  //------Second and third wave function transformations and initialize the trial function-----------------------
	inline void Initialtrialfunction_Left_to_Right(Parameter &para, Super &sup_space);
        inline void Initialtrialfunction_Right_to_Left(Parameter &para, Super &sup_space);

};
