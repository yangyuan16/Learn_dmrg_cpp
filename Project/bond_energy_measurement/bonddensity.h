#include<fstream>
#include<sstream>

class BondDensity {

public:

	BondDensity( const int &lsys, Parameter &para, Super &sup );

        ~BondDensity();

private:

	char trans_N, trans_T, side_R, side_L, uplo;
	double alpha, beta;

  //----Create filestream for saving data
        string filename;
        ofstream fout;

  //----Create spaces for matrices and functions
	double **electron_density;

	double **WaveFunction_block;
	double **WaveFunction_config_2;
        double **WaveFunction_config_3;

	inline void CreateSpace( const int &lsys, Parameter &para );

	inline void CreateFunction( Super &sup );

	//------A_N-matrices
	int * A_N_Block_Num_hole;

	int ** A_N_Num_hole_block;
	int ** A_N_Num_J_hole_block;

	int *** A_N_Value_J_block;
	int *** A_N_Dim_J_block;

	int **** A_N_Hole_blockOld;
	int **** A_N_J_blockOld;
	int **** A_N_Start;

	int * A_N_Block_T;
	int ** A_N_Index_block_T;
	int ** A_N_Num_block_T;
	int *** A_N_J_block_T_bra;
	int *** A_N_J_block_T_ket;

	double ***** A_six_j_T;
	double ***** A_six_j_S_Dia_old;
        double **** A_six_j_S_Dia_n;
        double ***** A_six_j_S_M_Dia_old;
        double ***** A_six_j_S_M_Dia_n;
	double ***** A_six_j_H;

	//------B_N-matrices
	int * B_N_Block_Num_hole;

	int ** B_N_Num_hole_block;
	int ** B_N_Num_J_hole_block;

	int *** B_N_Value_J_block;
	int *** B_N_Dim_J_block;

	int **** B_N_Hole_blockOld;
	int **** B_N_J_blockOld;
	int **** B_N_Start;

	int * B_N_Block_T;
	int ** B_N_Index_block_T;
	int ** B_N_Num_block_T;
	int *** B_N_J_block_T_bra;
	int *** B_N_J_block_T_ket;

	double ***** B_six_j_T;
	double ***** B_six_j_S_Dia_old;
        double **** B_six_j_S_Dia_n;
        double ***** B_six_j_S_M_Dia_old;
        double ***** B_six_j_S_M_Dia_n;
	double ***** B_six_j_H;

	//------A-matrices
	int * A_Block_Num_hole;
	int ** A_Num_hole_block;
	int ** A_Num_J_hole_block;
	int *** A_Value_J_block;
	int *** A_Dim_J_block;
	int *** A_density_dim;
	int ** A_Old_hole;
	int *** A_Old_J;
	double **** A;

	//------B-matrices
	int * B_Block_Num_hole;
	int ** B_Num_hole_block;
	int ** B_Num_J_hole_block;
	int *** B_Value_J_block;
	int *** B_Dim_J_block;
	int *** B_density_dim;
	int ** B_Old_hole;
	int *** B_Old_J;
	double **** B;

	inline void Density_sys( const int &lsys, Parameter &para, Super &sup );
	inline void Density_env( const int &lsys, Parameter &para, Super &sup );
	inline void Density_ns_ne( const int &lsys, Parameter &para, Super &sup );
	inline void Density_sys_env( const int &lsys, Parameter &para, Super &sup );

	double ***Ni_old_A, ***Ni_new_A; 
	double ***Ni_old_B, ***Ni_new_B;

	inline void A_Initial_Ni( const int &num );
	inline void B_Initial_Ni( const int &num );

	inline void A_New_Ni( const int &num );
	inline void B_New_Ni( const int &num );

	inline void A_Truncate_Ni( const int &num );
	inline void B_Truncate_Ni( const int &num );

        double ***SiSj_old_A, ***SiSj_new_A;
        double ***SiSj_old_B, ***SiSj_new_B;

        inline void Initial_A_SiSj(const int &num);
        inline void Initial_B_SiSj(const int &num);

        inline void New_A_SiSj(const int &num);
        inline void New_B_SiSj(const int &num);

        inline void Truncate_A_SiSj(const int &num);
        inline void Truncate_B_SiSj(const int &num);

	inline void DeleteFunction( Super &sup );
	inline void DeleteSpace( const int &lsys, Parameter &para );

};
