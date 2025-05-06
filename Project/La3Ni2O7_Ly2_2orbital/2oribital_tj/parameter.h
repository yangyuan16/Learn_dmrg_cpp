class Parameter {

public:

	Parameter(char &sign, const int &spin, const int &total_h, const int &total_j, const int &number_hole, const int &statenokept, const int &n_u, const int &n_x, const int &n_y, const int &total_n, const double &t_a_xx, const double &j_a_xx, const double &n_a_xx, const double &t_a_zz, const double &j_a_zz, const double &n_a_zz, const double &t_a_xz, const double &j_a_xz, const double &n_a_xz, const double &t_b_xx, const double &j_b_xx, const double &n_b_xx, const double &t_b_zz, const double &j_b_zz, const double &n_b_zz, const double &t_b_xz, const double &j_b_xz, const double &n_b_xz, const double &t_ab_xx, const double &j_ab_xx, const double &n_ab_xx, const double &t_ab_zz, const double &j_ab_zz, const double &n_ab_zz,	
			const int &desiredzeros, const int &seedvalue);

 //------Model and DMRG parameters
 	char Sign;
        int Spin, Total_h, Total_J, Number_hole, StateNoKept, N_u, N_x, N_y, Total_N, untruncated_site, total_site;
	double T_a_xx, J_a_xx, N_a_xx, T_a_zz, J_a_zz, N_a_zz, T_a_xz, J_a_xz, N_a_xz;
	double T_b_xx, J_b_xx, N_b_xx, T_b_zz, J_b_zz, N_b_zz, T_b_xz, J_b_xz, N_b_xz;
	double T_ab_xx, J_ab_xx, N_ab_xx, T_ab_zz, J_ab_zz, N_ab_zz;
	int Desiredzeros, Seedvalue;

	int * QuantumNumber_hole;
	int * QuantumNumber_J;	
        
	int * randomList; 
 //------Table that stores the sites with interactions
 	//T coupling
        int *Table_T;
 
        int *Table_T_sys;
        int **Table_T_sys_site;
        int *Table_T_env;
        int **Table_T_env_site;

        double *Interaction_T; 

	//J coupling
        int *Table_J;
 
        int *Table_J_sys;
        int **Table_J_sys_site;

        int *Table_J_env;
        int **Table_J_env_site;

        double *Interaction_J; 

	//N coupling
 	int *Table_N;
 
        int *Table_N_sys;
        int **Table_N_sys_site;

        int *Table_N_env;
        int **Table_N_env_site;

        double *Interaction_N; 

        ~Parameter();

private:

	int TotalSite;

	inline void Parameter_Initial();
	inline void QuantumNumber_Initial();
        inline void Interaction_Table_Initial();

};
