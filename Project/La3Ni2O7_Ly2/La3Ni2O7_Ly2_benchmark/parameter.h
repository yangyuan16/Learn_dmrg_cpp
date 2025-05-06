class Parameter {

public:

	Parameter(char &sign, const int &spin, const int &total_h, const int &total_j, const int &number_hole, const int &statenokept, const int &n_u, const int &n_x, const int &n_y, const int &total_n, const double &t_n, const double &t_nn, const double &j_n, const double &j_nn, const double &n_n, const double &n_nn, const int &desiredzeros, const int &seedvalue);

 //------Model and DMRG parameters
 	char Sign;
        int Spin, Total_h, Total_J, Number_hole, StateNoKept, N_u, N_x, N_y, Total_N, untruncated_site, total_site;
        double T_n, T_nn, J_n, J_nn, N_n, N_nn;
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
