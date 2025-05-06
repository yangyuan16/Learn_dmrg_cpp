//This class create and delete spaces for the operators

class SubCommon {

public:

        SubCommon();  //Initialize the system and environment

        SubCommon( const int &i, const int &lsys );  //Read angular space from hard disk


	int IndexNo;
	int TotSiteNo;		//Total site number
	int Block_Num_hole;	//Total hole number

	void CreateSpace_hole_block();
	int * Num_hole_block;
	int * Num_J_hole_block;

	void CreateSpace_J_block();
	int ** Value_J_block;
	int ** Dim_J_block;
	int *** Hole_blockOld;
	int *** J_blockOld;
	int *** Start;

	char ope_sign;  //Number of sites that operators need to be stored
	int operator_number_T;
	int operator_number_J; 
	
	double *** H;
	void Create_H();

	double **** S_Dia;
	double **** S_M_Dia;
	double **** NN;
	void Create_S_Dia();
	void Create_S_M_Dia();
	void Create_NN();

	int * Num_block_T;
	int ** J_block_T_bra;
	int ** J_block_T_ket;
	void CreateSpace_T();

	double **** T;
	void Create_T();

        void CreateSpace3();            //Create space for truncation
        double *** dm_eig;
        double *** dm_wave;

	int * Old_hole;
	int ** Old_J;

  //------Save to disk
        void Print_space( const int &i, const int &lsys );
	void Print_H( const int &i, const int &lsys );
	void Print_S_Dia( const int &i, const int &lsys );
	void Print_S_M_Dia( const int &i, const int &lsys );
	void Print_NN( const int &i, const int &lsys );
	void Print_T( const int &i, const int &lsys );

  //------Delete SubCommon functions
        ~SubCommon();                   //Delete the SubCommon functions!!!

private:

	inline void FreeSpace_T_operator();
	inline void FreeSpace_H();
	inline void FreeSpace_S_Dia();
	inline void FreeSpace_S_M_Dia();
	inline void FreeSpace_NN();
	inline void FreeSpace_hole_block();
	inline void FreeSpace_J_block();
	inline void FreeSpace_T();

};
