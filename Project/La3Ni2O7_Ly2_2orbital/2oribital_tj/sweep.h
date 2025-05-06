class Sweep {
public:

        Sweep(Parameter &para);

private:

	char sign;
        int lsys, lns, step, max_num, num, number, initialstate, site_start, site_end, optimalstate_sweep;
        double ground_state_energy, excited_state_energy, precision;

        inline void Initial(Parameter &para);           //Initialize the system and environment blocks
        inline void Infinite(Parameter &para, int &lsys, const int &lns);       //Infinite iteration process
        inline void Infinite_test(Parameter &para, int &lsys, const int &lns);       //Infinite iteration process

        inline void WriteEnergy(Parameter &para, const double &ground_state_energy, const double &excited_state_energy, const int &lsys, const int &site);
        inline void WriteEnergy_sweep( Parameter &para, const int &site_num, const int &step_num, const double &energy );
        inline void WriteWave_Function( const int &BlockNumber, int *Dim_block, double **WaveFunction );
        inline void WriteWave_Function_Shifted( const int &BlockNumber, int *Dim_block, double **WaveFunction );

//------Finite_Sweep
        inline void Finite_Sweep(Parameter &para, const int &lsys_start, const int &lsys_end, char &D, const int &optimalstate, const int &step, const double &precision);

        inline void Finite_Sweep_LeftToRight(Parameter &para, const int &lsys_start, const int &lsys_end, const int &optimalstate, const int &step, const double &precision);
       	inline void Finite_Sweep_RightToLeft(Parameter &para, const int &lsys_start, const int &lsys_end, const int &optimalstate, const int &step, const double &precision);
        inline void Finite_Sweep_LeftToRight_test(Parameter &para, const int &lsys_start, const int &lsys_end, const int &optimalstate, const int &step, const double &precision);
       	inline void Finite_Sweep_RightToLeft_test(Parameter &para, const int &lsys_start, const int &lsys_end, const int &optimalstate, const int &step, const double &precision);

        inline void Measurement( Parameter &para, const int &step, const double &precision );
	inline void Measurement_LeftToRight( Parameter &para, const double &precision );
	inline void Measurement_LeftToRight_shifted( Parameter &para, const double &precision );
	inline void Measurement_RightToLeft( Parameter &para, const double &precision );

};
