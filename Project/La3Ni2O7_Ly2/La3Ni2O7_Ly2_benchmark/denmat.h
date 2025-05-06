class DenMat {
public:
        DenMat(const int &f, Parameter &para, Super &sup, Sub &sysnew);
        ~DenMat();

private:
	
	int QuantumNumber_hole,  QuantumNumber_J;
	
    //--parameters for BLAS subroutine
        char jobz, uplo, trans_N, trans_T;
        int BlockNumber, Dim_den, n, lda, lwork, info, x, y;
        double beta, factor, den_sum;

        int *H_sysnew, *H_envnew, *J_sysnew, *J_envnew, *Dim_block;

        double *a;
        double *work;
        double *w;
 
        double **wavefunc_new;

        inline void CreateSpace(Parameter &para, Super &sup);

        inline void DenMat_Jacobi_sysnew(Super &sup, Sub &sysnew);
        inline void DenMat_Jacobi_envnew(Super &sup, Sub &envnew);

        inline void Find_aa_sysnew(Super &sup, int &jn, int &js);
        inline void Find_aa_envnew(Super &sup, int &jn, int &je);

        inline void Initial_Lapack(int &Dim);

        inline void FreeSpace();

};
