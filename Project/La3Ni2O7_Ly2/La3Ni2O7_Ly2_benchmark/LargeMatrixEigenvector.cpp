#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "DCppBlas.h"
#include "ZCppBlas.h"
#include "LargeMatrixEigenvector.h"

using namespace std;

int LargeMatrixEigenvector::Davidson(DATATYPE*Ritzvector,const int&dim, const int&N_Min, const int&N_Max,
	   const int&N_Wanted, double*RitzValue, const double&tol, const int&maxiteration)
{
#ifdef DEBUG_DAVIDSON_TIME
   timeval time_begin;
   gettimeofday(&time_begin,0);
#endif 

#ifdef DEBUG_DAVIDSON_DIAG
   totalProcessNum=0.1;
   nonDavidsonNum=0;
#endif

   //assert(dim>N_Max);	
   assert(N_Max>N_Min);
   assert(N_Min>N_Wanted);
   assert(N_Wanted>0);
   assert(maxiteration>50);

/*   
   DATATYPE*V=new DATATYPE[dim*N_Max];
   assert(V);

   DATATYPE*HV=new DATATYPE[dim*N_Max];
   assert(HV);

   DATATYPE*submatrix=new DATATYPE[N_Max*N_Max];
   assert(submatrix);

   DATATYPE*msubmatrix=new DATATYPE[N_Max*N_Max];
   assert(msubmatrix);

   DATATYPE*HRitzVector=new DATATYPE[dim];
   assert(HRitzVector);

   DATATYPE*HII=new DATATYPE[dim];
   assert(HII);
*/

   long int L_dim=dim*N_Max;
   long int LX=L_dim+L_dim+N_Max*N_Max+N_Max*N_Max+dim+dim;
   cout<<"\n LX="<<LX<<endl;

   DATATYPE*V=new DATATYPE[LX];
   assert(V);

   DATATYPE*HV=V+L_dim;

   DATATYPE*submatrix=HV+L_dim;

   DATATYPE*msubmatrix=submatrix+N_Max*N_Max;

   DATATYPE*HRitzVector=msubmatrix+N_Max*N_Max;
  
   DATATYPE*HII=HRitzVector+dim;

   double*msubeigenvalue=new double[N_Max];
   assert(msubeigenvalue);

/*V                    :       Storing the vectors using to construct the reduced Hilbert space.       
 *HV                   :       Storing the vectors H*V
 *submatrix            :       Storing the expression of H in reduced Hilbert space.
 *HII                  :       Storing the digonal element of H 
*/

   getMatrixDiagElement(HII,dim);
   
/*The initial vectors are given outside the programs*/
   
   memcpy(V, Ritzvector, dim*N_Wanted*sizeof(DATATYPE));  

   initializeBaseVector(V,HV,dim,N_Min,submatrix,N_Max,N_Wanted);
   
   int step=N_Min;

   int iteration=N_Min;

   int outiteration=0;

/*Begin to calculate Nth eigenvector.Nth=0,1,2,N_Wanted-1 */
   
   int Nth=0;  

   do
   {
      for(int i=0;i<step;i++)
      {
         for(int j=0;j<step;j++)
         {
            msubmatrix[i*step+j]=submatrix[i*N_Max+j];
	 }
      }

      SSMED(msubmatrix, step, msubeigenvalue,'A');
      
      constructEigenvector(V,dim,msubmatrix,step,Ritzvector+Nth*dim);

      constructEigenvector(HV,dim,msubmatrix,step,HRitzVector);
    
#ifdef DOUBLEDATATYPE
      C_daxpy(dim, -msubeigenvalue[0], Ritzvector+Nth*dim, 1, HRitzVector, 1);
#else
      C_zaxpy(dim, -msubeigenvalue[0], Ritzvector+Nth*dim, 1, HRitzVector, 1);
#endif
     
      double singleresidual=0.0;	      
      for(int j=0;j<dim;j++)
      {
#ifdef DOUBLEDATATYPE

         double modu=HRitzVector[j]*HRitzVector[j];
#elif defined COMPLEXDATATYPE
         double modu=HRitzVector[j].real()*HRitzVector[j].real()+HRitzVector[j].imag()*HRitzVector[j].imag();
#else 
#error "You should define DATATYPE first! See globalconfig.h for more detail."
#endif 	    

         if(singleresidual<modu)singleresidual=modu;
      }

      cout.precision(10);

      cout<<"singleresidual="<<singleresidual<<endl;
	 
      cout<<"eigenvalue="<<msubeigenvalue[0]<<endl;
	
      if(fabs(msubeigenvalue[0])>0.01)singleresidual/=msubeigenvalue[0]*msubeigenvalue[0]; 
      if(singleresidual<tol*tol)
      {
/*Judge if the convergence condition is satisfied, if so, Nth eigenvalue is given to RitzValue[Nth]. */
 
	 RitzValue[Nth]=msubeigenvalue[0];
	 Nth++;

         if(Nth<N_Wanted)
         {
            int NumB=step-1;
            if(NumB>N_Min)NumB=N_Min;

            DATATYPE*p_mv=new DATATYPE[dim*NumB];
            assert(p_mv);

            for(int i=0;i<NumB;i++)
            {
               constructEigenvector(V,dim,msubmatrix+(i+1)*step,step,p_mv+i*dim); 
            }

            //for(int i=0;i<dim*NumB;i++)V[i]=p_mv[i];
            memcpy(V, p_mv, dim*NumB*sizeof(DATATYPE));

            for(int i=0;i<NumB;i++)
            {
               constructEigenvector(HV,dim,msubmatrix+(i+1)*step,step,p_mv+i*dim);
            }

            //for(int i=0;i<dim*NumB;i++)HV[i]=p_mv[i];
            memcpy(HV, p_mv, dim*NumB*sizeof(DATATYPE));

            delete []p_mv;

            step=NumB;

            for(int i=0;i<NumB;i++)
            {
               for(int j=0;j<=i;j++)
               {
                  submatrix[i*N_Max+j]=0;
                  submatrix[j*N_Max+i]=0;
                  if(i==j)submatrix[i*N_Max+j]=msubeigenvalue[i+1];
               }
            } 
         }
      }
      else
      {
	 int ifsubtract=1;

#ifdef DEBUG_DAVIDSON_DIAG
         totalProcessNum++;
#endif

         for(int j=0;j<dim;j++)
	 {
	    if(fabs(msubeigenvalue[0]-HII[j])<0.001)
	    {
#ifdef DEBUG_DAVIDSON_DIAG
               nonDavidsonNum++;
#endif
	       ifsubtract=0;
               break;
	    }
	 }
	    
	 if(ifsubtract==1)
	 {
	    for(int j=0;j<dim;j++)
	    {
	       HRitzVector[j]=HRitzVector[j]/(msubeigenvalue[0]-HII[j]);
	    }
	 }

	 normalize(HRitzVector,dim);
   
         double kesa;
	 
         do
         {
            for(int j=0;j<Nth;j++)
            {
               orthogonalization(Ritzvector+j*dim,HRitzVector,dim);
            }
 
	    for(int j=0;j<step;j++)
	    {	       
	       orthogonalization(V+j*dim,HRitzVector,dim);
            }
               
	    kesa=normalize(HRitzVector,dim);
         }while(kesa<0.10);        
 
	 //for(int j=0;j<dim;j++)V[step*dim+j]=HRitzVector[j];
         memcpy(V+step*dim,HRitzVector,dim*sizeof(DATATYPE));

	 H_V(V+step*dim,HV+step*dim,dim);

	 iteration++;

/*If the iteration exceeds, jump out of the loop.*/
	 
	 if(iteration>maxiteration)break;

	 for(int j=0;j<=step;j++)
	 {
	    DATATYPE sum=innerProduct(V+j*dim,HV+step*dim,dim);
            
            submatrix[j*N_Max+step]=sum;

#ifdef DOUBLEDATATYPE
            submatrix[step*N_Max+j]=sum;
#elif defined COMPLEXDATATYPE
            submatrix[step*N_Max+j]=conj(sum);
#else
#error "You should define DATATYPE first! See globalconfig.h for detail."
#endif
	 }	  
	       
	 step++;

	 if(step==N_Max)
	 {
	    outiteration++;

            //for(int i=0;i<N_Max*N_Max;i++)msubmatrix[i]=submatrix[i];
            memcpy(msubmatrix,submatrix,N_Max*N_Max*sizeof(DATATYPE));

            SSMED(msubmatrix, step, msubeigenvalue, 'A');

	    if(outiteration%4==0)
	    {
               cout<<"outiteration%4==0 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<endl;
  
               constructProjectSubSpace(V,dim,msubmatrix,N_Max, 0, N_Min, N_Min,HV,Ritzvector,Nth);
               
               for(int i=0;i<N_Min;i++)
               {
                  for(int j=0;j<N_Min;j++)
                  {
                     submatrix[i*N_Max+j]=msubmatrix[i*N_Min+j];
                  }
               }

               iteration+=N_Min;
	    }
	    else
	    {
	       DATATYPE*p_mv=new DATATYPE[N_Min*dim];
               assert(p_mv);
                
	       for(int i=0;i<N_Min;i++)
	       {
		  constructEigenvector(V,dim,msubmatrix+i*N_Max,N_Max,p_mv+i*dim);
	       }
               
	       //for(int i=0;i<N_Min*dim;i++)V[i]=p_mv[i];
               memcpy(V, p_mv, N_Min*dim*sizeof(DATATYPE));

	       for(int i=0;i<N_Min;i++)
	       {
		  constructEigenvector(HV,dim,msubmatrix+i*N_Max,N_Max,p_mv+i*dim);
	       }

	       //for(int i=0;i<N_Min*dim;i++)HV[i]=p_mv[i];
               memcpy(HV, p_mv, N_Min*dim*sizeof(DATATYPE));

	       delete []p_mv;
    	    
	       for(int row=0;row<N_Min;row++)
	       {
		  for(int col=0;col<N_Min;col++)
		  {
		     if(row==col)submatrix[row*N_Max+col]=msubeigenvalue[row];
	             else submatrix[row*N_Max+col]=0.0;	      
		  }
	       }	   
	    }
	    
            step=N_Min;	  
	 }
      }
   }
   while(Nth<N_Wanted);

   cout<<"maxiteration in diagonalization is "<<iteration<<endl;
   cout<<"Converged eigenvectors are "<<Nth<<endl;
 
   FILE *e=fopen("e/iteration_step", "a+");
   fprintf(e, "%d\n", iteration);
   fclose(e);

   struct rusage selfusage;
   if(getrusage(RUSAGE_SELF,&selfusage)==0)
   {
      cout<<endl<<"Maximum resident memory used :"<<selfusage.ru_maxrss<<" kb"<<endl<<endl; 
   }

/* 
   delete []HII; 
   delete []HRitzVector;
   delete []msubmatrix;
   delete []submatrix;
   delete []HV;
   delete []V;
*/

   delete []msubeigenvalue;

   delete []V;

/*Here I will renormalize the ritzvectors*/

   for(int i=0; i<Nth; i++)
   {
      normalize(Ritzvector+i*dim,dim);
   }

#ifdef DEBUG_CHECK_ORTH

   for(int i=0;i<Nth-1;i++)
   {
      for(int j=i+1;j<Nth;j++)
      {
         DATATYPE sum=innerProduct(Ritzvector+i*dim,Ritzvector+j*dim,dim);
         assert(fabs(sum)<1.0e-10); 
      }
   }

#endif

#ifdef DEBUG_DAVIDSON_TIME

   timeval time_end;
   gettimeofday(&time_end,0);
   cout<<"Time for whole diagonalization is :"<<time_end.tv_sec-time_begin.tv_sec<<" seconds + "<<
	     time_end.tv_usec-time_begin.tv_usec<<" microseconds"<<endl;
 
#endif  

#ifdef DEBUG_DAVIDSON_DIAG
   cout<<"Ratio of nonDavidson process is : "<<double(nonDavidsonNum)/totalProcessNum<<endl;
#endif

/*Sucessful in finding several wanted eigenvectors! Return the number of converged eigenvectors.*/

   return Nth;
}
	 
int LargeMatrixEigenvector::initializeBaseVector(DATATYPE*V, DATATYPE*HV, const int&dim, const int&minibasenum, 
                                                 DATATYPE*submatrix, const int&maxbasenum, const int&requisitestatenum)
{
   int presentbasenum=0;   //present base number 

   DATATYPE*presentV=V; 
   DATATYPE*presentHV=HV;

   do
   {
      normalize(presentV,dim);

      H_V(presentV,presentHV,dim);

      presentbasenum++; 
      if(presentbasenum==minibasenum)break;
		
      presentV+=dim;
		
      if(presentbasenum>=requisitestatenum)memcpy(presentV,presentHV,dim*sizeof(DATATYPE));

      presentHV+=dim;

/*for the sake of precision and convergence
*I will do :normalize -->orthogonalize -->normalize -->
*orthogonalize -->test sum -->normalize
*/
	//1.normalize the new base

      normalize(presentV,dim);

  //2.Orhtogonalize the new base to other
      for(DATATYPE*previousV=V;previousV<presentV;previousV+=dim)
      {
         orthogonalization(previousV,presentV,dim);
      }

	//3.normalize the new base

      double sum=normalize(presentV,dim);

		//4.test sum
/*In this part,I want to judge if the new base is paralled to others
*which means V[k]=HV[k-1]-b0*V[0]-b1*V[1]-...-b(k-1)*V[k-1].If this
*happens,I will produce another new base by random to take the place
*of V[k]
*/
      while(sum<0.10)
      {
         for(int i=0;i<dim;i++)presentV[i]=drand48();//Produce the new base by random

         normalize(presentV,dim); 

         for(DATATYPE*previousV=V;previousV<presentV;previousV+=dim)orthogonalization(previousV,presentV,dim);

         sum=normalize(presentV,dim);
      }

      for(DATATYPE*previousV=V;previousV<presentV;previousV+=dim)
      {
         orthogonalization(previousV,presentV,dim);
      }
   }
   while(presentbasenum<minibasenum);

   for(int i=0;i<minibasenum;i++)
   {
      for(int j=0;j<=i;j++)
      {
         DATATYPE innpdt=innerProduct(V+i*dim, HV+j*dim, dim);

         submatrix[i*maxbasenum+j]=innpdt;
#ifdef DOUBLEDATATYPE
         submatrix[j*maxbasenum+i]=innpdt;
#elif defined COMPLEXDATATYPE
         submatrix[j*maxbasenum+i]=conj(innpdt);
#else
#error "You should define DATATYPE first! See globalconfig.h for detail."
#endif
      }
   }    
  
   return presentbasenum;
}

void LargeMatrixEigenvector::constructProjectSubSpace(DATATYPE*V, const int&dim, DATATYPE*sm, const int&smdim, 
              const int&shift, const int&AcNum, const int&ReNum, DATATYPE*HV, DATATYPE*Ritzvector, const int&Num)
{
   const double kepa=0.10; //control parameter , reorthogonalization  

   DATATYPE*sv=sm+shift*smdim;

   int InitNum;

   if(AcNum>ReNum)InitNum=ReNum;
   else InitNum=AcNum;

   for(int i=0;i<InitNum;i++)
   {
      constructEigenvector(V,dim,sv+i*smdim,smdim,HV+i*dim);
   } 

   //for(int i=0;i<dim*InitNum;i++)V[i]=HV[i];
   memcpy(V, HV, dim*InitNum*sizeof(DATATYPE));

   if(AcNum<ReNum)
   {
      for(int i=0;i<dim*(ReNum-AcNum);i++)V[dim*AcNum+i]=drand48();

      for(int i=AcNum;i<ReNum;i++)normalize(V+i*dim,dim);
   }

/*orthogonalization 
1. orthogonalized to "Num" RitzVectors
2. orthogonalized to previous V  
*/

   for(int i=0;i<ReNum;i++)
   {
      double resid;
      do
      {
         for(int j=0;j<Num;j++)
         {
            orthogonalization(Ritzvector+j*dim, V+i*dim,dim);
         }

         for(int j=0;j<i;j++)
         {
            orthogonalization(V+j*dim,V+i*dim,dim);
         }

         resid=normalize(V+i*dim,dim);

      }while(resid<kepa);
   }
   
   for(int i=0;i<ReNum;i++)
   {
      H_V(V+i*dim,HV+i*dim,dim);
   }
  
   for(int i=0;i<ReNum;i++)
   {
      for(int j=0;j<=i;j++)
      {
         DATATYPE sum=innerProduct(V+i*dim, HV+j*dim, dim);

         sm[i*ReNum+j]=sum;

#ifdef DOUBLEDATATYPE
         sm[j*ReNum+i]=sum;
#elif defined COMPLEXDATATYPE
         sm[j*ReNum+i]=conj(sum);
#else
#error "You should define DATATYPE first! See globalconfig.h for detail."
#endif
      }
   }
}
