#include "NurbsSurface.h"
#define __nIntPoint (_uDim+1)*(_vDim+1)
#define __nNode (_uDim)*(_vDim)
#define __elemDim 2
namespace Minilla3D {
	namespace Elements{

		double* fN(int _i,int _k,int _dim,int dim,double* knot)
		{
		    if (_dim==1)
			{
		        double *F=new double[dim]; 
				for (int i=0;i<dim;i++)
				{
					F[i]=0;
				}
		        if (_k==_i)
				{
					F[dim-1]=1;
				}
		        return F;
			}
		    double *S1=fN(_i,_k,_dim-1,dim,knot);
			double *S2=fN(_i,_k+1,_dim-1,dim,knot);
			double E1=knot[_k+_dim-2]-knot[_k-1];
			double E2=knot[_k+_dim-1]-knot[_k];
			double D1[2]={0,0};
			double D2[2]={0,0};
		    if (E1>0)
			{
				D1[0]=1./E1;
				D1[1]=-knot[_k-1]/E1;
			}
			if (E2>0)
			{
				D2[0]=-1./E2;
				D2[1]=knot[_k+_dim-1]/E2;
			}
		    double *F=new double[dim]; 
			for (int i=0;i<dim;i++)
			{
				F[i]=0;
			}
			for(int i=1;i<dim;i++)
			{
				F[i-1]=F[i-1]+S1[i]*D1[0];
				F[i]=F[i]+S1[i]*D1[1];
				F[i-1]=F[i-1]+S2[i]*D2[0];
				F[i]=F[i]+S2[i]*D2[1];
			}
			delete[] S1;
			delete[] S2;
			return F;
		}
		double* fM(int i,int dim,int ddim,double* knot){
			double *M=new double[dim*dim];
			for(int j=0;j<dim*dim;j++)
			{
				M[j]=0;
			}
		    for(int k=i;k<dim+i;k++)
			{
		        double *D=fN(i+ddim,k,dim,dim,knot);
				for (int n =0;n<dim;n++)
				{
					M[n*dim+k-i]=D[n];
				}
				delete[] D;
			}
			std::cout<<"M"<<std::endl;
			std::cout<<M[0]<<","<<M[1]<<","<<M[2]<<std::endl;
			std::cout<<M[3]<<","<<M[4]<<","<<M[5]<<std::endl;
			std::cout<<M[6]<<","<<M[7]<<","<<M[8]<<std::endl;

			double *S=new double[dim*dim];
			for(int j=0;j<dim*dim;j++)
			{
				S[j]=0;
			}
			for(int  n =1;n<dim+1;n++)
			{
				for (int k=1+n;k< dim+2;k++)
				{
					if (n==dim)
					{
						for (int t =0;t<n-1;t++)
						{
						   S[t*dim+n-1]=0;
						}
						S[(n-1)*dim+n-1]=1;
					}else
					{
						S[(k-2)*dim+n-1]=binominal(dim-n,dim+1-k)*std::pow(i-1,k-1-n);
					}
				}
			}
			std::cout<<"S"<<std::endl;
			std::cout<<S[0]<<","<<S[1]<<","<<S[2]<<std::endl;
			std::cout<<S[3]<<","<<S[4]<<","<<S[5]<<std::endl;
			std::cout<<S[6]<<","<<S[7]<<","<<S[8]<<std::endl;
			double *G=new double[dim*dim];
			for (int j=0;j<dim;j++)
			{
				for(int k=0;k<dim;k++)
				{
					double v=0;
					for(int l=0;l<dim;l++)
					{
						v+=S[j*dim+l]*M[l*dim+k];
					}
					G[j*dim+k]=v;
				}
			}
			delete[] S;
			delete[] M;
			return G;
		}

	}
}
#undef __nIntPoint
#undef __nNode
#undef __elemDim