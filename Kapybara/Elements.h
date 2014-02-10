#pragma once
#include<iostream>
#include"Material.h"
#include"IntegratingPoint.h"
#include"Kapybara.h"
using namespace System;
namespace Minilla3D {

	namespace Elements{
		public class iElement
		{
		public:
		    double Volume;
			double refVolume;
			___integratingPoint* iP;
			virtual void computeMetric()=0;
			virtual double computeVolume()=0;
			virtual void memoryMetric()=0;
			virtual void memoryVolume()=0;
			virtual void setupNodes(double* x)=0;
			virtual void setupIndex(int* f)=0;
			virtual void setupNodesFromList(double* x)=0;
			virtual void setRefVolume(double v)=0;
			virtual void computeHessEd()=0;
			virtual void computeGradient()=0;
			virtual void computeGlobalCoord()=0;
			virtual void copyGradient(double* ptr,int DOF)=0;
			virtual void copyGlobalCoord(double* ptr,int num)=0;
			virtual void copyStress(double* ptr,int n)=0;
			virtual void mergeGradient(double* ptr)=0;
			virtual void mergeHessian(ShoNS::Array::SparseDoubleArray ^hess,int dim)=0;
			virtual void computeBaseVectors()=0;
			virtual void computeEigenVectors()=0;
			virtual void copyBaseVectors(double* ptr,int num)=0;
			virtual double copyEigenVectors(double* ptr,double* ptr2,int num)=0;
			virtual double getNode(int i,int j)=0;
			virtual int getIndex(int i)=0;
			virtual void collectStar(array<System::Collections::Generic::List<array<double,1>^>^,1> ^f)=0;
		};

		template<int _nNode,int _elemDim,int _nIntPoint>
		public class element:public iElement
		{
		public:
			___integratingPoint _alias[_nIntPoint];
			double node[_nNode*__DIM];							//êﬂì_ç¿ïW
			int index[_nNode];
			integratingPoint<_nNode,_elemDim> intP[_nIntPoint];	//êœï™ì_
			double gradient[_nNode*__DIM];
			double hess[_nNode*__DIM*_nNode*__DIM];
			double force[_nNode*__DIM];
		
		public:
			element()
			{
				iP=_alias;
				for(int i=0;i<_nNode*__DIM;i++)
				{
					index[i]=i;
				}
			}
			element(int i[_nNode])
			{
				iP=_alias;
				for(int i=0;i<_nNode*__DIM;i++)
				{
					index[i]=i[i];
				}
			}
			virtual double getNode(int i,int j)
			{
				return node[i*__DIM+j];
			}
			virtual int getIndex(int i)
			{
				return index[i];
			}
			virtual void setupNodes(double* x)
			{
				memcpy(node,x,sizeof(node));
			}
			virtual void setupNodesFromList(double* x)
			{
				for(int i=0;i<_nNode;i++)
				{
					for(int j=0;j<__DIM;j++)
					{
						node[i*__DIM+j]=x[index[i]*__DIM+j];
					}
				}
			}
			virtual void setupIndex(int* f)
			{
				memcpy(index,f,sizeof(index));
			}
			virtual double _Volume(){
				return Volume;
			}
			virtual void computeGlobalCoord()
			{
				for(int i=0;i<_nIntPoint;i++)
				{
					intP[i].computeGlobalCoord(node);
				}
			}
			virtual void computeMetric(){
				for(int i=0;i<_nIntPoint;i++)
				{
					intP[i].computeMetric(node);
				}
			}
			virtual void computeBaseVectors()
			{
				for(int i=0;i<_nIntPoint;i++)
				{
					intP[i].computeBaseVectors(node);
				}
			}
			virtual void computeEigenVectors()
			{
				for(int i=0;i<_nIntPoint;i++)
				{
					intP[i].computeEigenVectors();
				}
			}
			virtual double computeVolume(){
				double v=0;
				for(int i=0;i<_nIntPoint;i++)
				{
					v+=intP[i].weight*intP[i].dv;
				}
				this->Volume=v;
				return v;
			}
			virtual void memoryVolume(){
				this->refVolume=this->Volume;
			}
			virtual void setRefVolume(double v){
				this->refVolume=v;
			}
			
			
			virtual void memoryMetric(){
				for(int i=0;i<_nIntPoint;i++)
				{
					intP[i].memoryMetric();
				}
			}
			virtual void computeHessEd()
			{
				int nd=_nNode*__DIM;
				int ee=_elemDim*_elemDim;
				for(int i=0;i<nd;i++)
				{
					for(int j=0;j<nd;j++)
					{
						hess[j+i*nd]=0;
					}
				}
				for(int i=0;i<_nIntPoint;i++)
				{
					double val=intP[i].weight*intP[i].refDv*0.5;
					for(int j=0;j<ee;j++)
					{
						for(int k=0;k<nd;k++)
						{
							for(int l=0;l<nd;l++)
							{
								double D=intP[i].refInvMetric[j];
								double S=intP[i].B[j][l+k*nd];
								hess[l+k*nd]+=val*D*S;
							}
						}
					}
				}
			}

			virtual void computeGradient()
			{
				int nd=_nNode*__DIM;
				int ee=_elemDim*_elemDim;

				for(int i=0;i<nd;i++)
				{
					gradient[i]=0;
				}
				for(int i=0;i<_nIntPoint;i++)
				{
					double val=intP[i].weight*intP[i].dv*0.5;
					for(int j=0;j<ee;j++)
					{
						for(int k=0;k<nd;k++)
						{
							for(int l=0;l<nd;l++)
							{
								double D=intP[i].Cauchy[j];
								double S=intP[i].B[j][l+k*nd];
								double E=node[l];
								gradient[k]+=val*D*S*E;
							}
						}
					}
				}
			}
			virtual void copyGradient(double* ptr,int DOF)
			{
				for(int i=0;i<DOF;i++)
				{
					ptr[i]=gradient[i];
				}
			}
			virtual void copyGlobalCoord(double* ptr,int num)
			{
				memcpy(ptr,intP[num].globalCoord,sizeof(double)*__DIM);
			}
			virtual void copyBaseVectors(double* ptr,int num)
			{
				memcpy(ptr,intP[num].baseVectors,sizeof(double)*__DIM*_elemDim);

			}
			virtual double copyEigenVectors(double* ptr,double* ptr2, int num)
			{
				memcpy(ptr,intP[num].eigenVectors,sizeof(double)*__DIM*_elemDim);
				memcpy(ptr2,intP[num].eigenValues,sizeof(double)*_elemDim);
				return intP[num].weight*intP[num].refDv;
			}
			virtual void copyStress(double* ptr,int n){
				for(int i=0;i<_elemDim*_elemDim;i++)
				{
					ptr[i]=intP[n].Cauchy[i];
				}
			}
			virtual void mergeGradient(double *ptr)
			{
				for(int i=0;i<_nNode;i++)
				{
					for(int j=0;j<__DIM;j++)
					{
						ptr[this->index[i]*__DIM+j]+=this->gradient[i*__DIM+j];
					}
				}
			}
			virtual void mergeHessian(ShoNS::Array::SparseDoubleArray ^hess,int dim)
			{
				int nd=__DIM*_nNode;
				for(int i=0;i<_nNode;i++)
				{
					for(int j=0;j<dim;j++)
					{
						for(int k=0;k<_nNode;k++)
						{
							for(int l=0;l<dim;l++)
							{
								hess[this->index[i]*dim+j,this->index[k]*dim+l]+=this->hess[(i*__DIM+j)*nd+(k*__DIM+l)];
							}
						}
					}
				}
			}
		};
	}
}