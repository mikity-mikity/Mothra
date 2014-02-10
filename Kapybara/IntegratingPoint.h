#pragma once
#include"Kapybara.h"
namespace Minilla3D
{
	namespace Elements
	{
		public class ___integratingPoint{
			public:
				int elemDim;
				double *weight;
				double *metric;
				double *refMetric;
				double *invMetric;
				double *refInvMetric;
				double *Cauchy;
				double *SPK;
				double *dv;
				double *refDv;
			public:
				void initialize(int _elemDim,double* _weight,double* _metric,double* _refMetric,double* _invMetric,double* _refInvMetric,double* _Chachy,double* _SPK,double *_dv,double *_refDv);
		};
		template<int _nNode,int _elemDim>
		public class integratingPoint{
		public:
			double weight;
			double *N;
			double (*C)[_nNode*__DIM*__DIM];
			double (*B)[_nNode*__DIM*_nNode*__DIM];
			double *localCoord;
			double globalCoord[__DIM];
			double baseVectors[_elemDim][__DIM];
			double eigenVectors[_elemDim][__DIM];
			double eigenValues[_elemDim];
			double metric[_elemDim*_elemDim];
			double strain[_elemDim*_elemDim];
			double refMetric[_elemDim*_elemDim];
			double invMetric[_elemDim*_elemDim];
			double refInvMetric[_elemDim*_elemDim];
			double Cauchy[_elemDim*_elemDim];
			double SPK[_elemDim*_elemDim];
			double dv;
			double refDv;
		public:
			integratingPoint(){
				refDv=1.0;
			}
			void makeAlias(___integratingPoint* iP)
			{
				iP->initialize(_elemDim,&weight,metric,refMetric,invMetric,refInvMetric,Cauchy,SPK,&dv,&refDv);
			}
			void computeGlobalCoord(double x[])
			{
				for(int i=0;i<__DIM;i++)
				{
					double val=0;
					for(int j=0;j<__DIM*_nNode;j++)
					{
						val+=x[j]*this->N[j+i*__DIM*_nNode];
					}
					this->globalCoord[i]=val;
				}
			}
			void computeBaseVectors(double x[])
			{
				for(int i=0;i<_elemDim;i++)
				{
					for(int j=0;j<__DIM;j++)
					{
						double val=0;
						for(int k=0;k<_nNode*__DIM;k++)
						{
							val+=this->C[i][k+j*_nNode*__DIM]*x[k];
						}
						this->baseVectors[i][j]=val;
					}
				}
			}
				
			void computeEigenVectors()
			{
				if(_elemDim==2)
				{

					if(this->dv==0)
					{
						this->eigenValues[0]=0;
						this->eigenValues[1]=0;
						for(int i=0;i<_elemDim;i++)
						{
							for(int k=0;k<__DIM;k++)
							{
								this->eigenVectors[i][k]=0;
							}
						}
					}else
					{
						double a=refInvMetric[0]*metric[0]+refInvMetric[1]*metric[2];
						double b=refInvMetric[0]*metric[1]+refInvMetric[1]*metric[3];
						double c=refInvMetric[2]*metric[0]+refInvMetric[3]*metric[2];
						double d=refInvMetric[2]*metric[1]+refInvMetric[3]*metric[3];
						//1. EigenValues
						double f=a+d;double s=std::sqrt(f*f-4*(a*d-b*c));
						double l1=(f+s)/2;
						double l2=(f-s)/2;

						//2. EigenVectors
						double P11=0,P21=0,P12=0,P22=0;
						if(c!=0)
						{
							P11=l1-d;
							P21=c;

							P12=l2-d;
							P22=c;
						}else if(b!=0)
						{
							P11=b;
							P21=l1-a;

							P12=b;
							P22=l2-a;
						}else
						{
							P11=1;
							P21=0;

							P12=0;
							P22=1;
						}
					
						double norm;
						norm=std::sqrt(metric[0]*P11*P11+2*metric[1]*P11*P21+metric[3]*P21*P21);
						if(norm!=0)
						{
							P11/=norm;
							P21/=norm;
						}
						norm=std::sqrt(metric[0]*P12*P12+2*metric[1]*P12*P22+metric[3]*P22*P22);
						if(norm!=0)
						{
							P12/=norm;
							P22/=norm;
						}
						this->eigenValues[0]=l1;
						this->eigenValues[1]=l2;
						for(int i=0;i<_elemDim;i++)
						{
							for(int k=0;k<__DIM;k++)
							{
								this->eigenVectors[i][k]=0;
							}
						}

						for(int k=0;k<__DIM;k++)
						{
							this->eigenVectors[0][k]+=P11*this->baseVectors[0][k];
							this->eigenVectors[0][k]+=P21*this->baseVectors[1][k];
							this->eigenVectors[1][k]+=P12*this->baseVectors[0][k];
							this->eigenVectors[1][k]+=P22*this->baseVectors[1][k];
						}
					}
				}else
				{
				}
			}
			void computeMetric(double x[])
			{
				for(int i=0;i<_elemDim*_elemDim;i++)
				{
					double val=0;
					for(int k=0;k<_nNode*__DIM;k++)
					{
						for(int l=0;l<_nNode*__DIM;l++)
						{
							val+=0.5*x[k]*B[i][l+k*_nNode*__DIM]*x[l];
						}
					}
					metric[i]=val;
					strain[i]=val-refMetric[i];

				}
				switch(_elemDim)
				{
				case 1:
					_inv1(metric,invMetric);
					dv=std::sqrt(_det1(metric));
					break;
				case 2:
					_inv2(metric,invMetric);
					dv=std::sqrt(_det2(metric));

					break;
				case 3:
					_inv3(metric,invMetric);
					dv=std::sqrt(_det3(metric));
					break;
				default:
					throw gcnew System::NotImplementedException("Sorry, only 1x1, 2x2, and 3x3 matrices can be inversed.");
				}
			}
			void memoryMetric()
			{
				memcpy(refMetric,metric,sizeof(metric));
				memcpy(refInvMetric,invMetric,sizeof(invMetric));
				refDv=dv;
			}
			static void ___CtoB(double (*C)[_nNode*__DIM*__DIM],double (*B)[_nNode*__DIM*_nNode*__DIM])
			{
				for(int i=0;i<_elemDim;i++)
				{
					for(int j=0;j<_elemDim;j++)
					{
						for(int u=0;u<_nNode*__DIM;u++)
						{
							for(int v=0;v<_nNode*__DIM;v++)
							{
								B[j+i*_elemDim][v+u*_nNode*__DIM]=0;
								for(int r=0;r<__DIM;r++)
								{
									//‘ÎÌ‰»
									B[j+i*_elemDim][v+u*_nNode*__DIM]+=C[i][u+r*_nNode*__DIM]*C[j][v+r*_nNode*__DIM];
									B[j+i*_elemDim][v+u*_nNode*__DIM]+=C[j][u+r*_nNode*__DIM]*C[i][v+r*_nNode*__DIM];
								}
							}
						}
					}
				}
			}

		private:
			inline void _inv1(double* from, double* to)
			{
				if(from[0]<=0)
				{
						to[0]=0;
						from[0]=0;
				}else
				{
					to[0]=1/from[0];
				}
			}
			inline void _inv2(double* from, double* to)
			{
				double det=_det2(from);
				if(det<=0)
				{
					to[0]=0;
					to[1]=0;
					to[2]=0;
					to[3]=0;
					from[0]=0;
					from[1]=0;
					from[2]=0;
					from[3]=0;
				}else{
					to[0]=from[3]/det;
					to[3]=from[0]/det;
					to[1]=-from[1]/det;
					to[2]=-from[2]/det;
				}
			}
			inline void _inv3(double* from, double* to)
			{
#define _a0 from[0]
#define _a1 from[3]
#define _a2 from[6]
#define _b0 from[1]
#define _b1 from[4]
#define _b2 from[7]
#define _c0 from[2]
#define _c1 from[5]
#define _c2 from[8]
				double det = _det3(from);
				if(det<=0)
				{
					to[0]=0;
					to[1]=0;
					to[2]=0;
					to[3]=0;
					to[4]=0;
					to[6]=0;
					to[7]=0;
					from[0]=0;
					from[1]=0;
					from[2]=0;
					from[3]=0;
					from[4]=0;
					from[5]=0;
					from[6]=0;
					from[7]=0;

				}else
				{
					to[0]=(_b1*_c2-_b2*_c1)/det;
					to[1]=(_c1*_a2-_c2*_a1)/det;
					to[2]=(_a1*_b2-_a2*_b1)/det;
					to[3]=(_b2*_c0-_b0*_c2)/det;
					to[4]=(_c2*_a0-_c0*_a2)/det;
					to[5]=(_a2*_b0-_a0*_b2)/det;
					to[6]=(_b0*_c1-_b1*_c0)/det;
					to[7]=(_c0*_a1-_c1*_a0)/det;
					to[8]=(_a0*_b1-_a1*_b0)/det;
				}

			}
			inline double _det1(double* m)
			{
				return m[0];
			}
			inline double _det2(double* m)
			{
				return m[0]*m[3]-m[1]*m[2];
			}
			inline double _det3(double* m)
			{
				return m[0] * m[4] * m[8]
				+m[3] * m[7] * m[2]
				+m[6] * m[1] * m[5]
				-m[6] * m[4] * m[2]
				-m[3] * m[1] * m[8]
				-m[0] * m[7] * m[5];
			}
		};
	}
}