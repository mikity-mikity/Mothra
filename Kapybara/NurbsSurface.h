#pragma once
#include<iostream>
#include"Material.h"
#include"IntegratingPoint.h"
#include"Kapybara.h"
#include"Elements.h"
#define __nIntPoint (_uDim+1)*(_vDim+1)
#define __nNode (_uDim)*(_vDim)
#define __elemDim 2
using namespace System;
namespace Minilla3D {
	namespace Elements{
		double* fM(int i,int dim,int ddim,double* knot);
	    inline int Factorial(int x) 
		{
			if(x==0)return 1;
			if(x==1)return 1;
			if(x==2)return 2;
			if(x==3)return 6;
			if(x==4)return 24;
			if(x==5)return 120;
			int val=1;
			for(int i=2;i<=x;i++)
			{
				val*=i;
			}
			return val;
		}

		inline double binominal(int N,int k)
		{
			return Factorial(N)/Factorial(N-k)/Factorial(k);
		}

		template<int _uDim,int _vDim>
		public class nurbsSurfaceElement:public element<__nNode,__elemDim,__nIntPoint>{
		private:
			static double _cu[__elemDim][_uDim+1];
			static double _pu[__elemDim][_uDim+1];
			static double _weight[__nIntPoint];
			static double _localCoord[__nIntPoint][__elemDim];
			double _N[__nIntPoint][__nNode*__DIM*__DIM];
			double _C[__nIntPoint][__elemDim][__nNode*__DIM*__DIM];
			double _B[__nIntPoint][__elemDim*__elemDim][__nNode*__DIM*__nNode*__DIM];
		public:
			virtual void collectStar(array<System::Collections::Generic::List<array<double,1>^>^,1> ^f)
			{
			}
			nurbsSurfaceElement(int uNum,int vNum,double *uKnot,double *vKnot)
			{
				static int ss[__nIntPoint][__elemDim];		//積分点インデクス
				static int dd[__nNode][__elemDim];			//節点インデクス
				static int dim[2]={_uDim,_vDim};
				static double* hh[__elemDim];
				static double* tt[__elemDim];
				static bool isInitialized=false;
				if(!isInitialized)
				{
					isInitialized=true;
					//多項式生成用
					for(int i=0;i<__elemDim;i++)
					{
						hh[i]=new double[dim[i]];
						tt[i]=new double[dim[i]];
					}
					
					//重み分布初期化
					//座標分布初期化
					if(_uDim==5)
					{
						_cu[0][0]=(-0.9324695142031521)/2.+0.5;
						_cu[0][1]=(-0.6612093864662645)/2.+0.5;
						_cu[0][2]=(-0.2386191860831969)/2.+0.5;
						_cu[0][3]=(0.2386191860831969)/2.+0.5;
						_cu[0][4]=(0.6612093864662645)/2.+0.5;
						_cu[0][5]=(0.9324695142031521)/2.+0.5;
						_pu[0][0]=0.1713244923791704/2.;
						_pu[0][1]=0.3607615730481386/2.;
						_pu[0][2]=0.4679139345726910/2.;
						_pu[0][3]=0.4679139345726910/2.;
						_pu[0][4]=0.3607615730481386/2.;
						_pu[0][5]=0.1713244923791704/2.;
					}
					if(_vDim==5)
					{
						_cu[1][0]=(-0.9324695142031521)/2.+0.5;
						_cu[1][1]=(-0.6612093864662645)/2.+0.5;
						_cu[1][2]=(-0.2386191860831969)/2.+0.5;
						_cu[1][3]=(0.2386191860831969)/2.+0.5;
						_cu[1][4]=(0.6612093864662645)/2.+0.5;
						_cu[1][5]=(0.9324695142031521)/2.+0.5;
						_pu[1][0]=0.1713244923791704/2.;
						_pu[1][1]=0.3607615730481386/2.;
						_pu[1][2]=0.4679139345726910/2.;
						_pu[1][3]=0.4679139345726910/2.;
						_pu[1][4]=0.3607615730481386/2.;
						_pu[1][5]=0.1713244923791704/2.;
					}
					if(_uDim==4)
					{
						_cu[0][0]=(-0.9061798459386640)/2.+0.5;
						_cu[0][1]=(-0.5384693101056831)/2.+0.5;
						_cu[0][2]=(0.0000000000000000)/2.+0.5;
						_cu[0][3]=(0.5384693101056831)/2.+0.5;
						_cu[0][4]=(0.9061798459386640)/2.+0.5;
						_pu[0][0]=0.2369268850561891/2.;
						_pu[0][1]=0.4786286704993665/2.;
						_pu[0][2]=0.5688888888888889/2.;
						_pu[0][3]=0.4786286704993665/2.;
						_pu[0][4]=0.2369268850561891/2.;
					}
					if(_vDim==4)
					{
						_cu[1][0]=(-0.9061798459386640)/2.+0.5;
						_cu[1][1]=(-0.5384693101056831)/2.+0.5;
						_cu[1][2]=(0.0000000000000000)/2.+0.5;
						_cu[1][3]=(0.5384693101056831)/2.+0.5;
						_cu[1][4]=(0.9061798459386640)/2.+0.5;
						_pu[1][0]=0.2369268850561891/2.;
						_pu[1][1]=0.4786286704993665/2.;
						_pu[1][2]=0.5688888888888889/2.;
						_pu[1][3]=0.4786286704993665/2.;
						_pu[1][4]=0.2369268850561891/2.;
					}
					if(_uDim==3)
					{
						_cu[0][0]=(-0.8611363115940526)/2.+0.5;
						_cu[0][1]=(-0.3399810435848563)/2.+0.5;
						_cu[0][2]=(0.3399810435848563)/2.+0.5;
						_cu[0][3]=(0.8611363115940526)/2.+0.5;
						_pu[0][0]=0.3478548451374538/2.;
						_pu[0][1]=0.6521451548625461/2.;
						_pu[0][2]=0.6521451548625461/2.;
						_pu[0][3]=0.3478548451374538/2.;
					}
					if(_vDim==3)
					{
						_cu[1][0]=(-0.8611363115940526)/2.+0.5;
						_cu[1][1]=(-0.3399810435848563)/2.+0.5;
						_cu[1][2]=(0.3399810435848563)/2.+0.5;
						_cu[1][3]=(0.8611363115940526)/2.+0.5;
						_pu[1][0]=0.3478548451374538/2.;
						_pu[1][1]=0.6521451548625461/2.;
						_pu[1][2]=0.6521451548625461/2.;
						_pu[1][3]=0.3478548451374538/2.;
					}
					if(_uDim==2)
					{
						_cu[0][0]=(-0.7745966692414834)/2.+0.5;
						_cu[0][1]=(0.0000000000000000)/2.+0.5;
						_cu[0][2]=(0.7745966692414834)/2.+0.5;
						_pu[0][0]=0.5555555555555556/2.;
						_pu[0][1]=0.8888888888888888/2.;
						_pu[0][2]=0.5555555555555556/2.;
					}
					if(_vDim==2)
					{
						_cu[1][0]=(-0.7745966692414834)/2.+0.5;
						_cu[1][1]=(0.0000000000000000)/2.+0.5;
						_cu[1][2]=(0.7745966692414834)/2.+0.5;
						_pu[1][0]=0.5555555555555556/2.;
						_pu[1][1]=0.8888888888888888/2.;
						_pu[1][2]=0.5555555555555556/2.;
					}
					//積分点インデクス作成
					for(int i=0;i<__elemDim;i++)
					{
						ss[0][i]=0;
					}
					for(int i=1;i<__nIntPoint;i++)
					{
						for(int j=0;j<__elemDim;j++)
						{
							ss[i][j]=ss[i-1][j];
						}
						for(int j=0;j<__elemDim;j++)
						{

							if(ss[i][j]<dim[j])
							{
								ss[i][j]++;
								for(int k=0;k<j;k++)
								{
									ss[i][k]=0;
								}
								break;
							}
						}
					}
					//節点インデクス作成
					for(int i=0;i<__elemDim;i++)
					{
						dd[0][i]=0;
					}
					for(int i=1;i<__nNode;i++)
					{
						for(int j=0;j<__elemDim;j++)
						{
							dd[i][j]=dd[i-1][j];
						}
						for(int j=0;j<__elemDim;j++)
						{
							if(dd[i][j]<dim[j]-1)
							{
								dd[i][j]++;
								for(int k=0;k<j;k++)
								{
									dd[i][k]=0;
								}
								break;
							}
						}
					}
					//積分点係数計算
					for(int i=0;i<__nIntPoint;i++)
					{
						_weight[i]=1.0;
						for(int j=0;j<__elemDim;j++)
						{
							//重み
							_localCoord[i][j]=_cu[j][ss[i][j]];
							_weight[i]*=_pu[j][ss[i][j]];
						}
					}
				}
				//コピー
				for(int i=0;i<__nIntPoint;i++)
				{
					intP[i].weight=_weight[i];
					intP[i].localCoord=_localCoord[i];
					intP[i].makeAlias(&iP[i]);
					intP[i].N=_N[i];
					intP[i].C=_C[i];
					intP[i].B=_B[i];
				    
				}
				double *M[2];
				M[0]=fM(uNum,_uDim,_uDim-1,uKnot);
				M[1]=fM(vNum,_vDim,_vDim-1,vKnot);
				std::cout<<"M[0]"<<std::endl;
			    std::cout<<M[0][0]<<","<<M[0][1]<<","<<M[0][2]<<std::endl;
			    std::cout<<M[0][3]<<","<<M[0][4]<<","<<M[0][5]<<std::endl;
			    std::cout<<M[0][6]<<","<<M[0][7]<<","<<M[0][8]<<std::endl;
				std::cout<<"M[1]"<<std::endl;
				std::cout<<M[1][0]<<","<<M[1][1]<<","<<M[1][2]<<std::endl;
			    std::cout<<M[1][3]<<","<<M[1][4]<<","<<M[1][5]<<std::endl;
			    std::cout<<M[1][6]<<","<<M[1][7]<<","<<M[1][8]<<std::endl;
				//形状関数係数
				for(int i=0;i<__nIntPoint;i++)
				{
					for(int j=0;j<__elemDim;j++)
					{
						double t=intP[i].localCoord[j];
						for(int k=0;k<dim[j];k++)
						{
							hh[j][k]=std::pow(t,(dim[j]-k-1));
						}
						for(int k=0;k<dim[j];k++)
						{
							double val=0;
							for(int l=0;l<dim[j];l++)
							{
								val+=hh[j][l]*M[j][l*dim[j]+k];
							}
							tt[j][k]=val;
						}
					}
					for(int j=0;j<__nNode*__DIM*__DIM;j++)
					{
						intP[i].N[j]=0;
					}
					for(int k=0;k<__nNode;k++)
					{
						//形状関数
						double N=1.0;
						for(int j=0;j<__elemDim;j++)
						{
							N*=tt[j][dd[k][j]];
						}
						for(int j=0;j<__DIM;j++)
						{
							intP[i].N[k*__DIM+j+j*__nNode*__DIM]=N;
						}
					}
					for (int m=0;m<__elemDim;m++)
					{
						for(int j=0;j<__elemDim;j++)
						{
							double t=intP[i].localCoord[j];
							if(j!=m)
							{
								for(int k=0;k<dim[j];k++)
								{
									hh[j][k]=std::pow(t,(dim[j]-k-1));
								}
							}else
							{
								for(int k=0;k<dim[j]-1;k++)
								{
									hh[j][k]=(dim[j]-k-1)*std::pow(t,(dim[j]-k-2));
								}
								hh[j][dim[j]-1]=0;
							}
							for(int k=0;k<dim[j];k++)
							{
								double val=0;
								for(int l=0;l<dim[j];l++)
								{
									val+=hh[j][l]*M[j][l*dim[j]+k];
								}
								tt[j][k]=val;
							}
						}
						for(int j=0;j<__nNode*__DIM*__DIM;j++)
						{
							intP[i].C[m][j]=0;
						}
						for(int k=0;k<__nNode;k++)
						{
							//形状関数
							double C=1.0;
							for(int j=0;j<__elemDim;j++)
							{
								C*=tt[j][dd[k][j]];
							}
							for(int j=0;j<__DIM;j++)
							{
								intP[i].C[m][k*__DIM+j+j*__nNode*__DIM]=C;
							}
						}
					}
					//計量用インデクス更新
					integratingPoint<__nNode,__elemDim>::___CtoB(intP[i].C,intP[i].B);
				}
				for(int i=0;i<__elemDim;i++)
				{
					delete[] M[i];
				}
			}
		};
		
	}
}
#undef __nIntPoint
#undef __nNode
#undef __elemDim