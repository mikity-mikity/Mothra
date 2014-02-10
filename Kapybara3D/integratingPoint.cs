﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Minilla3D.Elements
{    
    public class integratingPoint
    {
        const int __DIM=3;
        int nNode,elemDim;
        int nDV;
        public double weight;
		public double[,] N;
        public double[,,]C;
        public double[, , ,] B;
        public double[, , ,] D;
        public double[] localCoord;
        public double[] globalCoord;
        public double[,] baseVectors;
        public double[, ,] Gamma;
        public double[,] F;
        public double[,] f;
        public double[,] hessUV;
        public double[,] hessXY;
        public double[][] eigenVectors;
        public double[] eigenValues;
		double[,] metric;
		double[,] refMetric;
		double[,] invMetric;
		double[,] refInvMetric;
        public double[,] SPK
        {
            protected set;
            get;
        }
        public double[,] Cauchy
        {
            protected set;
            get;
        }
        public double refDv
        {
            protected set;
            get;
        }
        public double dv
        {
            protected set;
            get;
        }
        public integratingPoint(int _nNode, int _elemDim)
        {
            nNode = _nNode;
            elemDim = _elemDim;
            nDV=nNode*__DIM;
            N=new double[__DIM,nDV];                        //Global coordinate *coefficient*
            C=new double[elemDim,__DIM,nDV];                //Base vectors *coefficient*
            B=new double[elemDim,elemDim,nDV,nDV];          //Metric *coefficient*
            D = new double[elemDim, elemDim, __DIM, nDV];   //Hessian coefficient
            localCoord=new double[elemDim];                 //Local coordinate
            globalCoord=new double[__DIM];                  //Global coordinate
            baseVectors = new double[elemDim, __DIM];       //covariant base vectors
            Gamma = new double[2, 2, 2];  //Connection coefficient
            F = new double[2, 2];                    //Transform matrix (u,v)->(x,y)
            f = new double[2, 3];                    //Transform matrix (x,y)->(u,v)
            hessUV = new double[elemDim, elemDim];          //hessian of airy function with respect to uv
            hessXY = new double[elemDim, elemDim];          //hessian of airy function with respect to xy
            eigenVectors=new double[2][]{new double[3],new double[3]};
            eigenValues=new double[elemDim];
            
            metric=new double[elemDim,elemDim];
            refMetric=new double[elemDim,elemDim];
            invMetric=new double[elemDim,elemDim];
            refInvMetric=new double[elemDim,elemDim];
            SPK = new double[elemDim, elemDim];
            Cauchy = new double[elemDim, elemDim];
            weight = 1.0;
            refDv = 1.0;
            dv=1.0;
        }
        public void computeAiryFunction(double[] x)
        {
            for (int n = 0; n < 2; n++)
            {
                for (int m = 0; m < 2; m++)
                {
                    double val = 0;
                    for (int i = 0; i < nDV; i++)
                    {
                        val += D[n, m, 2, i] * x[i];
                    }
                    hessUV[n, m]=val;
                }
            }
            //covariant base vectors
            for (int n = 0; n < 2; n++)
            {
                double fx=0, fy=0, fphi=0;
                for (int i = 0; i < nDV; i++)
                {
                    fx+= C[n, 0, i] * x[i];
                    fy+= C[n, 1, i] * x[i];
                    fphi+= C[n, 2, i] * x[i];
                }
                f[n, 0] = fx;
                f[n, 1] = fy;
                f[n, 2] = fphi; //first order derivative of phi
            }
            for (int n = 0; n < 2; n++)
            {
                for (int m = 0; m < 2; m++)
                {
                    metric[n, m] = f[n, 0] * f[m, 0] + f[n, 1] * f[m, 1];
                }
            }
            _inv2(metric, invMetric);
            dv = Math.Sqrt(_det2(metric));
            refDv = dv;
            double dv2 = dv * dv;
            //contravatiant base vectors
            for (int n = 0; n < 2; n++)
            {
                F[n, 0] = f[0, 0] * invMetric[0, n] + f[1, 0] * invMetric[1, n];
                F[n, 1] = f[0, 1] * invMetric[0, n] + f[1, 1] * invMetric[1, n];
            }
            
            //Connection coefficients
            for (int n = 0; n < 2; n++)
            {
                for (int m = 0; m < 2; m++)
                {
                    double gx = 0, gy = 0;
                    for (int i = 0; i < nDV; i++)
                    {
                        gx += D[n, m, 0, i] * x[i];
                        gy += D[n, m, 1, i] * x[i];
                    }
                    for(int k=0;k<elemDim;k++)
                    {
                        Gamma[n, m, k] = gx * F[k, 0] + gy * F[k, 1];
                    }
                }
            }
            //Overrite hessuv with connection coefficients
            for (int m = 0; m < 2; m++)
            {
                for (int n = 0; n < 2; n++)
                {
                    double val = 0;
                    for (int k = 0; k < 2; k++)
                    {
                        val += Gamma[n, m, k] * f[k, 2];
                    }
                    hessUV[n, m] -= val;
                }
            }

            //Hodge star on curvelinear coordinate
            double temp = hessUV[1, 1];
            SPK[1, 1] = hessUV[0, 0] / dv2;
            SPK[0, 0] = temp / dv2;
            SPK[0, 1] = -hessUV[0, 1] / dv2;
            SPK[1, 0] = -hessUV[1, 0] / dv2;
            
            
            /*//Coordinate transformation
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    double val = 0;
                    for (int n = 0; n < 2; n++)
                    {
                        for (int m = 0; m < 2; m++)
                        {
                            val += hessUV[n, m] * F[n, i] * F[m, j];
                        }
                    }
                    hessXY[i, j] = val;
                }
            }*/
            //Hodge star?????
            /*double temp = hessXY[1, 1];
            hessXY[1, 1] = hessXY[0, 0];
            hessXY[0, 0] = temp;
            hessXY[0, 1] = -hessXY[0, 1];
            hessXY[1, 0] = -hessXY[1, 0];
            */
            
            
            /*//Coordinate transformation again! ->SPK
            for (int n = 0; n < 2; n++)
            {
                for (int m = 0; m < 2; m++)
                {
                    double val = 0;
                    for (int i = 0; i < 2; i++)
                    {
                        for (int j = 0; j < 2; j++)
                        {
                            val += hessXY[i, j] * F[n, i] * F[m, j];
                        }
                    }
                    SPK[n, m] = val;
                }
            }*/

        }
        public void computeGlobalCoord(double[] x)
		{
			for(int i=0;i<__DIM;i++)
			{
				double val=0;
				for(int j=0;j<nDV;j++)
				{
					val+=this.N[i,j]*x[j];
				}
				this.globalCoord[i]=val;
			}
		}
		public void computeBaseVectors(double[] x)
		{
			for(int i=0;i<elemDim;i++)
			{
				for(int j=0;j<__DIM;j++)
				{
					double val=0;
					for(int k=0;k<nDV;k++)
					{
						val+=this.C[i,j,k]*x[k];
					}
					this.baseVectors[i,j]=val;
				}
			}
		}
				
		//compute Eigen Vectors of SPK
        public void computeEigenVectors()
		{
			if(elemDim!=2)return;

			if(this.dv==0)
			{
				this.eigenValues[0]=0;
				this.eigenValues[1]=0;
				for(int i=0;i<elemDim;i++)
				{
					for(int k=0;k<__DIM;k++)
					{
						this.eigenVectors[i][k]=0;
					}
				}
			}else
			{
				double a=SPK[0,0]*metric[0,0]+SPK[0,1]*metric[1,0];
				double c=SPK[0,0]*metric[0,1]+SPK[0,1]*metric[1,1];
				double b=SPK[1,0]*metric[0,0]+SPK[1,1]*metric[1,0];
				double d=SPK[1,0]*metric[0,1]+SPK[1,1]*metric[1,1];

                //1. EigenValues
				double f=a+d;double s=Math.Sqrt(f*f-4*(a*d-b*c));
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
				norm=Math.Sqrt(invMetric[0,0]*P11*P11+2*invMetric[0,1]*P11*P21+invMetric[1,1]*P21*P21);
				if(norm!=0)
				{
					P11/=norm;
					P21/=norm;
				}
				norm=Math.Sqrt(invMetric[0,0]*P12*P12+2*invMetric[0,1]*P12*P22+invMetric[1,1]*P22*P22);
				if(norm!=0)
				{
					P12/=norm;
					P22/=norm;
				}
				this.eigenValues[0]=l1;
				this.eigenValues[1]=l2;

				for(int i=0;i<elemDim;i++)
				{
					for(int k=0;k<3;k++)
					{
						this.eigenVectors[i][k]=0;
					}
				}

				for(int k=0;k<2;k++)
				{
					this.eigenVectors[0][k]+=P11*this.F[0,k];
					this.eigenVectors[0][k]+=P21*this.F[1,k];
					this.eigenVectors[1][k]+=P12*this.F[0,k];
					this.eigenVectors[1][k]+=P22*this.F[1,k];
				}
			}
		}
		public void computeMetric(double[] x)
		{
			for(int i=0;i<elemDim;i++)
			{
                for(int j=0;j<elemDim;j++)
				{
                    double val=0;
				    for(int k=0;k<nDV;k++)
				    {
					    for(int l=0;l<nDV;l++)
					    {
						    val+=0.5*x[k]*B[i,j,k,l]*x[l];
					    }
				    }
				    metric[i,j]=val;
				    //strain[i,j]=val-refMetric[i];
                }
			}
			switch(elemDim)
			{
			case 1:
				_inv1(metric,invMetric);
				dv=Math.Sqrt(_det1(metric));
				break;
			case 2:
				_inv2(metric,invMetric);
				dv=Math.Sqrt(_det2(metric));

				break;
			case 3:
				_inv3(metric,invMetric);
				dv=Math.Sqrt(_det3(metric));
				break;
			default:
				throw new System.NotImplementedException("Sorry, only 1x1, 2x2, and 3x3 matrices can be inversed.");
			}
		}
		public void memoryMetric()
		{
            Array.Copy(metric,refMetric,elemDim*elemDim);
            Array.Copy(refMetric,refInvMetric,elemDim*elemDim);
			refDv=dv;
		}
	    public void CtoB()
		{
			for(int i=0;i<elemDim;i++)
			{
				for(int j=0;j<elemDim;j++)
				{
					for(int u=0;u<nDV;u++)
					{
						for(int v=0;v<nDV;v++)
						{
							B[i,j,u,v]=0;
							for(int r=0;r<__DIM;r++)
							{
								//対称化
								B[i,j,u,v]+=C[i,r,u]*C[j,r,v];
								B[i,j,u,v]+=C[j,r,u]*C[i,r,v];
							}
						}
					}
				}
			}
		}
		private void _inv1(double[,] from, double[,] to)
		{
			if(from[0,0]<=0)
			{
					to[0,0]=0;
					from[0,0]=0;
			}else
			{
				to[0,0]=1/from[0,0];
			}
		}
		private void _inv2(double[,] from, double[,] to)
		{
			double det=_det2(from);
			if(det<=0)
			{
				to[0,0]=0;
				to[1,0]=0;
				to[0,1]=0;
				to[1,1]=0;
				from[0,0]=0;
				from[1,0]=0;
				from[0,1]=0;
				from[1,1]=0;
			}else{
				to[0,0]=from[1,1]/det;
				to[1,1]=from[0,0]/det;
				to[1,0]=-from[1,0]/det;
				to[0,1]=-from[0,1]/det;
			}
		}
	    private void _inv3(double[,] from, double[,] to)
		{
			double det = _det3(from);
			if(det<=0)
			{
                to[0, 0] = 0;
                to[0, 1] = 0;
                to[0, 2] = 0;
                to[1, 0] = 0;
                to[1, 1] = 0;
                to[1, 2] = 0;
                to[1, 3] = 0;
                to[2, 0] = 0;
                to[2, 1] = 0;
                to[2, 2] = 0;
                from[0, 0] = 0;
                from[0, 1] = 0;
                from[0, 2] = 0;
                from[1, 0] = 0;
                from[1, 1] = 0;
                from[1, 2] = 0;
                from[1, 3] = 0;
                from[2, 0] = 0;
                from[2, 1] = 0;
                from[2, 2] = 0;
			}else
			{
				to[0,0]=(from[1,1]*from[2,2]-from[1,2]*from[2,1])/det;
				to[1,0]=(from[2,1]*from[0,2]-from[2,2]*from[0,1])/det;
				to[2,0]=(from[0,1]*from[1,2]-from[0,2]*from[1,1])/det;
				to[0,1]=(from[1,2]*from[2,0]-from[1,0]*from[2,2])/det;
				to[1,1]=(from[2,2]*from[0,0]-from[2,0]*from[0,2])/det;
				to[2,1]=(from[0,2]*from[1,0]-from[0,0]*from[1,2])/det;
				to[0,2]=(from[1,0]*from[2,1]-from[1,1]*from[2,0])/det;
				to[1,2]=(from[2,0]*from[0,1]-from[2,1]*from[0,0])/det;
				to[2,2]=(from[0,0]*from[1,1]-from[0,1]*from[1,0])/det;
			}
		}
		private double _det1(double[,] m)
		{
			return m[0,0];
		}
		private double _det2(double[,] m)
		{
			return m[0,0]*m[1,1]-m[0,1]*m[1,0];
		}
		private double _det3(double[,] m)
		{
			return m[0,0] * m[1,1] * m[2,2]
			+m[1,0] * m[2,1] * m[0,2]
			+m[2,0] * m[0,1] * m[1,2]
			-m[2,0] * m[1,1] * m[0,2]
			-m[1,0] * m[0,1] * m[2,2]
			-m[0,0] * m[2,1] * m[1,2];
		}

    }
}