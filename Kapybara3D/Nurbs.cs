using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Minilla3D.Elements
{
    public class nurbsElement:element
    {
        double[,] _cu;
		double[,] _pu;
        public int uDim,vDim;

        double[] fN(int _i, int _k, int _dim, int dim, double[] knot)
		{
		    if (_dim==1)
			{
		        double[] F=new double[dim]; 
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
		    double[] S1=fN(_i,_k,_dim-1,dim,knot);
			double[] S2=fN(_i,_k+1,_dim-1,dim,knot);
			double E1=knot[_k+_dim-2]-knot[_k-1];
			double E2=knot[_k+_dim-1]-knot[_k];
			double[] D1=new double[2]{0,0};
			double[] D2=new double[2]{0,0};
		    if (E1>0)
			{
				D1[0]=1d/E1;
				D1[1]=-knot[_k-1]/E1;
			}
			if (E2>0)
			{
				D2[0]=-1d/E2;
				D2[1]=knot[_k+_dim-1]/E2;
			}
		    double[] F2=new double[dim]; 
			for (int i=0;i<dim;i++)
			{
				F2[i]=0;
			}
			for(int i=1;i<dim;i++)
			{
				F2[i-1]=F2[i-1]+S1[i]*D1[0];
				F2[i]=F2[i]+S1[i]*D1[1];
				F2[i-1]=F2[i-1]+S2[i]*D2[0];
				F2[i]=F2[i]+S2[i]*D2[1];
			}
			return F2;
		}
        double[,] fM(int shift, int dim, int ddim, double[] knot){
			double[,] M=new double[dim,dim];
            for (int i = 0; i < dim; i++)
            {
                for (int j = 0; j < dim; j++)
                {
                    M[i,j] = 0;
                }
            }
		    for(int k=shift;k<dim+shift;k++)
			{
		        double[] D=fN(shift+ddim,k,dim,dim,knot);
				for (int n =0;n<dim;n++)
				{
					M[n,k-shift]=D[n];
				}
			}

			double[,] S=new double[dim,dim];
            for(int i=0;i<dim;i++)
			{
                for(int j=0;j<dim;j++)
			    {
				    S[i,j]=0;
			    }
            }
			for(int  n =1;n<dim+1;n++)
			{
				for (int k=1+n;k<dim+2;k++)
				{
					if (n==dim)
					{
						for (int t =0;t<n-1;t++)
						{
						   S[t,n-1]=0;
						}
						S[(n-1),n-1]=1;
					}else
					{
						S[(k-2),n-1]=binominal(dim-n,dim+1-k)*Math.Pow(shift-1,k-1-n);
					}
				}
			}
			double[,] G=new double[dim,dim];
			for (int j=0;j<dim;j++)
			{
				for(int k=0;k<dim;k++)
				{
					double v=0;
					for(int l=0;l<dim;l++)
					{
						v+=S[j,l]*M[l,k];
					}
					G[j,k]=v;
				}
			}
			return G;
		}
        private static int Factorial(int x)
        {
            if (x == 0) return 1;
            if (x == 1) return 1;
            if (x == 2) return 2;
            if (x == 3) return 6;
            if (x == 4) return 24;
            if (x == 5) return 120;
            int val = 1;
            for (int i = 2; i <= x; i++)
            {
                val *= i;
            }
            return val;
        }

        public static double binominal(int N, int k)
        {
            return Factorial(N) / Factorial(N - k) / Factorial(k);
        }


        public nurbsElement(int _uDim,int _vDim,int _elemDim,int[]_index,int uNum,int vNum,double[] uKnot,double[] vKnot):base(_index,_uDim*_vDim,_elemDim,(_uDim+1)*(_vDim+1)){
            uDim=_uDim;
            vDim=_vDim;
            _cu=new double[elemDim,uDim+1];
            _pu=new double[elemDim,uDim+1];
			int[,] ss=new int[nIntPoint,elemDim];		//Indeces for integrating points
			int[,] dd=new int[nNode,elemDim];		    //Indeces for nodes
            int[] dim=new int[2]{uDim,vDim};
			double[][] hh=new double[elemDim][];
		    double[][] tt=new double[elemDim][];

            //For polynominal
			for(int i=0;i<elemDim;i++)
			{
				hh[i]=new double[dim[i]];
				tt[i]=new double[dim[i]];
			}
					
			//Weight coefficient distribution
			//Coordinates distribution
			if(uDim==5)
			{
				_cu[0,0]=(-0.9324695142031521)*0.5+0.5;
				_cu[0,1]=(-0.6612093864662645)*0.5+0.5;
				_cu[0,2]=(-0.2386191860831969)*0.5+0.5;
				_cu[0,3]=(0.2386191860831969)*0.5+0.5;
				_cu[0,4]=(0.6612093864662645)*0.5+0.5;
				_cu[0,5]=(0.9324695142031521)*0.5+0.5;
				_pu[0,0]=0.1713244923791704*0.5;
				_pu[0,1]=0.3607615730481386*0.5;
				_pu[0,2]=0.4679139345726910*0.5;
				_pu[0,3]=0.4679139345726910*0.5;
				_pu[0,4]=0.3607615730481386*0.5;
				_pu[0,5]=0.1713244923791704*0.5;
			}
			if(vDim==5)
			{
				_cu[1,0]=(-0.9324695142031521)*0.5+0.5;
				_cu[1,1]=(-0.6612093864662645)*0.5+0.5;
				_cu[1,2]=(-0.2386191860831969)*0.5+0.5;
				_cu[1,3]=(0.2386191860831969)*0.5+0.5;
				_cu[1,4]=(0.6612093864662645)*0.5+0.5;
				_cu[1,5]=(0.9324695142031521)*0.5+0.5;
				_pu[1,0]=0.1713244923791704*0.5;
				_pu[1,1]=0.3607615730481386*0.5;
				_pu[1,2]=0.4679139345726910*0.5;
				_pu[1,3]=0.4679139345726910*0.5;
				_pu[1,4]=0.3607615730481386*0.5;
				_pu[1,5]=0.1713244923791704*0.5;
			}
			if(uDim==4)
			{
				_cu[0,0]=(-0.9061798459386640)*0.5+0.5;
				_cu[0,1]=(-0.5384693101056831)*0.5+0.5;
				_cu[0,2]=(0.0000000000000000)*0.5+0.5;
				_cu[0,3]=(0.5384693101056831)*0.5+0.5;
				_cu[0,4]=(0.9061798459386640)*0.5+0.5;
				_pu[0,0]=0.2369268850561891*0.5;
				_pu[0,1]=0.4786286704993665*0.5;
				_pu[0,2]=0.5688888888888889*0.5;
				_pu[0,3]=0.4786286704993665*0.5;
				_pu[0,4]=0.2369268850561891*0.5;
			}
			if(vDim==4)
			{
				_cu[1,0]=(-0.9061798459386640)*0.5+0.5;
				_cu[1,1]=(-0.5384693101056831)*0.5+0.5;
				_cu[1,2]=(0.0000000000000000)*0.5+0.5;
				_cu[1,3]=(0.5384693101056831)*0.5+0.5;
				_cu[1,4]=(0.9061798459386640)*0.5+0.5;
				_pu[1,0]=0.2369268850561891*0.5;
				_pu[1,1]=0.4786286704993665*0.5;
				_pu[1,2]=0.5688888888888889*0.5;
				_pu[1,3]=0.4786286704993665*0.5;
				_pu[1,4]=0.2369268850561891*0.5;
			}
			if(uDim==3)
			{
				_cu[0,0]=(-0.8611363115940526)*0.5+0.5;
				_cu[0,1]=(-0.3399810435848563)*0.5+0.5;
				_cu[0,2]=(0.3399810435848563)*0.5+0.5;
				_cu[0,3]=(0.8611363115940526)*0.5+0.5;
				_pu[0,0]=0.3478548451374538*0.5;
				_pu[0,1]=0.6521451548625461*0.5;
				_pu[0,2]=0.6521451548625461*0.5;
				_pu[0,3]=0.3478548451374538*0.5;
			}
			if(vDim==3)
			{
				_cu[1,0]=(-0.8611363115940526)*0.5+0.5;
				_cu[1,1]=(-0.3399810435848563)*0.5+0.5;
				_cu[1,2]=(0.3399810435848563)*0.5+0.5;
				_cu[1,3]=(0.8611363115940526)*0.5+0.5;
				_pu[1,0]=0.3478548451374538*0.5;
				_pu[1,1]=0.6521451548625461*0.5;
				_pu[1,2]=0.6521451548625461*0.5;
				_pu[1,3]=0.3478548451374538*0.5;
			}
			if(uDim==2)
			{
				_cu[0,0]=(-0.7745966692414834)*0.5+0.5;
				_cu[0,1]=(0.0000000000000000)*0.5+0.5;
				_cu[0,2]=(0.7745966692414834)*0.5+0.5;
				_pu[0,0]=0.5555555555555556*0.5;
				_pu[0,1]=0.8888888888888888*0.5;
				_pu[0,2]=0.5555555555555556*0.5;
			}
			if(vDim==2)
			{
				_cu[1,0]=(-0.7745966692414834)*0.5+0.5;
				_cu[1,1]=(0.0000000000000000)*0.5+0.5;
				_cu[1,2]=(0.7745966692414834)*0.5+0.5;
				_pu[1,0]=0.5555555555555556*0.5;
				_pu[1,1]=0.8888888888888888*0.5;
				_pu[1,2]=0.5555555555555556*0.5;
			}

            //Indeces for integrating points
			for(int i=0;i<elemDim;i++)
			{
				ss[0,i]=0;
			}
			for(int i=1;i<nIntPoint;i++)
			{
				for(int j=0;j<elemDim;j++)
				{
					ss[i,j]=ss[i-1,j];
				}
				for(int j=0;j<elemDim;j++)
				{

					if(ss[i,j]<dim[j])
					{
						ss[i,j]++;
						for(int k=0;k<j;k++)
						{
							ss[i,k]=0;
						}
						break;
					}
				}
			}

			//Indices for nodes
            for(int i=0;i<elemDim;i++)
			{
				dd[0,i]=0;
			}
			for(int i=1;i<nNode;i++)
			{
				for(int j=0;j<elemDim;j++)
				{
					dd[i,j]=dd[i-1,j];
				}
				for(int j=0;j<elemDim;j++)
				{
					if(dd[i,j]<dim[j]-1)
					{
						dd[i,j]++;
						for(int k=0;k<j;k++)
						{
							dd[i,k]=0;
						}
						break;
					}
				}
			}
			//weight coefficients for integrating points
			for(int i=0;i<nIntPoint;i++)
			{
				intP[i].weight=1.0;
				for(int j=0;j<elemDim;j++)
				{
					intP[i].localCoord[j]=_cu[j,ss[i,j]];
					intP[i].weight*=_pu[j,ss[i,j]];
				}
			}
		
		    double[][,] M=new double[2][,];
		    M[0]=fM(uNum,_uDim,_uDim-1,uKnot);
		    M[1]=fM(vNum,_vDim,_vDim-1,vKnot);
			//Shape functions  [N] (for global coordinate)
			for(int i=0;i<nIntPoint;i++)
			{
				for(int j=0;j<elemDim;j++)
				{
					double t=intP[i].localCoord[j];
					for(int k=0;k<dim[j];k++)
					{
						hh[j][k]=Math.Pow(t,(dim[j]-k-1));
					}
					for(int k=0;k<dim[j];k++)
					{
						double val=0;
						for(int l=0;l<dim[j];l++)
						{
							val+=hh[j][l]*M[j][l,k];
						}
						tt[j][k]=val;
					}
				}
				for(int j=0;j<__DIM;j++)
				{
                    for(int k=0;k<nDV;k++)
                    {
					    intP[i].N[j,k]=0;
                    }
				}
				for(int k=0;k<nNode;k++)
				{
					//Shape functinos
					double N=1.0;
					for(int j=0;j<elemDim;j++)
					{
						N*=tt[j][dd[k,j]];
					}
					for(int j=0;j<__DIM;j++)
					{
						intP[i].N[j,k*__DIM+j]=N;
					}
				}

                //Create [C]  (for base vectors)
				for (int m=0;m<elemDim;m++)
				{
					for(int j=0;j<elemDim;j++)
					{
						double t=intP[i].localCoord[j];
						if(j!=m)
						{
							for(int k=0;k<dim[j];k++)
							{
								hh[j][k]=Math.Pow(t,(dim[j]-k-1));
							}
						}else
						{
							for(int k=0;k<dim[j]-1;k++)
							{
								hh[j][k]=(dim[j]-k-1)*Math.Pow(t,(dim[j]-k-2));
							}
							hh[j][dim[j]-1]=0;
						}
						for(int k=0;k<dim[j];k++)
						{
							double val=0;
							for(int l=0;l<dim[j];l++)
							{
								val+=hh[j][l]*M[j][l,k];
							}
							tt[j][k]=val;
						}
					}
                    for(int jj=0;jj<__DIM;jj++)
                    {
					    for(int j=0;j<nDV;j++)
					    {
						    intP[i].C[m,jj,j]=0;
					    }
                    }
					for(int k=0;k<nNode;k++)
					{
						//[C]
						double C=1.0;
						for(int j=0;j<elemDim;j++)
						{
							C*=tt[j][dd[k,j]];
						}
						for(int j=0;j<__DIM;j++)
						{
							intP[i].C[m,j,k*__DIM+j]=C;
						}
					}
                }
				//Create [B]  (for metric)
                intP[i].CtoB();

                //Create [D] (for second derivative)
                for (int m = 0; m < elemDim; m++)
                {
                    for (int n = 0; n < elemDim; n++)
                    {
                        for (int j = 0; j < elemDim; j++)
                        {
                            double t = intP[i].localCoord[j];
                            if (j != m&&j!=n)
                            {
                                for (int k = 0; k < dim[j]; k++)
                                {
                                    hh[j][k] = Math.Pow(t, (dim[j] - k - 1));
                                }
                            }
                            if((j!=m&&j==n)||(j==m&&j!=n))
                            {
                                for (int k = 0; k < dim[j] - 1; k++)
                                {
                                    hh[j][k] = (dim[j] - k - 1) * Math.Pow(t, (dim[j] - k - 2));
                                }
                                hh[j][dim[j] - 1] = 0;
                            }
                            if (j == m && j == n)
                            {
                                for (int k = 0; k < dim[j] - 1; k++)
                                {
                                    hh[j][k] = (dim[j] - k - 1) *(dim[j] - k - 2) * Math.Pow(t, (dim[j] - k - 3));
                                }
                                hh[j][dim[j] - 1] = 0;
                                hh[j][dim[j] - 2] = 0;
                            }

                            for (int k = 0; k < dim[j]; k++)
                            {
                                double val = 0;
                                for (int l = 0; l < dim[j]; l++)
                                {
                                    val += hh[j][l] * M[j][l, k];
                                }
                                tt[j][k] = val;
                            }
                        }
                        for (int jj = 0; jj < __DIM; jj++)
                        {
                            for (int j = 0; j < nDV; j++)
                            {
                                intP[i].D[m, n,jj, j] = 0;
                            }
                        }
                        for (int k = 0; k < nNode; k++)
                        {
                            //[D]
                            double D = 1.0;
                            for (int j = 0; j < elemDim; j++)
                            {
                                D *= tt[j][dd[k, j]];
                            }
                            for (int j = 0; j < __DIM; j++)
                            {
                                intP[i].D[m,n, j, k * __DIM + j] = D;
                            }
                        }
                    }
                }
            }
        }
    }
}
