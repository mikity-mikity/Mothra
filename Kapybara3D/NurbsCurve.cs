using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Minilla3D.Elements
{
    public class nurbsCurve:element
    {
        double[] _cu;
		double[] _pu;
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
        bool exceptional()
        {
            if (this.typeOfBorder == (border.Left | border.Top))
            {
                return true;
            }
            if (this.typeOfBorder == (border.Left | border.Bottom))
            {
                return true;
            }
            if (this.typeOfBorder == (border.Right | border.Top))
            {
                return true;
            }
            if (this.typeOfBorder == (border.Right | border.Bottom))
            {
                return true;
            }
            return false;
        }
        public void setPlane(double a, double b, double c, double d)
        {
            foreach (var p in intP)
            {
                p.refIntP.setPlane(a, b, c, d);
            }
        }
        public int numberOfConstraintConditions()
        {
            return nIntPoint;
        }
        public int mergeResidual(ShoNS.Array.DoubleArray residual, int i)
        {
            if (exceptional()) return i;
            for (int j = 0; j < nBIntPoint; j++)
            {
                for (int k = 0; k < 2; k++)
                {
                    double resid = bIntP[j].getResidualOfBoundaryCondition(node,k);
                    residual[i + j * 2 + k] = resid;
                }
            }
            return i + nBIntPoint * 2;
        }
        public int mergeJacobian(ShoNS.Array.SparseDoubleArray jacobian, int i)
        {
            if (exceptional()) return i;
            for (int j = 0; j < nBIntPoint; j++)
            {
                for (int k = 0; k < 2; k++)
                {
                    var grad = bIntP[j].getGradientOfBoundaryCondition(k);
                    for(int f=0;f<nNode;f++)
                    {
                        jacobian[i + j*2 + k, index[f]] = grad[f];
                    }
                }
            }
            return i + nBIntPoint * 2;
        }
        public int mergeJacobianOfPosition(ShoNS.Array.SparseDoubleArray jacob, int num)
        {
            for (int i = 0; i < nIntPoint; i++)
            {
                var grad = intP[i].N;
                for (int j = 0; j < nNode; j++)
                {
                    jacob[num + i, index[j]] = grad[2, j * 3 + 2];
                }
            }
            return num + nIntPoint;
        }

        double __cu(int n, int dim)
        {
            switch (dim)
            {
                case 6:
                    switch (n)
                    {
                        case 0:
                            return (-0.9491079123427585245261897) * 0.5 + 0.5;
                        case 1:
                            return (-0.7415311855993944398638648) * 0.5 + 0.5;
                        case 2:
                            return (-0.4058451513773971669066064) * 0.5 + 0.5;
                        case 3:
                            return 0.5;
                        case 4:
                            return (0.4058451513773971669066064) * 0.5 + 0.5;
                        case 5:
                            return (0.7415311855993944398638648) * 0.5 + 0.5;
                        case 6:
                            return (0.9491079123427585245261897) * 0.5 + 0.5;
                        default:
                            return 0;
                    }
                case 5:
                    switch (n)
                    {
                        case 0:
                            return (-0.9324695142031521) * 0.5 + 0.5;
                        case 1:
                            return (-0.6612093864662645) * 0.5 + 0.5;
                        case 2:
                            return (-0.2386191860831969) * 0.5 + 0.5;
                        case 3:
                            return (0.2386191860831969) * 0.5 + 0.5;
                        case 4:
                            return (0.6612093864662645) * 0.5 + 0.5;
                        case 5:
                            return (0.9324695142031521) * 0.5 + 0.5;
                        default:
                            return 0;
                    }
                case 4:
                    switch (n)
                    {
                        case 0:
                            return (-0.9061798459386640)*0.5+0.5;
                        case 1:
                            return (-0.5384693101056831)*0.5+0.5;
                        case 2:
                            return (0.0000000000000000)*0.5+0.5;
                        case 3:
                            return (0.5384693101056831)*0.5+0.5;
                        case 4:
                            return (0.9061798459386640)*0.5+0.5;
				        default:
                            return 0;
                    }
                case 3:
                    switch (n)
                    {
                        case 0:
                            return (-0.8611363115940526)*0.5+0.5;
                        case 1:
                            return (-0.3399810435848563)*0.5+0.5;
                        case 2:
                            return (0.3399810435848563)*0.5+0.5;
                        case 3:
                            return (0.8611363115940526)*0.5+0.5;
                        default:
                            return 0;
                    }
                case 2:
                    switch (n)
                    {
                        case 0:
                            return (-0.7745966692414834) * 0.5 + 0.5;
                        case 1:
                            return (0.0000000000000000) * 0.5 + 0.5;
                        case 2:
                            return (0.7745966692414834) * 0.5 + 0.5;
                        default:
                            return 0;
                    }
                default:
                    return 0;

            }
        }
        double __pu(int n, int dim)
        {
            switch (dim)
            {
                case 6:
                    switch (n)
                    {
                        case 0:
                            return 0.1294849661688696932706114 * 0.5;
                        case 1:
                            return 0.2797053914892766679014678 * 0.5;;
                        case 2:
                            return 0.3818300505051189449503698 * 0.5;;
                        case 3:
                            return 0.4179591836734693877551020 * 0.5;;
                        case 4:
                            return 0.3818300505051189449503698 * 0.5;;
                        case 5:
                            return 0.2797053914892766679014678 * 0.5;;
                        case 6:
                            return 0.1294849661688696932706114 * 0.5;;
                        default:
                            return 0;
                    }
                case 5:
                    switch (n)
                    {
                        case 0:
                            return 0.1713244923791704 * 0.5;
                        case 1:
                            return 0.3607615730481386 * 0.5;
                        case 2:
                            return 0.4679139345726910 * 0.5;
                        case 3:
                            return 0.4679139345726910 * 0.5;
                        case 4:
                            return 0.3607615730481386 * 0.5;
                        case 5:
                            return 0.1713244923791704 * 0.5;
                        default:
                            return 0;
                    }
                case 4:
                    switch (n)
                    {
                        case 0:
                            return 0.2369268850561891 * 0.5;
                        case 1:
                            return 0.4786286704993665 * 0.5;
                        case 2:
                            return 0.5688888888888889 * 0.5;
                        case 3:
                            return 0.4786286704993665 * 0.5;
                        case 4:
                            return 0.2369268850561891 * 0.5;
                        default:
                            return 0;
                    }
                case 3:
                    switch (n)
                    {
                        case 0:
                            return 0.3478548451374538*0.5;
                        case 1:
                            return 0.6521451548625461*0.5;
                        case 2:
                            return 0.6521451548625461*0.5;
                        case 3:
                            return 0.3478548451374538 * 0.5;
                        default:
                            return 0;
                    }
                case 2:
                    switch (n)
                    {
                        case 0:
                            return 0.5555555555555556*0.5;
                        case 1:
                            return 0.8888888888888888*0.5;
                        case 2:
                            return 0.5555555555555556 * 0.5;
                        default:
                            return 0;
                    }
                default:
                    return 0;
            }
        }
        public nurbsCurve(int _uDim, int[] _index, int uNum, double[] uKnot, border _border = border.None)
            : base(_index, _uDim, 1, (_uDim + 1))
        {
            uDim = _uDim;
            _cu = new double[uDim + 1];
            _pu = new double[uDim + 1];
            int[] ss = new int[nIntPoint];		//Indeces for integrating points
			int[] dd=new int[nNode];		    //Indeces for nodes
            int dim=uDim;
			double[] hh=new double[uDim];
		    double[] tt=new double[uDim];

					
			//Weight coefficient distribution
			//Coordinates distribution
            for (int i = 0; i < uDim+1; i++)
            {
                _cu[i] = __cu(i, uDim);
                _pu[i] = __pu(i, uDim);
            }
            

            //Indeces for integrating points
            ss[0] = 0;
            for (int i = 1; i < nIntPoint; i++)
			{
				ss[i]=ss[i-1];
				if(ss[i]<dim)
				{
					ss[i]++;
				}
			}

			//Indices for nodes
		    dd[0]=0;
			for(int i=1;i<nNode;i++)
			{
			    dd[i]=dd[i-1];
				if(dd[i]<dim-1)
				{
					dd[i]++;
				}
			}
			//weight coefficients and loacl coordinates for integrating points
			for(int i=0;i<nIntPoint;i++)
			{
				intP[i].weight=1.0;
				intP[i].localCoord[0]=_cu[ss[i]];
				intP[i].weight*=_pu[ss[i]];
			}
		    double[,] M=null;
		    M=fM(uNum,_uDim,_uDim-1,uKnot);
			//Shape functions  [N] (for global coordinate)
			for(int i=0;i<nIntPoint;i++)
			{
				double t=intP[i].localCoord[0];
				for(int k=0;k<dim;k++)
				{
					hh[k]=Math.Pow(t,(dim-k-1));
				}
				for(int k=0;k<dim;k++)
				{
					double val=0;
					for(int l=0;l<dim;l++)
					{
						val+=hh[l]*M[l,k];
					}
					tt[k]=val;
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
					N*=tt[dd[k]];
					for(int j=0;j<__DIM;j++)
					{
						intP[i].N[j,k*__DIM+j]=N;
					}
				}

                //Create [C]  (for base vectors)
				t=intP[i].localCoord[0];
				{
					for(int k=0;k<dim-1;k++)
					{
						hh[k]=(dim-k-1)*Math.Pow(t,(dim-k-2));
					}
					hh[dim-1]=0;
				}
				for(int k=0;k<dim;k++)
				{
					double val=0;
					for(int l=0;l<dim;l++)
					{
						val+=hh[l]*M[l,k];
					}
					tt[k]=val;
				}
                for(int jj=0;jj<__DIM;jj++)
                {
					for(int j=0;j<nDV;j++)
					{
						intP[i].C[0,jj,j]=0;
					}
                }
				for(int k=0;k<nNode;k++)
				{
					//[C]
					double C=1.0;
					C*=tt[dd[k]];
					for(int j=0;j<__DIM;j++)
					{
						intP[i].C[0,j,k*__DIM+j]=C;
					}
				}
				//Create [B]  (for metric)
                intP[i].CtoB();

                //Create [D] (for second derivative)
                t = intP[i].localCoord[0];
                for (int k = 0; k < dim - 1; k++)
                {
                    hh[k] = (dim - k - 1) *(dim - k - 2) * Math.Pow(t, (dim - k - 3));
                }
                hh[dim - 1] = 0;
                hh[dim - 2] = 0;

                for (int k = 0; k < dim; k++)
                {
                    double val = 0;
                    for (int l = 0; l < dim; l++)
                    {
                        val += hh[l] * M[l, k];
                    }
                    tt[k] = val;
                }
                        
                for (int jj = 0; jj < __DIM; jj++)
                {
                    for (int j = 0; j < nDV; j++)
                    {
                        intP[i].D[0, 0,jj, j] = 0;
                    }
                }
                for (int k = 0; k < nNode; k++)
                {
                    //[D]
                    double D = 1.0;
                    D *= tt[dd[k]];
                    for (int j = 0; j < __DIM; j++)
                    {
                        intP[i].D[0, 0, j, k * __DIM + j] = D;
                    }
                }
            }
        }

    }
}
