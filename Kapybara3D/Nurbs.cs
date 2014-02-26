using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Minilla3D.Elements
{
    public class nurbsElement:element
    {
        double[][] _cu;
		double[][] _pu;
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
        public int numberOfConstraintConditions()
        {
            if (exceptional()) return 0;
            return nBIntPoint * 2;
        }
        public int numberOfConstraintConditions2(bool T)
        {
            if (T&&exceptional()) return 0;
            return nIntPoint * 3;
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
        public int mergeJacobianOfCurvature(ShoNS.Array.SparseDoubleArray jacobH, int num,bool T)
        {
            if (T&&exceptional()) return num;
            for (int i = 0; i < nIntPoint; i++)
            {
                var grad = intP[i].getGradientOfH(0, 0);
                for (int f = 0; f < nNode; f++)
                {
                    jacobH[num + i * 3 + 0, index[f]] = grad[f];
                }
                grad = intP[i].getGradientOfH(1, 1);
                for (int f = 0; f < nNode; f++)
                {
                    jacobH[num + i * 3 + 1, index[f]] = grad[f];
                }
                grad = intP[i].getGradientOfH(0, 1);
                for (int f = 0; f < nNode; f++)
                {
                    jacobH[num + i * 3 + 2, index[f]] = grad[f];
                }
            }
            return num + nIntPoint*3;
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
        public void stitch(params nurbsCurve[] elems)
        {
            int count = 0;
            foreach (var e in elems)
            {
                for (int i = 0; i < e.nIntPoint; i++, count++)
                {
                    bIntP[count].refIntP = e.intP[i];
                    e.intP[i].refIntP = bIntP[count];
                }
            }
        }
        public void setPlane(double a, double b, double c, double d)
        {
            foreach (var p in intP)
            {
                p.refIntP.setPlane(a, b, c, d);
            }
        }
        public void computeAngle()
        {
            foreach (var p in bIntP)
            {
                p.computeAngle(this.node);
                p.transferAngleToTension();
            }
        }
        public nurbsElement(int _uDim, int _vDim, int[] _index, int uNum, int vNum, double[] uKnot, double[] vKnot, border _border = border.None)
            : base(_index, _uDim * _vDim, 2, (_uDim + 1) * (_vDim + 1))
        {

            this.typeOfBorder = _border;
            int bVdim = _vDim +1;
            int bUdim = _uDim +1;
            switch (_border)
            {
                case border.Left | border.Top|border.Right:
                    this.bIntP = new integratingPoint[bVdim * 2 + bUdim];
                    nBIntPoint = bVdim * 2 + bUdim;
                    break;
                case border.Left | border.Bottom | border.Right:
                    this.bIntP = new integratingPoint[bVdim * 2 + bUdim];
                    nBIntPoint = bVdim * 2 + bUdim;
                    break;
                case border.Left | border.Top | border.Bottom:
                    this.bIntP = new integratingPoint[bUdim * 2 + bVdim];
                    nBIntPoint = bUdim * 2 + bVdim;
                    break;
                case border.Right | border.Top | border.Bottom:
                    this.bIntP = new integratingPoint[bUdim * 2 + bVdim];
                    nBIntPoint = bUdim * 2 + bVdim;
                    break;
                case border.Left | border.Right:
                    this.bIntP = new integratingPoint[bVdim + bVdim];
                    nBIntPoint = bVdim + bVdim;
                    break;
                case border.Top | border.Bottom:
                    this.bIntP = new integratingPoint[bUdim + bUdim];
                    nBIntPoint = bUdim + bUdim;
                    break;
                case border.Left | border.Top:
                    this.bIntP = new integratingPoint[bUdim + bVdim];
                    nBIntPoint = bUdim + bVdim;
                    break;
                case border.Left | border.Bottom:
                    this.bIntP = new integratingPoint[bUdim + bVdim];
                    nBIntPoint = bUdim + bVdim;
                    break;
                case border.Right | border.Top:
                    this.bIntP = new integratingPoint[bUdim + bVdim];
                    nBIntPoint = bUdim + bVdim;
                    break;
                case border.Right | border.Bottom:
                    this.bIntP = new integratingPoint[bUdim + bVdim];
                    nBIntPoint = bUdim + bVdim;
                    break;
                case border.Left:
                    this.bIntP = new integratingPoint[bVdim];
                    nBIntPoint = bVdim;
                    break;
                case border.Right:
                    this.bIntP = new integratingPoint[bVdim];
                    nBIntPoint = bVdim;
                    break;
                case border.Top:
                    this.bIntP = new integratingPoint[bUdim];
                    nBIntPoint = bUdim;
                    break;
                case border.Bottom:
                    this.bIntP = new integratingPoint[bUdim];
                    nBIntPoint = bUdim;
                    break;
                default:
                    this.bIntP = new integratingPoint[0];
                    nBIntPoint = 0;
                    break;
            }
            for (int i = 0; i < nBIntPoint; i++)
            {
                bIntP[i] = new integratingPoint(nNode, elemDim);
            }
            uDim = _uDim;
            vDim=_vDim;
            _cu = new double[elemDim][];
            _cu[0] = new double[uDim + 1];
            _cu[1] = new double[vDim + 1];
            _pu = new double[elemDim][];
            _pu[0] = new double[uDim + 1];
            _pu[1] = new double[vDim + 1];
            int[,] ss = new int[nIntPoint, elemDim];		//Indeces for integrating points
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
            for (int i = 0; i < uDim+1; i++)
            {
                _cu[0][i] = __cu(i, uDim);
                _pu[0][i] = __pu(i, uDim);
            }
            for (int i = 0; i < vDim+1; i++)
            {
                _cu[1][i] = __cu(i, vDim);
                _pu[1][i] = __pu(i, vDim);
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
			//weight coefficients and loacl coordinates for integrating points
			for(int i=0;i<nIntPoint;i++)
			{
				intP[i].weight=1.0;
				for(int j=0;j<elemDim;j++)
				{
					intP[i].localCoord[j]=_cu[j][ss[i,j]];
					intP[i].weight*=_pu[j][ss[i,j]];
				}
			}
            //local coordinates for integrating points on border
            int count = 0;
            switch (_border)
            {
                case border.Left | border.Top | border.Right:
                    count = 0;
                    for (int i = 0; i < bVdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = 0;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, -1 };
                    }
                    for (int i = 0; i < bUdim; i++, count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 0;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { 1, 0 };
                    }
                    for (int i = 0; i < bVdim; i++, count++)
                    {
                        bIntP[count].localCoord[0] = 1;
                        bIntP[count].localCoord[1] = __cu(i, bVdim - 1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, 1 };
                    }
                    break;
                case border.Left | border.Bottom | border.Right:
                    count = 0;
                    for (int i = 0; i < bVdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = 0;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, -1 };
                    }
                    for (int i = 0; i < bUdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 1;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { -1, 0 };
                    }
                    for (int i = 0; i < bVdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = 1;
                        bIntP[count].localCoord[1] = __cu(i,bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, 1 };
                    }
                    break;
                case border.Left | border.Top | border.Bottom:
                    count = 0;
                    for (int i = 0; i < bVdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = 0;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, -1 };
                    }
                    for (int i = 0; i < bUdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 0;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { 1, 0 };
                    }
                    for (int i = 0; i <bUdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 1;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { -1, 0 };
                    }
                    break;
                case border.Right | border.Top | border.Bottom:
                    count = 0;
                    for (int i = 0; i < bVdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = 1;
                        bIntP[count].localCoord[1] = __cu(i, bVdim - 1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, 1 };
                    }
                    for (int i = 0; i < bUdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim - 1);
                        bIntP[count].localCoord[1] = 0;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { 1, 0 };
                    }
                    for (int i = 0; i < bUdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 1;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { -1, 0 };
                    }
                    break;
                case border.Left | border.Right:
                    count = 0;
                    for (int i = 0; i < bVdim; i++, count++)
                    {
                        bIntP[count].localCoord[0] = 0;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, -1 };
                    }
                    for (int i = 0; i < bVdim; i++, count++)
                    {
                        bIntP[count].localCoord[0] = 1;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, 1 };
                    }
                    break;
                case border.Top | border.Bottom:
                    count = 0;
                    for (int i = 0; i < bUdim; i++, count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 0;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { 1, 0 };
                    }
                    for (int i = 0; i < bUdim; i++, count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 1;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { -1, 0 };
                    }
                    break;
                case border.Left | border.Top:
                    count = 0;
                    for (int i = 0; i < bVdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = 0;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, -1 };
                    }
                    for (int i = 0; i < bUdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 0;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { 1, 0 };
                    }
                    break;
                case border.Left | border.Bottom:
                    count=0;
                    for (int i = 0; i < bVdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = 0;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, -1 };
                    }
                    for (int i = 0; i < bUdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 1;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { -1, 0 };
                    }
                    break;
                case border.Right | border.Top:
                    count=0;
                    for (int i = 0; i < bVdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = 1;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, 1 };
                    }
                    for (int i = 0; i < bUdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 0;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { 1, 0 };
                    }
                    break;
                case border.Right | border.Bottom:
                    count=0;
                    for (int i = 0; i < bVdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = 1;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, 1 };
                    }
                    for (int i = 0; i < bUdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 1;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { -1, 0 };
                    }
                    break;
                case border.Left:
                    count=0;
                    for (int i = 0; i < bVdim; i++, count++)
                    {
                        bIntP[count].localCoord[0] = 0;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, -1 };
                    }
                    break;
                case border.Right:
                    count=0;
                    for (int i = 0; i < bVdim; i++, count++)
                    {
                        bIntP[count].localCoord[0] = 1;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, 1 };
                    }
                    break;
                case border.Top:
                    count=0;
                    for (int i = 0; i < bUdim; i++, count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 0;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { 1, 0 };
                    }
                    break;
                case border.Bottom:
                    count=0;
                    for (int i = 0; i < bUdim; i++, count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 1;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { -1, 0 };
                    }
                    break;
                default:
                    break;
            }
            List<Minilla3D.Elements.integratingPoint> allIntP = new List<integratingPoint>();
            allIntP.AddRange(intP);
            allIntP.AddRange(bIntP);
            int nAllIntP = nIntPoint + nBIntPoint;
		    double[][,] M=new double[2][,];
		    M[0]=fM(uNum,_uDim,_uDim-1,uKnot);
		    M[1]=fM(vNum,_vDim,_vDim-1,vKnot);
			//Shape functions  [N] (for global coordinate)
			for(int i=0;i<nAllIntP;i++)
			{
				for(int j=0;j<elemDim;j++)
				{
					double t=allIntP[i].localCoord[j];
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
					    allIntP[i].N[j,k]=0;
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
						allIntP[i].N[j,k*__DIM+j]=N;
					}
				}

                //Create [C]  (for base vectors)
				for (int m=0;m<elemDim;m++)
				{
					for(int j=0;j<elemDim;j++)
					{
						double t=allIntP[i].localCoord[j];
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
						    allIntP[i].C[m,jj,j]=0;
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
							allIntP[i].C[m,j,k*__DIM+j]=C;
						}
					}
                }
				//Create [B]  (for metric)
                allIntP[i].CtoB();

                //Create [D] (for second derivative)
                for (int m = 0; m < elemDim; m++)
                {
                    for (int n = 0; n < elemDim; n++)
                    {
                        for (int j = 0; j < elemDim; j++)
                        {
                            double t = allIntP[i].localCoord[j];
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
                                allIntP[i].D[m, n,jj, j] = 0;
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
                                allIntP[i].D[m,n, j, k * __DIM + j] = D;
                            }
                        }
                    }
                }
            }
        }

    }
}
