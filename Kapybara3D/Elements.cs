using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Minilla3D.Elements
{
    
    public class element
	{
        enum type{
            Cauchy,SPK
        };
        type typeOfStress;
        protected const int __DIM=3;
        public int nIntPoint,nBIntPoint,nNode,elemDim,nDV;                
		public Minilla3D.Elements.integratingPoint[] intP;     //Integrating points
        public Minilla3D.Elements.integratingPoint[] bIntP;     //Integrating Points on border
        double[] node;							        //Nodal coordinate (global)
		protected int[] index;                                    //indeces of the nodes
        double[] gradient;                              //internal force(equivalent nodal force of stress field)
        double[,] hess;                                 //Hessian  (Geometric Stiffness only)
        double[] force;                                 //external force(equivalent nodal force of gravity)
		
		public element(int _nNode,int _elemDim,int _nIntPoint)
		{
            var _index=new int[nDV];
			for(int i=0;i<nDV;i++)
			{
				_index[i]=i;
			}
            initialize(_index, _nNode, _elemDim, _nIntPoint);
		}
        public element(int[] _index,int _nNode,int _elemDim,int _nIntPoint)
		{
            initialize(_index,_nNode,_elemDim,_nIntPoint);
		}
        private void initialize(int[] _index,int _nNode,int _elemDim,int _nIntPoint)
        {
            typeOfStress = type.Cauchy;   //default value
            nIntPoint=_nIntPoint;
            nBIntPoint = 0;
            nNode=_nNode;
            elemDim=_elemDim;
            nDV=nNode*__DIM;
            intP=new integratingPoint[nIntPoint];
            this.bIntP = new integratingPoint[0];
            for (int i = 0; i < nIntPoint; i++)
            {
                intP[i]=new integratingPoint(nNode,elemDim);
            }
            node=new double[nDV];
            index=new int[nNode];
		    gradient=new double[nDV];
		    hess=new double [nDV,nDV];
		    force=new double[nDV];
            setupIndex(_index);
        }
        public void setupIndex(int[] _index)
		{
            Array.Copy(_index,index,nNode);
		}
        public double[] getIntPoint(int i)
        {
            return intP[i].globalCoord;
        }
        public double[] getBIntPoint(int i)
        {
            return bIntP[i].globalCoord;
        }
        public double[] getNode(int i)
		{
			return new double[3]{node[i*__DIM+0],node[i*__DIM+1],node[i*__DIM+2]};
	    }
		public int getIndex(int i)
		{
			return index[i];
		}
		public void setupNodes(double[,] x)
		{
            Array.Copy(x,node,nDV);
		}
		public void setupNodesFromList(double[,] x)
		{
			for(int i=0;i<nNode;i++)
			{
				for(int j=0;j<__DIM;j++)
				{
					node[i*__DIM+j]=x[index[i],j];
				}
			}
		}

		public double Volume{
	        get;
            protected set;
        }
		public double refVolume{
	        get;
            protected set;
        }
        public void computeAiryFunction()
        {
            for (int i = 0; i < nIntPoint; i++)
            {
                intP[i].computeAiryFunction(node);
            }
            for (int i = 0; i < nBIntPoint; i++)
            {
                bIntP[i].computeAiryFunction(node);
            }
        }
        public void precompute()
        {
            for (int i = 0; i < nIntPoint; i++)
            {
                intP[i].precompute(node);
            }
            for (int i = 0; i < nBIntPoint; i++)
            {
                bIntP[i].precompute(node);
            }
        }
        public void computeEigenVectors()
        {
            for(int i=0;i<nIntPoint;i++)
            {
                intP[i].computeEigenVectors();
            }
        }
        public void getEigenVectors(double[][] vec, double[] val, int num)
        {
            for (int i = 0; i < 2; i++)
            {
                vec[i][0] = intP[num].eigenVectors[i][0];
                vec[i][1] = intP[num].eigenVectors[i][1];
                vec[i][2] = intP[num].eigenVectors[i][2];
                val[i] = intP[num].eigenValues[i];
            }
        }
        public void computeGlobalCoord()
		{
            for (int i = 0; i < nIntPoint; i++)
            {
                intP[i].computeGlobalCoord(node);
            }
            for (int i = 0; i < nBIntPoint; i++)
            {
                bIntP[i].computeGlobalCoord(node);
            }
        }
		public void computeMetric(){
			for(int i=0;i<nIntPoint;i++)
			{
				intP[i].computeMetric(node);
			}
		}
		public void computeBaseVectors()
		{
			for(int i=0;i<nIntPoint;i++)
			{
				intP[i].computeBaseVectors(node);
			}
		}
		public double computeVolume(){
			double v=0;
			for(int i=0;i<nIntPoint;i++)
			{
				v+=intP[i].weight*intP[i].dv;
			}
			this.Volume=v;
			return v;
		}
		public void memoryVolume(){
			this.refVolume=this.Volume;
		}
		public void setRefVolume(double v){
			this.refVolume=v;
		}
		public void memoryMetric(){
			for(int i=0;i<nIntPoint;i++)
			{
				intP[i].memoryMetric();
			}
		}
		public void computeHessian()
		{
			for(int i=0;i<nDV;i++)
			{
				for(int j=0;j<nDV;j++)
				{
					hess[i,j]=0;
				}
			}
			for(int i=0;i<nIntPoint;i++)
			{
				double val=intP[i].weight*intP[i].refDv*0.5;
				for(int j=0;j<elemDim;j++)
				{
				    for(int jj=0;jj<elemDim;jj++)
				    {
				        for(int k=0;k<nDV;k++)
				        {
					        for(int l=0;l<nDV;l++)
					        {
						        double D=intP[i].SPK[j,jj];
						        double S=intP[i].B[j,jj,k,l];
						        hess[k,l]+=val*D*S;
					        }
				        }
                    }
				}
			}
		}

		public void computeGradient()
		{
			
			for(int i=0;i<nDV;i++)
			{
				gradient[i]=0;
			}
            if(typeOfStress==type.Cauchy){
			    for(int i=0;i<nIntPoint;i++)
			    {
				    double val=intP[i].weight*intP[i].dv*0.5;
				    for(int j=0;j<elemDim;j++)
				    {
                        for(int jj=0;jj<elemDim;jj++)
                        {
					        for(int k=0;k<nDV;k++)
					        {
						        for(int l=0;l<nDV;l++)
						        {
							        double D=intP[i].Cauchy[j,jj];
							        double S=intP[i].B[j,jj,k,l];
							        double E=node[l];
							        gradient[k]+=val*D*S*E;
						        }
					        }
                        }
				    }
			    }
            }else{
			    for(int i=0;i<nIntPoint;i++)
			    {
				    double val=intP[i].weight*intP[i].refDv*0.5;
				    for(int j=0;j<elemDim;j++)
				    {
                        for(int jj=0;jj<elemDim;jj++)
                        {
					        for(int k=0;k<nDV;k++)
					        {
						        for(int l=0;l<nDV;l++)
						        {
							        double D=intP[i].SPK[j,jj];
							        double S=intP[i].B[j,jj,k,l];
							        double E=node[l];
							        gradient[k]+=val*D*S*E;
						        }
					        }
                        }
				    }
			    }
            }
		}
		/*virtual void copyGradient(double* ptr,int DOF)
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
		}*/
        /*
		virtual void mergeGradient(double *ptr)
		{
			for(int i=0;i<_nNode;i++)
			{
				for(int j=0;j<__DIM;j++)
				{
					ptr[this->index[i]*__DIM+j]+=this->gradient[i*__DIM+j];
				}
			}
		}*/
		public void mergeHessian(ShoNS.Array.SparseDoubleArray _hess)
		{
			for(int i=0;i<nNode;i++)
			{
				for(int j=0;j<__DIM;j++)
				{
					for(int k=0;k<nNode;k++)
					{
						for(int l=0;l<__DIM;l++)
						{
                            _hess[this.index[i] * __DIM + j, this.index[k] * __DIM + l] += this.hess[i * __DIM + j, k * __DIM + l];
						}
					}
				}
			}
		}
	}
}
