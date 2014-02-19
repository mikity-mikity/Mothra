using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using ShoNS.Array;
using System.IO;
using System.Reflection;
//using System.Reactive.Linq;
namespace Minilla3D
{
	namespace Objects{
        
		public interface iObject
		{
		}
        public class masonry:iObject{

            public List<Minilla3D.Elements.nurbsElement> elemList=new List<Elements.nurbsElement>();
            public void memoryMetric()
            {
                foreach (Minilla3D.Elements.element e in elemList)
                {
                    e.memoryMetric();       //metric->refMetric,invMetric->invRefMetric,dv->refDv
                }
            }
            public void Add(Minilla3D.Elements.nurbsElement e)
            {
                elemList.Add(e);
            }
            public void AddRange(IEnumerable<Minilla3D.Elements.nurbsElement> collection)
            {
                elemList.AddRange(collection);
            }
            public void Clear()
            {
                elemList.Clear();
            }
            public void computeAiryFunction()
            {
                /*Parallel.ForEach(elemList, (e) =>
                    e.computeAiryFunction()
                );*/
                foreach(var e in elemList){
                    e.computeAiryFunction();
                }
            }
            public void precompute()
            {
                /*Parallel.ForEach(elemList, (e) =>
                    e.computeAiryFunction()
                );*/
                foreach (var e in elemList)
                {
                    e.precompute();
                }
            }
            public void computeEigenVectors()
            {
                //Parallel.ForEach(elemList, (e) =>
                 //   e.computeEigenVectors()
                //    );

                //foreach (var e in elemList)
                for(int i=0;i<elemList.Count;i++)
                {
                    var e = elemList[i];
                    e.computeEigenVectors();
                }
            }
            public void setupNodesFromList(double[,] x)
            {
                Parallel.ForEach(elemList, (e) =>
                    e.setupNodesFromList(x)
                    );
            }
            public void computeHessian()
            {
                Parallel.ForEach(elemList, (e) =>
                    e.computeHessian()
                );
            }
            public void getResidual(ShoNS.Array.DoubleArray resid)
            {
                int num = 0;
                foreach (var e in elemList)
                {
                    num = e.mergeResidual(resid, num);
                }
            }
            public void getJacobian(ShoNS.Array.SparseDoubleArray jacob)
            {
                int num = 0;
                foreach (var e in elemList)
                {
                    num=e.mergeJacobian(jacob, num);
                }
            }
            public void GetJacobianOfCurvature(SparseDoubleArray jacobH)
            {
                int num = 0;
                foreach (var e in elemList)
                {
                    num = e.mergeJacobianOfCurvature(jacobH, num);
                }
            }

            public void getHessian(ShoNS.Array.SparseDoubleArray hess)
            {
                foreach (var e in elemList)
                {
                    e.mergeHessian(hess);
                }
            }

            public void computeGlobalCoord()
            {
                Parallel.ForEach(elemList,(e) =>
                    e.computeGlobalCoord()                
                );
            }

            public int totalNumberOfIntPnts()
            {
                int num = 0;
                foreach (var e in elemList)
                {
                    num += e.nIntPoint;
                }
                return num;
            }
            public int totalNumberOfBconst()
            {
                int num = 0;
                foreach (var e in elemList)
                {
                    num += e.numberOfConstraintConditions();
                }
                return num;
            }
            public int totalNumberOfIconst()
            {
                int num = 0;
                foreach (var e in elemList)
                {
                    num += e.numberOfConstraintConditions2();
                }
                return num;
            }

        }
        /*public class generalSpring : iObject
        {
            private List<Minilla3D.Elements.managedElement> elemList;
            private double[] _grad;
            unsafe public void getGrad(double[,] acc)
            {
                fixed (double* _ptr1 = &acc[0, 0], _ptr2 = &_grad[0])
                {
                    double* ptr1 = _ptr1;
                    double* ptr2 = _ptr2;
                    for (int i = 0; i < acc.GetLength(0); i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            *ptr1 = *ptr2;
                            ptr1++;
                            ptr2++;
                        }
                    }
                }
            }
            unsafe public void getGrad(DoubleArray acc)
            {
                fixed (double* _ptr2 = &_grad[0])
                {
                    double* ptr2 = _ptr2;
                    for (int i = 0; i < acc.size1; i++)
                    {
                        acc[i] = *ptr2;
                        ptr2++;
                    }
                }
            }
                
            unsafe public void getGrad(double[] acc)
            {
                fixed (double* _ptr1 = &acc[0], _ptr2 = &_grad[0])
                {
                    double* ptr1 = _ptr1;
                    double* ptr2 = _ptr2;
                    for (int i = 0; i < _grad.Length; i++)
                    {
                        *ptr1 = *ptr2;
                        ptr1++;
                        ptr2++;
                    }
                }
            }
            unsafe void mergeGradient()
            {
                fixed (double* _ptr = &_grad[0])
                {
                    for (int i = 0; i < _grad.Length; i++)
                    {
                        _grad[i] = 0;
                    }
                    foreach (Minilla3D.Elements.managedElement e in elemList)
                    {
                        e.mergeGradient(_ptr);
                    }
                }
            }

            public void memoryMetric()
            {
                foreach (Minilla3D.Elements.managedElement e in elemList)
                {
                    e.memoryMetric();
                }
            }
            public void memoryVolume()
            {
                foreach (Minilla3D.Elements.managedElement e in elemList)
                {
                    e.memoryVolume();
                }
            }
            public generalSpring()
            {
                elemList = new List<Minilla3D.Elements.managedElement>();
            }
            public void Add(Minilla3D.Elements.managedElement e)
            {
                elemList.Add(e);
            }
            public void AddRange(IEnumerable<Minilla3D.Elements.managedElement> collection)
            {
                elemList.AddRange(collection);
            }
            public void Clear()
            {
                elemList.Clear();
            }
            public void initialize(int nParticles)
            {
                _grad = new double[nParticles * 3];
            }
            public void setMaterial(Minilla3D.Materials.material _m)
            {
                foreach (Minilla3D.Elements.managedElement e in elemList)
                {
                    e.setMaterial(_m);
                }
            }
            public void computeEigenVectors()
            {
                Parallel.For(0, elemList.Count, (i) =>
                {
                    elemList[i].computeEigenVectors();
                });
            }
            unsafe public void computeHessEd(double[,] x)
            {
                for (int i = 0; i < elemList.Count; i++)
                {
                    elemList[i].computeHessEd();
                }

            }
            unsafe public void getHessian(ShoNS.Array.SparseDoubleArray hess, int dim)
            {
                //hess.Clear();
                //hess.RemoveZeros();
                //hess.FillValue(0);
                //hess.RemoveZeros();
                //z座標に関するものを落とす。
                foreach (var e in elemList)
                {
                    e.mergeHessian(hess, dim);
                }
            }

            unsafe public void computeAll(double[,] x)
            {
                fixed (double* ptr = &x[0, 0])
                {
                    double* ptr1 = ptr;
                    Parallel.For(0,elemList.Count,(i)=>
                        {
                            elemList[i].setupNodesFromList(ptr1);
                            elemList[i].computeGlobalCoord();
                            elemList[i].computeBaseVectors();
                            elemList[i].computeMetric();
                            elemList[i].computeVolume();
                            elemList[i].computeStress();
                            elemList[i].computeGradient();
                        }
                    );
                }
                if (_grad == null)
                {
                    _grad = new double[x.Length];
                }
                else if (_grad.Length != x.Length)
                {
                    _grad = new double[x.Length];
                }
                mergeGradient();
            }
            unsafe public void computeVolume(double[,] x)
            {
                fixed (double* ptr = &x[0, 0])
                {
                    double* ptr1 = ptr;
                    Parallel.For(0, elemList.Count, (i) =>
                    {
                        elemList[i].setupNodesFromList(ptr1);
                        elemList[i].computeMetric();
                        elemList[i].computeVolume();
                    }
                    );
                }
            }
            public double getTotalVolume()
            {
                double value=0;
                for (int i = 0; i < elemList.Count; i++)
                {
                    value += elemList[i].getVolume();
                }

                return value;
            }

            public void getStar(List<double[]>[] star)
            {
                foreach (var e in elemList)
                {
                    e.createStar(star);
                }
            }
        }

		public class constraineVolume:iObject
		{
			private List<Minilla3D.Elements.managedElement> elemList;
			private double[] _grad;
            private double _refVolume,_volume;
            unsafe public void getGrad(double[,] acc)
            {
                fixed (double* _ptr1 = &acc[0, 0],_ptr2=&_grad[0])
                {
                    double* ptr1 = _ptr1;
                    double* ptr2 = _ptr2;
                    for (int i = 0; i < _grad.Length; i++)
                    {
                        *ptr1 = *ptr2;
                        ptr1++;
                        ptr2++;
                    }
                }
            }
            unsafe public void getGrad(double[] acc)
            {
                fixed (double* _ptr1 = &acc[0], _ptr2 = &_grad[0])
                {
                    double* ptr1 = _ptr1;
                    double* ptr2 = _ptr2;
                    for (int i = 0; i < _grad.Length; i++)
                    {
                        *ptr1 = *ptr2;
                        ptr1++;
                        ptr2++;
                    }
                }
            }
            unsafe public void getGrad(double[,] J,int index)
            {
                fixed (double* _ptr1 = &J[index,0], _ptr2 = &_grad[0])
                {
                    double* ptr1 = _ptr1;
                    double* ptr2 = _ptr2;
                    for (int i = 0; i < _grad.Length; i++)
                    {
                        *ptr1 = *ptr2;
                        ptr1++;
                        ptr2++;
                    }
                }
            }
            public void getResidual(out double r)
            {
                r = _volume - _refVolume;
            }
			unsafe void mergeGradient()
			{
                fixed(double*_ptr=&_grad[0])
                {
				    for(int i=0;i<_grad.Length;i++)
				    {
					    _grad[i]=0;
				    }
				    foreach (Minilla3D.Elements.managedElement e in elemList)
				    {
					    e.mergeGradient(_ptr);
				    }
                }
            }
            public void memoryMetric()
            {
                foreach (Minilla3D.Elements.managedElement e in elemList)
                {
                    e.memoryMetric();
                }
            }
            public void setRefVolume(double v)
            {
                _refVolume=v;
            }
			public constraineVolume()
			{
				elemList=new List<Minilla3D.Elements.managedElement>();
			}
    		public void Add(Minilla3D.Elements.managedElement e)
			{
				elemList.Add(e);
                Materials.defaultMaterial dM=new Materials.defaultMaterial();
                foreach(Minilla3D.Elements.managedElement _e in elemList)
                {
                    _e.setMaterial(dM.getMaterial());
                }
			}
    		public void AddRange(IEnumerable<Minilla3D.Elements.managedElement> collection)
	    	{
		    	elemList.AddRange(collection);
                Materials.defaultMaterial dM=new Materials.defaultMaterial();
                foreach(Minilla3D.Elements.managedElement _e in elemList)
                {
                    _e.setMaterial(dM.getMaterial());
                }
			}
            public void Clear()
            {
                elemList.Clear();
            }
            public void initialize(int nParticles)
            {
                _grad = new double[nParticles * 3];
            }
			unsafe public void computeAll(double[,] x)
			{
                fixed(double* ptr=&x[0,0])
                {
                    double* ptr1=ptr;
                    for (int i = 0; i < (int)(elemList.Count); i++)
                    {
                        elemList[i].setupNodesFromList(ptr1);
                        elemList[i].computeMetric();
                        elemList[i].computeVolume();
                        elemList[i].computeStress();
                        elemList[i].computeGradient();
                    }
                }
    			if(_grad==null)
	    		{
		    		_grad=new double[x.Length];
			    }else if(_grad.Length!=x.Length)
				{
					_grad=new double[x.Length];
				}
				mergeGradient();
                double v = 0;
                for (int i = 0; i < elemList.Count; i++)
                {
                    v += elemList[i].getVolume();
                }
                this._volume = v;
            }
		}*/
	}
}
