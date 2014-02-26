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

            public List<Minilla3D.Elements.nurbsElement> elemList = new List<Elements.nurbsElement>();
            public List<Minilla3D.Elements.nurbsCurve> edgeList = new List<Elements.nurbsCurve>();
            public List<Minilla3D.Elements.nurbsCurve> topEdge = new List<Minilla3D.Elements.nurbsCurve>();
            public List<Minilla3D.Elements.nurbsCurve> bottomEdge = new List<Minilla3D.Elements.nurbsCurve>();
            public List<Minilla3D.Elements.nurbsCurve> rightEdge = new List<Minilla3D.Elements.nurbsCurve>();
            public List<Minilla3D.Elements.nurbsCurve> leftEdge = new List<Minilla3D.Elements.nurbsCurve>();
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
            public void AddEdge(Minilla3D.Elements.nurbsCurve e)
            {
                edgeList.Add(e);
            }
            public void AddRange(IEnumerable<Minilla3D.Elements.nurbsElement> collection)
            {
                elemList.AddRange(collection);
            }
            public void Clear()
            {
                elemList.Clear();
                edgeList.Clear();
                topEdge.Clear();
                bottomEdge.Clear();
                rightEdge.Clear();
                leftEdge.Clear();
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
                foreach (var e in edgeList)
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
                Parallel.ForEach(edgeList, (e) =>
                    e.setupNodesFromList(x)
                    );
            }
            public void giveEdgeTension(double t)
            {
                foreach (var e in edgeList)
                {
                    e.giveTension(t);
                }
            }
            public void computeHessian()
            {
                Parallel.ForEach(elemList, (e) =>
                    e.computeHessian()
                );
                Parallel.ForEach(edgeList, (e) =>
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
            public void GetJacobianOfCurvature(SparseDoubleArray jacobH,bool T)
            {
                int num = 0;
                foreach (var e in elemList)
                {
                    num = e.mergeJacobianOfCurvature(jacobH, num,T);
                }
            }
            public int totalNumbrOfTopEdgeIntPoint()
            {
                int num = 0;
                foreach (var e in topEdge)
                {
                    num += e.numberOfConstraintConditions();
                }
                return num;
            }
            public int totalNumbrOfBottomEdgeIntPoint()
            {
                int num = 0;
                foreach (var e in bottomEdge)
                {
                    num += e.numberOfConstraintConditions();
                }
                return num;
            }
            public int totalNumbrOfLeftEdgeIntPoint()
            {
                int num = 0;
                foreach (var e in leftEdge)
                {
                    num += e.numberOfConstraintConditions();
                }
                return num;
            }
            public int totalNumbrOfRightEdgeIntPoint()
            {
                int num = 0;
                foreach (var e in rightEdge)
                {
                    num += e.numberOfConstraintConditions();
                }
                return num;
            }
            public void GetJacobianOfTopEdge(SparseDoubleArray jacob)
            {
                int num = 0;
                foreach (var e in topEdge)
                {
                    num = e.mergeJacobianOfPosition(jacob, num);
                }
            }
            public void GetJacobianOfBottomEdge(SparseDoubleArray jacob)
            {
                int num = 0;
                foreach (var e in bottomEdge)
                {
                    num = e.mergeJacobianOfPosition(jacob, num);
                }
            }
            public void GetJacobianOfRightEdge(SparseDoubleArray jacob)
            {
                int num = 0;
                foreach (var e in rightEdge)
                {
                    num = e.mergeJacobianOfPosition(jacob, num);
                }
            }
            public void GetJacobianOfLeftEdge(SparseDoubleArray jacob)
            {
                int num = 0;
                foreach (var e in leftEdge)
                {
                    num = e.mergeJacobianOfPosition(jacob, num);
                }
            }
            public void getHessian(ShoNS.Array.SparseDoubleArray hess)
            {
                foreach (var e in elemList)
                {
                    e.mergeHessian(hess);
                }
                foreach (var e in edgeList)
                {
                    e.mergeHessian(hess);
                }
            }

            public void computeGlobalCoord()
            {
                Parallel.ForEach(elemList, (e) =>
                    e.computeGlobalCoord()
                );
                Parallel.ForEach(edgeList, (e) =>
                    e.computeGlobalCoord()
                );
            }
            public void computeAngle()
            {
                foreach (var e in elemList)
                {
                    e.computeAngle();
                }
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
            public int totalNumberOfIconst(bool T)
            {
                int num = 0;
                foreach (var e in elemList)
                {
                    num += e.numberOfConstraintConditions2(T);
                }
                return num;
            }

        }

	}
}
