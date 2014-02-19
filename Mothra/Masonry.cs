using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ShoNS.Array;
using Rhino.Geometry;
using System.IO;
using System.Reflection;
using Gurobi;
using Microsoft.SolverFoundation.Common;
using Microsoft.SolverFoundation.Solvers;
namespace mikity.ghComponents
{

    /// <summary>
    /// Construct a point array using isoparametric shape functions.
    /// </summary>
    public class Mothra : Grasshopper.Kernel.GH_Component
    {
        
        Rhino.Geometry.NurbsSurface inputNurbs, airyNurbs, xyNurbs,outputNurbs;
        Minilla3D.Objects.masonry myMasonry = new Minilla3D.Objects.masonry();
        double[,] x;
        bool initialized = false;
        int nU, nV;
        List<Guid> fixedPointGuids = null;
        Rhino.Geometry.Mesh alternativeMesh = null;
        List<int> boundaryIndex=null;
        bool __update = false;
        System.Windows.Forms.Timer timer=null;
        //List<Minilla3D.Elements.managedElement> elemList = new List<Minilla3D.Elements.managedElement>();

        /*        protected override System.Drawing.Bitmap Icon
                {
                    get
                    {
                        //現在のコードを実行しているAssemblyを取得
                        System.Reflection.Assembly myAssembly =
                            System.Reflection.Assembly.GetExecutingAssembly();

                        System.IO.Stream st = myAssembly.GetManifestResourceStream("mikity.ghComponents.icons.icon46.bmp");
                        //指定されたマニフェストリソースを読み込む
                        System.Drawing.Bitmap bmp = new System.Drawing.Bitmap(st);
                        return bmp;
                    }
                }
        */
        public Mothra()
            : base("Mothra", "Mothra", "Mothra", "Kapybara3D", "Computation")
        {
        }
        public override Guid ComponentGuid
        {
            get { return new Guid("13a84560-1b45-401b-bad8-b8487a55d51e"); }
        }

        protected override void RegisterInputParams(Grasshopper.Kernel.GH_Component.GH_InputParamManager pManager)
        {
            
            pManager.AddGenericParameter("Surface", "S", "inputSurface", Grasshopper.Kernel.GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(Grasshopper.Kernel.GH_Component.GH_OutputParamManager pManager)
        {
        }
        void timer_Tick(object sender, EventArgs e)
        {
            if (__update)
            {
                __update = false;
                update();
                this.ExpirePreview(true);
            }
        }
        //System.Windows.Forms.Timer timer;
        public override void AddedToDocument(Grasshopper.Kernel.GH_Document document)
        {
            base.AddedToDocument(document);
            Rhino.RhinoDoc.ReplaceRhinoObject += RhinoDoc_ReplaceRhinoObject;
            timer = new System.Windows.Forms.Timer();
            timer.Tick += timer_Tick;
            timer.Enabled = true;
            timer.Interval = 30;
        }
        public override void RemovedFromDocument(Grasshopper.Kernel.GH_Document document)
        {
            base.RemovedFromDocument(document);
            deleteFixedPoints();
        }
        public override void DocumentContextChanged(Grasshopper.Kernel.GH_Document document, Grasshopper.Kernel.GH_DocumentContext context)
        {
            if (context == Grasshopper.Kernel.GH_DocumentContext.Unloaded)
            {
                deleteFixedPoints();
            }
            base.DocumentContextChanged(document, context);
        }
        void RhinoDoc_ReplaceRhinoObject(object sender, Rhino.DocObjects.RhinoReplaceObjectEventArgs e)
        {
            if (initialized && fixedPointGuids != null)
            {
                for (int i = 0; i < nV; i++)
                {
                    for (int j = 0; j < nU; j++)
                    {
                        Guid gi = fixedPointGuids[i*nU+j];
                        if (e.ObjectId.CompareTo(gi) == 0)
                        {
                            var PO = e.NewRhinoObject as Rhino.DocObjects.PointObject;
                            var P1 = PO.PointGeometry.Location;
                            var P2 = airyNurbs.Points.GetControlPoint(j, i).Location;
                            if (P1.X != P2.X || P1.Y != P2.Y)
                            {
                                Rhino.RhinoDoc.ActiveDoc.Objects.Replace(gi, new Rhino.Geometry.Point3d(P2.X, P2.Y, P1.Z));
                                airyNurbs.Points.SetControlPoint(j, i, new Rhino.Geometry.ControlPoint(P2.X, P2.Y, P1.Z));
                                __update = true;
//                                this.ExpirePreview(true);
                                return;
                            }else
                            {
                                airyNurbs.Points.SetControlPoint(j, i, new Rhino.Geometry.ControlPoint(P1));
                                __update = true;
//                                this.ExpirePreview(true);
                                return;
                            }
                        }
                    }
                }
            }
        }
        bool isBoundary(int n)
        {
            for (int i = 0; i < boundaryIndex.Count(); i++)
            {
                if (boundaryIndex[i] == n)
                {
                    return true;
                }
            }
            return false;
        }
        public void update()
        {
            Nurbs2x(airyNurbs, x);
            double norm = 0;
            double minZ = Double.MaxValue;
            for (int i = 0; i < nU * nV; i++)
            {
                if (x[i, 2] < minZ) minZ = x[i, 2];
            }
            for (int i = 0; i < nU * nV; i++)
            {
                x[i, 2]-=minZ;
            }
            for (int i = 0; i < nU * nV; i++)
            {
                norm += x[i, 2] * x[i, 2];
            }
            norm = Math.Sqrt(norm);
            for (int i = 0; i < nU * nV; i++)
            {
                x[i, 2] /= norm/500;
            }            
            x2Nurbs(x, airyNurbs);
            myMasonry.setupNodesFromList(x);
            myMasonry.precompute();
            int nConstraints2 = myMasonry.totalNumberOfIconst();
            int nConstraints = myMasonry.totalNumberOfBconst();
            ShoNS.Array.SparseDoubleArray jacobian = new SparseDoubleArray(nConstraints, nU * nV);
            ShoNS.Array.DoubleArray residual = new DoubleArray(nConstraints, 1);
            
            ShoNS.Array.SparseDoubleArray jacobianH = new SparseDoubleArray(nConstraints2, nU * nV);
            myMasonry.GetJacobianOfCurvature(jacobianH);
            //myMasonry.GetJacobianOfGaussianCurvature(jacobianK);

            myMasonry.getJacobian(jacobian);
            myMasonry.getResidual(residual);
            System.Windows.Forms.MessageBox.Show(residual.Norm().ToString());
            
            try
            {
                ////////////////////////////GUROBI///////////////////
                GRBEnv env = new GRBEnv("qcp.log");
                GRBModel model = new GRBModel(env);
                double[] lb = new double[nU * nV];
                double[] ub = new double[nU * nV];
                char[] type = new char[nU * nV];
                string[] name = new string[nU * nV];
                int count = 0;
                for (int j = 0; j < nV; j++)
                {
                    for (int i = 0; i < nU; i++, count++)
                    {
                        lb[count] = 0;
                        ub[count] = GRB.INFINITY;
                        type[count] = GRB.CONTINUOUS;
                        name[count] = "N[" + i.ToString("g") + "," + j.ToString("g") + "]";
                    }
                }
                lb[0] = 250;
                ub[0] = 250;
                lb[nU - 1] = 250;
                ub[nU - 1] = 250;
                lb[nU * nV - nU] = 250;
                ub[nU * nV - nU] = 250;
                lb[nU * nV - 1] = 250;
                ub[nU * nV - 1] = 250;
                

                lb[(nU / 2) + (nV / 2) * nU] = 0;
                ub[(nU / 2) + (nV / 2) * nU] = 0;
                GRBVar[] vars = model.AddVars(lb, ub, null, type, name);
                model.Update();
                
                
                lb = new double[nConstraints2];
                ub = new double[nConstraints2];
                type = new char[nConstraints2];
                name = new string[nConstraints2];

                for (int j = 0; j < nConstraints2; j+=3)
                {
                    lb[j] = 0;
                    ub[j] = GRB.INFINITY;
                    type[j] = GRB.CONTINUOUS;
                    name[j] = "P" + j.ToString("g") + "-" + "phi_{0,0}";
                    lb[j+1] = 0;
                    ub[j+1] = GRB.INFINITY;
                    type[j+1] = GRB.CONTINUOUS;
                    name[j+1] = "P" + j.ToString("g") + "-" + "phi_{1,1}";
                    lb[j + 2] = -GRB.INFINITY;
                    ub[j + 2] = GRB.INFINITY;
                    type[j+2] = GRB.CONTINUOUS;
                    name[j + 2] = "P" + j.ToString("g") + "-" + "phi_{0,1}";
                }
                GRBVar[] phis = model.AddVars(lb, ub, null, type, name);

                model.Update();
                

                //Add Hessian entries to variables
                GRBLinExpr[] exprs = new GRBLinExpr[nConstraints2];
                double[] rhs = new double[nConstraints2];
                char[] senses = new char[nConstraints2];

                for (int i = 0; i < nConstraints2; i++)
                {
                    exprs[i] = new GRBLinExpr();
                    rhs[i] = 0;
                    senses[i] = GRB.EQUAL;
                }
                foreach (var e in jacobianH.Elements)
                {
                    exprs[e.Row].AddTerm(e.Value, vars[e.Col]);
                }
                for (int i = 0; i < nConstraints2; i++)
                {
                    exprs[i].AddTerm(-1, phis[i]);
                }
                model.AddConstrs(exprs, senses, rhs, null);
                
                //Quadratic constraints  detH>0
                for (int i = 0; i < nConstraints2; i+=3)
                {
                    GRBQuadExpr Qlhs = new GRBQuadExpr(phis[i +2] * phis[i+ 2]);
                    GRBQuadExpr Qrhs = new GRBQuadExpr(phis[i +0] * phis[i + 1]);
                    model.AddQConstr(Qlhs, GRB.LESS_EQUAL, Qrhs, null);
                }
                


                //Boundary variable
                lb = new double[nConstraints];
                ub = new double[nConstraints];
                type = new char[nConstraints];
                name = new string[nConstraints];
                for (int j = 0; j < nConstraints; j++)
                {
                    lb[j] = -100;// -GRB.INFINITY;
                    ub[j] = 100;// GRB.INFINITY;
                    type[j] = GRB.CONTINUOUS;
                    name[j] = "B" + j.ToString("g");
                }
                GRBVar[] bVars = model.AddVars(lb, ub, null, type, name);
                model.Update();
                GRBLinExpr[] bConds = new GRBLinExpr[nConstraints];
                double[] bRhs = new double[nConstraints];
                char[] bSenses = new char[nConstraints];
                for (int i = 0; i < nConstraints; i++)
                {
                    bConds[i] = new GRBLinExpr();
                    bRhs[i] = 0;
                    bSenses[i] = GRB.EQUAL;
                }
                for (int i = 0; i < nConstraints; i++)
                {
                    for (int j = 0; j < nU * nV; j++)
                    {
                        if (jacobian[i, j] != 0)
                        {
                            bConds[i].AddTerm(jacobian[i, j], vars[j]);
                        }
                    }
                    bConds[i].AddTerm(-1, bVars[i]);
                }
                
                model.AddConstrs(bConds, bSenses, bRhs, null);
                //Objective
                GRBQuadExpr obj = new GRBQuadExpr();
                for (int i = 0; i < nConstraints; i++)
                {
                    obj.AddTerm(1, bVars[i], bVars[i]);
                }
                model.SetObjective(obj);
                
                model.Optimize();
                switch (model.Get(GRB.IntAttr.Status))
                {
                    case GRB.Status.OPTIMAL:
                        System.Windows.Forms.MessageBox.Show("Congraturation");
                        break;
                    case GRB.Status.UNBOUNDED:
                        System.Windows.Forms.MessageBox.Show("UNBOUNDED");
                        break;
                    case GRB.Status.CUTOFF:
                        System.Windows.Forms.MessageBox.Show("CUTOFF");
                        break;
                    case GRB.Status.INF_OR_UNBD:
                        System.Windows.Forms.MessageBox.Show("INF_OR_UNBD");
                        break;
                    case GRB.Status.INFEASIBLE:
                        System.Windows.Forms.MessageBox.Show("INFEASIBLE");
                        break;
                    case GRB.Status.INPROGRESS:
                        System.Windows.Forms.MessageBox.Show("INPROGRESS");
                        break;
                    case GRB.Status.ITERATION_LIMIT:
                        System.Windows.Forms.MessageBox.Show("ITERATION_LIMIT");
                        break;
                    case GRB.Status.INTERRUPTED:
                        System.Windows.Forms.MessageBox.Show("INTERRUPTED");
                        break;
                    case GRB.Status.SUBOPTIMAL:
                        System.Windows.Forms.MessageBox.Show("SUBOPTIMAL");
                        break;
                    case GRB.Status.SOLUTION_LIMIT:
                        System.Windows.Forms.MessageBox.Show("SOLUTION_LIMIT");
                        break;
                    case GRB.Status.NODE_LIMIT:
                        System.Windows.Forms.MessageBox.Show("NODE_LIMIT");
                        break;                    
                    default:
                        System.Windows.Forms.MessageBox.Show("Something elese");
                        break;

                }
                for (int i = 0; i < 3; i++)
                {
                    if (model.Get(GRB.IntAttr.Status) == GRB.Status.OPTIMAL)
                    {
                        for (int j = 0; j < nU * nV; j++)
                        {
                            x[j, 2] = vars[j].Get(GRB.DoubleAttr.X);
                        }
                        for (int j = 0; j < nConstraints; j += 2)
                        {
                            double NORM = Math.Sqrt(bVars[j].Get(GRB.DoubleAttr.X) * bVars[j].Get(GRB.DoubleAttr.X) + bVars[j + 1].Get(GRB.DoubleAttr.X) * bVars[j + 1].Get(GRB.DoubleAttr.X));
                            if (NORM < 2)
                            {
                                lb[j] = 0;
                                ub[j] = 0;
                                lb[j + 1] = 0;
                                ub[j + 1] = 0;
                            }
                        }
                        model.Set(GRB.DoubleAttr.LB, bVars, lb);
                        model.Set(GRB.DoubleAttr.UB, bVars, ub);
                        model.Optimize();
                    }
                    else break;
                }
                model.Dispose();
                env.Dispose();

            }
            catch (GRBException e)
            {
                AddRuntimeMessage(Grasshopper.Kernel.GH_RuntimeMessageLevel.Error, "GUROBI ERROR!" + e.ErrorCode + ". " + e.Message);
                return;
            }
            double[] residual2 = new double[nConstraints];
            for (int i = 0; i < nConstraints; i++)
            {
                double val = 0;
                for (int j = 0; j < nU * nV; j++)
                {
                    val+= jacobian[i, j] * x[j, 2];
                }
                residual2[i]=val;
            }
            x2Nurbs(x, airyNurbs);
            createAiryPoints(airyNurbs);
            myMasonry.setupNodesFromList(x);
            myMasonry.computeAiryFunction();
            myMasonry.getResidual(residual);
            System.Windows.Forms.MessageBox.Show(residual.Norm().ToString());
            //Solve
            myMasonry.computeEigenVectors();
            Nurbs2x(xyNurbs, x);
            myMasonry.setupNodesFromList(x);
            myMasonry.computeGlobalCoord();
            ShoNS.Array.SparseDoubleArray hess = new SparseDoubleArray(nU * nV * 3, nU * nV * 3);
            myMasonry.computeHessian();
            myMasonry.getHessian(hess);
            Nurbs2x(outputNurbs, x);
            int nParticles = nU * nV;
            //Current configuration
            var origX = DoubleArray.Zeros(nParticles*3, 1);
            for (int i = 0; i < nParticles; i++)
            {
                origX[i * 3 + 0, 0] = x[i, 0];
                origX[i * 3 + 1, 0] = x[i, 1];
                origX[i * 3 + 2, 0] = x[i, 2];
            }
            List<int> shift = new List<int>();
            int T1 = 0;
            int T2 = 0;
            for (int i = 0; i < nParticles; i++)
            {
                shift.Add(i);
            }
            T1 = (nParticles - boundaryIndex.Count()) * 3 - 1;
            T2 = nParticles * 3 - 1;
            int C1 = 0;
            int C2 = nParticles - boundaryIndex.Count();
            for (int i = 0; i < nParticles; i++)
            {
                if (isBoundary(i))
                {
                    shift[i] = C2;
                    C2++;
                }
                else
                {
                    shift[i] = C1;
                    C1++;
                }
            }
            var shiftArray = new SparseDoubleArray(nParticles * 3, nParticles * 3);
            for (int i = 0; i < nParticles; i++)
            {
                shiftArray[i * 3, shift[i] * 3] = 1;
                shiftArray[i * 3 + 1, shift[i] * 3 + 1] = 1;
                shiftArray[i * 3 + 2, shift[i] * 3 + 2] = 1;
            }
            var ED = shiftArray.T.Multiply(hess) as SparseDoubleArray;
            ED = ED.Multiply(shiftArray) as SparseDoubleArray;
            var slice1 = new SparseDoubleArray(T1 + 1, T2 + 1);
            var slice2 = new SparseDoubleArray(T2 + 1, T2 - T1);
            for (int i = 0; i < T1 + 1; i++)
            {
                slice1[i, i] = 1;
            }
            for (int i = 0; i < T2 - T1; i++)
            {
                slice2[i + T1 + 1, i] = 1;
            }
            var DIB = (slice1.Multiply(ED) as SparseDoubleArray).Multiply(slice2) as SparseDoubleArray;
            var DII = (slice1.Multiply(ED) as SparseDoubleArray).Multiply(slice1.T) as SparseDoubleArray;
            var solver = new SparseLU(DII);
            origX = shiftArray.T * origX;
            var fixX = origX.GetSlice(T1 + 1, T2, 0, 0);
            var B = -DIB * fixX;
            var force = DoubleArray.Zeros(nParticles * 3, 1);
            for(int i=0;i<nParticles;i++)
            {
                force[i * 3 + 2, 0] = 0;//-0.3;
            }
            force = (shiftArray.T*force).GetSlice(0, T1, 0, 0);
            B = B + force;
            var dx = solver.Solve(B);

            var ret = DoubleArray.Zeros(nParticles *3, 1);
            for (int i = 0; i < T1 + 1; i++)
            {
                ret[i, 0] = dx[i, 0];
            }
            for (int i = T1 + 1; i <= T2; i++)
            {
                ret[i, 0] = fixX[i - T1 - 1, 0];
            }
            var xx = shiftArray * ret;
            
            for (int i = 0; i < nParticles; i++)
            {
                x[i, 0] = xx[i * 3, 0];
                x[i, 1] = xx[i * 3 + 1, 0];
                x[i, 2] = xx[i * 3 + 2, 0];
            }
            x2Nurbs(x, outputNurbs);
        }
        public override void DrawViewportWires(Grasshopper.Kernel.IGH_PreviewArgs args)
        {
            if (Hidden)
            {
                return;
            }

            args.Display.DrawSurface(airyNurbs, System.Drawing.Color.Blue, 1);
            args.Display.DrawSurface(outputNurbs, System.Drawing.Color.Red, 1);
            args.Display.DrawSurface(xyNurbs, System.Drawing.Color.Green, 1);
            double[][] vec = new double[2][] { new double[3], new double[3] };
            double[] val = new double[2];
            double[] node=null;
            double S = 1000;
            Nurbs2x(xyNurbs, x);
            myMasonry.setupNodesFromList(x);
            myMasonry.computeGlobalCoord();

            foreach (var e in myMasonry.elemList)
            {
                for(int i=0;i<e.nIntPoint;i++)
                {
                    e.getEigenVectors(vec, val, i);
                    node=e.getIntPoint(i);
                    System.Drawing.Color color;
                    double S1 = S * val[0];
                    double S2 = S * val[1];
                    if (val[0] > 0)
                    { color = System.Drawing.Color.Cyan; }
                    else { color = System.Drawing.Color.Magenta; }
                    args.Display.DrawLine(new Rhino.Geometry.Point3d(node[0], node[1], node[2]), new Rhino.Geometry.Point3d(node[0] + vec[0][0] * S1, node[1] + vec[0][1] * S1, node[2] + vec[0][2] * S1), color, 1);
                    args.Display.DrawLine(new Rhino.Geometry.Point3d(node[0], node[1], node[2]), new Rhino.Geometry.Point3d(node[0] - vec[0][0] * S1, node[1] - vec[0][1] * S1, node[2] - vec[0][2] * S1), color, 1);
                    if (val[1] > 0)
                    { color = System.Drawing.Color.Cyan; }
                    else { color = System.Drawing.Color.Magenta; }
                    args.Display.DrawLine(new Rhino.Geometry.Point3d(node[0], node[1], node[2]), new Rhino.Geometry.Point3d(node[0] + vec[1][0] * S2, node[1] + vec[1][1] * S2, node[2] + vec[1][2] * S2), color, 1);
                    args.Display.DrawLine(new Rhino.Geometry.Point3d(node[0], node[1], node[2]), new Rhino.Geometry.Point3d(node[0] - vec[1][0] * S2, node[1] - vec[1][1] * S2, node[2] - vec[1][2] * S2), color, 1);
                }
                for(int i=0;i<e.nBIntPoint;i++)
                {
                    e.getBEigenVectors(vec, val, i);
                    node = e.getBIntPoint(i);
                    System.Drawing.Color color;
                    double S1 = S * val[0];
                    double S2 = S * val[1];
                    if (val[0] > 0)
                    { color = System.Drawing.Color.Cyan; }
                    else { color = System.Drawing.Color.Magenta; }
                    args.Display.DrawLine(new Rhino.Geometry.Point3d(node[0], node[1], node[2]), new Rhino.Geometry.Point3d(node[0] + vec[0][0] * S1, node[1] + vec[0][1] * S1, node[2] + vec[0][2] * S1), color, 1);
                    args.Display.DrawLine(new Rhino.Geometry.Point3d(node[0], node[1], node[2]), new Rhino.Geometry.Point3d(node[0] - vec[0][0] * S1, node[1] - vec[0][1] * S1, node[2] - vec[0][2] * S1), color, 1);
                    if (val[1] > 0)
                    { color = System.Drawing.Color.Cyan; }
                    else { color = System.Drawing.Color.Magenta; }
                    args.Display.DrawLine(new Rhino.Geometry.Point3d(node[0], node[1], node[2]), new Rhino.Geometry.Point3d(node[0] + vec[1][0] * S2, node[1] + vec[1][1] * S2, node[2] + vec[1][2] * S2), color, 1);
                    args.Display.DrawLine(new Rhino.Geometry.Point3d(node[0], node[1], node[2]), new Rhino.Geometry.Point3d(node[0] - vec[1][0] * S2, node[1] - vec[1][1] * S2, node[2] - vec[1][2] * S2), color, 1);

                }
            }
            List<Rhino.Geometry.Point3d> xyP = new List<Point3d>();
            //List<Rhino.Geometry.Point3d> airyP = new List<Point3d>();
            List<Rhino.Geometry.Point3d> outputP = new List<Point3d>();
            Nurbs2x(xyNurbs, x);
            myMasonry.setupNodesFromList(x);
            myMasonry.computeGlobalCoord();
            foreach (var e in myMasonry.elemList)
            {
                for (int i = 0; i < e.nIntPoint; i++)
                {
                    var d = e.getIntPoint(i);
                    var d2 = new Rhino.Geometry.Point3d(d[0], d[1], d[2]);
                    xyP.Add(d2);
                    args.Display.DrawPoint(d2, Rhino.Display.PointStyle.ControlPoint, 2, System.Drawing.Color.Red);
                }
                for (int i = 0; i < e.nBIntPoint; i++)
                {
                    var d = e.getBIntPoint(i);
                    var d2 = new Rhino.Geometry.Point3d(d[0], d[1], d[2]);
                    args.Display.DrawPoint(d2, Rhino.Display.PointStyle.X, 2, System.Drawing.Color.Red);
                }
                
            }
            /*Nurbs2x(airyNurbs, x);
            myMasonry.computeAll(x);
            foreach (var e in myMasonry.elemList)
            {
                for (int i = 0; i < e.nIntPoint; i++)
                {
                    var d = e.getIntPoint(i);
                    var d2 = new Rhino.Geometry.Point3d(d[0], d[1], d[2]);
                    airyP.Add(d2);
                    args.Display.DrawPoint(d2, Rhino.Display.PointStyle.ControlPoint, 2, System.Drawing.Color.Red);
                }
            }*/
            Nurbs2x(outputNurbs, x);
            myMasonry.setupNodesFromList(x);
            myMasonry.computeGlobalCoord();
            foreach (var e in myMasonry.elemList)
            {
                for (int i = 0; i < e.nIntPoint; i++)
                {
                    var d = e.getIntPoint(i);
                    var d2 = new Rhino.Geometry.Point3d(d[0], d[1], d[2]);
                    outputP.Add(d2);
                    args.Display.DrawPoint(d2, Rhino.Display.PointStyle.ControlPoint, 2, System.Drawing.Color.Red);
                }
            }
            /*foreach (var i in xyP.Zip(outputP, (f, s) => new Rhino.Geometry.Line(f, s)))
            {
                args.Display.DrawLine(i, System.Drawing.Color.LightGreen);
            }*/
            base.DrawViewportWires(args);
        }
        protected override void SolveInstance(Grasshopper.Kernel.IGH_DataAccess DA)
        {
            Object inputGeometry = null;
            if (!DA.GetData(0, ref inputGeometry)) { return; }
            if (inputGeometry is Grasshopper.Kernel.Types.GH_Surface)
            {
                inputNurbs = (inputGeometry as Grasshopper.Kernel.Types.GH_Surface).Value.Surfaces[0].ToNurbsSurface();
                initialized=initializeNurbs(inputNurbs);
            }
            __update = true;
        }
        void deleteFixedPoints()
        {
            if (fixedPointGuids != null)
            {
                foreach (Guid g in fixedPointGuids)
                {
                    Rhino.RhinoDoc.ActiveDoc.Objects.Delete(g, true);
                }
            }
            fixedPointGuids = null;

        }
        void createAiryPoints(Rhino.Geometry.NurbsSurface S)
        {
            if (fixedPointGuids != null)
            {
                deleteFixedPoints();
            }
            fixedPointGuids = new List<Guid>();
            for (int i = 0; i < nV; i++)
            {
                for (int j = 0; j < nU; j++)
                {
                    fixedPointGuids.Add(Rhino.RhinoDoc.ActiveDoc.Objects.AddPoint(S.Points.GetControlPoint(j,i).Location));
                }
            }
        }
        void createAiryInitialSurface(Rhino.Geometry.NurbsSurface S)
        {
            //find the center point
            for (int i = 0; i < nV; i++)
            {
                for (int j = 0; j < nU; j++)
                {
                    var P = S.Points.GetControlPoint(j, i);
                    double Z = (i - ((nV-1) / 2d)) * (i - ((nV-1) / 2d)) + (j - ((nU-1) / 2d)) * (j - ((nU-1) / 2d));
                    airyNurbs.Points.SetControlPoint(j, i, new Rhino.Geometry.ControlPoint(P.Location.X, P.Location.Y, -Z));
                }
            }
        }
        void x2Nurbs(double[,] _x,Rhino.Geometry.NurbsSurface S)
        {
            for (int i = 0; i < nV; i++)
            {
                for (int j = 0; j < nU; j++)
                {
                    var P=new Rhino.Geometry.Point3d(x[(i * nU + j) , 0],x[(i * nU + j) , 1] ,x[(i * nU + j) ,2] );
                    S.Points.SetControlPoint(j, i, P);
                }
            }

        }
        void Nurbs2x(Rhino.Geometry.NurbsSurface S,double[,] _x)
        {
            for (int i = 0; i < nV; i++)
            {
                for (int j = 0; j < nU; j++)
                {
                    x[(i*nU+j),0]=S.Points.GetControlPoint(j, i).Location.X;
                    x[(i*nU+j),1]=S.Points.GetControlPoint(j, i).Location.Y;
                    x[(i*nU+j),2]=S.Points.GetControlPoint(j, i).Location.Z;
                }
            }
        }
        bool initializeNurbs(Rhino.Geometry.NurbsSurface S)
        {
            var X = Rhino.Geometry.Transform.Translation(250, 0, 0);
            alternativeMesh = new Rhino.Geometry.Mesh();
            nU = S.Points.CountU;
            nV = S.Points.CountV;

            for (int j = 0; j < nV; j++)
            {
                for (int i = 0; i < nU; i++)
                {
                    Rhino.Geometry.ControlPoint cPoint = S.Points.GetControlPoint(i, j);
                    alternativeMesh.Vertices.Add(cPoint.Location);
                }
            }
            for (int j = 0; j < nV - 1; j++)
            {
                for (int i = 0; i < nU - 1; i++)
                {
                    int D1 = j * nU + i;
                    int D2 = j * nU + i + 1;
                    int D3 = (j + 1) * nU + i + 1;
                    int D4 = (j + 1) * nU + i;
                    alternativeMesh.Faces.AddFace(new Rhino.Geometry.MeshFace(D1, D2, D3, D4));
                }
            }
            //Get Boundary
            boundaryIndex = new List<int>();
            /*Rhino.Geometry.Polyline[] boundary = alternativeMesh.GetNakedEdges();
            
            foreach (var p in boundary)
            {
                for (int i = 0; i < p.Count - 1; i++)
                {
                    var f = alternativeMesh.Vertices.Select((e,index)=>new{e,index}).Where(e=>e.e==p[i]).Select(e=>e.index).First();
                    boundaryIndex.Add(f);
                }
            }*/
            boundaryIndex.Add(0);
            boundaryIndex.Add(nU-1);
            boundaryIndex.Add(nU * nV - 1-nU+1);
            boundaryIndex.Add(nU * nV - 1);
            
            airyNurbs = S.Duplicate() as Rhino.Geometry.NurbsSurface;
            outputNurbs = S.Duplicate() as Rhino.Geometry.NurbsSurface;
            xyNurbs = S.Duplicate() as Rhino.Geometry.NurbsSurface;

            airyNurbs.Transform(X);
            outputNurbs.Transform(X);
            xyNurbs.Transform(X);

            for (int i = 0; i < nV; i++)
            {
                for (int j = 0; j < nU; j++)
                {
                    var P = xyNurbs.Points.GetControlPoint(j, i);
                    airyNurbs.Points.SetControlPoint(j, i, new Rhino.Geometry.ControlPoint(P.Location.X, P.Location.Y, 0));
                    xyNurbs.Points.SetControlPoint(j, i, new Rhino.Geometry.ControlPoint(P.Location.X, P.Location.Y, 0));
                }
            }
            createAiryInitialSurface(airyNurbs);
            createAiryPoints(airyNurbs);
            createNurbsElements(inputNurbs);

            x = new double[nU * nV , 3];
            Nurbs2x(airyNurbs,x);
            myMasonry.setupNodesFromList(x);
            myMasonry.computeAiryFunction();
            myMasonry.computeEigenVectors();
            Nurbs2x(xyNurbs, x);
            myMasonry.setupNodesFromList(x);
            myMasonry.computeGlobalCoord();

            return true;
        }
        void createNurbsElements(Rhino.Geometry.NurbsSurface S)
        {
            double[] uKnot;
            double[] vKnot;

            int N = nU * nV;
            int uDim = S.OrderU;
            int vDim = S.OrderV;
            int uDdim = S.OrderU - 1;
            int vDdim = S.OrderV - 1;


            uKnot = new double[nU - uDdim + 1 + uDdim * 2];
            vKnot = new double[nV - vDdim + 1 + vDdim * 2];
            for (int i = 0; i < uDdim; i++)
            {
                uKnot[i] = 0;
            }
            for (int i = 0; i < vDdim; i++)
            {
                vKnot[i] = 0;
            }
            for (int i = 0; i < nU - uDdim + 1; i++)
            {
                uKnot[i + uDdim] = i;
            }
            for (int i = 0; i < nV - vDdim + 1; i++)
            {
                vKnot[i + vDdim] = i;
            }
            for (int i = 0; i < uDdim; i++)
            {
                uKnot[i + nU + 1] = nU - uDdim;
            }
            for (int i = 0; i < vDdim; i++)
            {
                vKnot[i + nV + 1] = nV - vDdim;
            }
            myMasonry.elemList.Clear();
            myMasonry.Clear();
            for (int j = 1; j < nV - vDdim + 1; j++)
            {
                for (int i = 1; i < nU - uDdim + 1; i++)
                {
                    int[] index = new int[uDim * vDim];
                    for (int k = 0; k < vDim; k++)
                    {
                        for (int l = 0; l < uDim; l++)
                        {
                            index[k * uDim + l] = (j - 1 + k) * nU + i - 1 + l;
                        }
                    }
                    //judge if on the border
                    Minilla3D.Elements.nurbsElement.border border=Minilla3D.Elements.nurbsElement.border.None;
                    if (j == 1)
                    {
                        border = border | Minilla3D.Elements.nurbsElement.border.Top;
                    }
                    if (i == 1)
                    {
                        border = border | Minilla3D.Elements.nurbsElement.border.Left;
                    }
                    if (j == nV - vDdim)
                    {
                        border = border | Minilla3D.Elements.nurbsElement.border.Bottom;
                    }
                    if (i == nU - uDdim)
                    {
                        border = border | Minilla3D.Elements.nurbsElement.border.Right;
                    }
                    myMasonry.elemList.Add(new Minilla3D.Elements.nurbsElement(uDim, vDim, 2, index, i, j, uKnot, vKnot, border));
                }
            }
        }
    }
    /*
    public partial class Masonry : Grasshopper.Kernel.GH_Component
    {

        static Mothra()
        {
            var dir = Grasshopper.Folders.DefaultAssemblyFolder;

            string assemblyDir = Path.GetDirectoryName(dir);
            
            if (File.Exists(Path.Combine(assemblyDir, "Minilla.proxy.dll"))
                || !File.Exists(Path.Combine(assemblyDir, "Minilla.x86.dll"))
                || !File.Exists(Path.Combine(assemblyDir, "Minilla.x64.dll")))
            {
                throw new InvalidOperationException("Found Kapybara.dll which cannot exist. "
                    + "Must instead have Kapybara.x86.dll and Kapybara.x64.dll. Check your build settings.");
            }

            AppDomain.CurrentDomain.AssemblyResolve += (_, e) =>
            {
                if (e.Name.StartsWith("Minilla.proxy", StringComparison.OrdinalIgnoreCase))
                {
                    string fileName = Path.Combine(assemblyDir,
                        string.Format("Minilla.{0}.dll", (IntPtr.Size == 4) ? "x86" : "x64"));
                    return Assembly.LoadFile(fileName);
                }
                return null;
            };
        }


        Func<double, double> Drift0 = (v) => { return v / 20d + 0.95; };
        static float Shifting = 500;
        mikity.visualize.FigureUI3 controlPanel;
            

        public Mothra()
            : base("Mothra", "Mothra", "Mothra", "Kapybara3D", "Computation")
        {
        }
        ~Mothra()
        {
            keyboardHook.Uninstall();
        }
        RamGecTools.MouseHook mouseHook = new RamGecTools.MouseHook();
        RamGecTools.KeyboardHook keyboardHook = new RamGecTools.KeyboardHook();
        bool reset = false;
        void activate()
        {
            full.activate();
        }
        void deactivate()
        {
            full.deactivate();
        }
        bool keyboardHook_KeyUp(RamGecTools.KeyboardHook.VKeys key)
        {
            if (key == RamGecTools.KeyboardHook.VKeys.KEY_R)
            {
                if (_RF)
                {
                    _RF = false;
                    full.offRF();
                }
                else
                {
                    _RF = true;
                    full.onRF();
                }
                return true;
            }
            if (key == RamGecTools.KeyboardHook.VKeys.ESCAPE)
            {
                full.resetGo();
                _go = false;
                Clock = -1;
                timer.Enabled = false;
                isInitialized = false;
                reset = true;
                vel.FillValue(0);
                vel2.FillValue(0);
                ExpireSolution(true);
                return true;
            }
            if (key == RamGecTools.KeyboardHook.VKeys.KEY_G)
            {
                if (_go)
                {
                    full.pauseGo();
                    _go = false;
                    timer.Enabled = false;
                }
                else
                {
                    full.onGo();
                    _go = true;
                    timer.Enabled = true;
                    setupMaterial(mat);
                }
                return true;
            }
            if (key == RamGecTools.KeyboardHook.VKeys.KEY_I)
            {
                if (_intPoint)
                {
                    _intPoint = false;
                }
                else
                {
                    _intPoint = true;
                }
                reset = true;
                this.ExpireSolution(true);
                return true;
            }
            if (key == RamGecTools.KeyboardHook.VKeys.KEY_B)
            {
                if (_base)
                {
                    _base = false;
                }
                else
                {
                    _base = true;
                }
                reset = true;
                this.ExpireSolution(true);
                return true;
            }
            if (key == RamGecTools.KeyboardHook.VKeys.KEY_D)
            {
                if (stress)
                {
                    stress = false;
                }
                else
                {
                    stress = true;
                }
                reset = true;
                this.ExpireSolution(true);
                return true;
            }
            return false;
        }

        bool keyboardHook_KeyDown(RamGecTools.KeyboardHook.VKeys key)
        {
            if (key == RamGecTools.KeyboardHook.VKeys.KEY_M)
            {
                return true;
            }
            if (key == RamGecTools.KeyboardHook.VKeys.KEY_R)
            {
                return true;
            }
            if (key == RamGecTools.KeyboardHook.VKeys.KEY_F)
            {
                return true;
            }
            if (key == RamGecTools.KeyboardHook.VKeys.KEY_A)
            {
                return true;
            }
            if (key == RamGecTools.KeyboardHook.VKeys.ESCAPE)
            {
                return true;
            }
            if (key == RamGecTools.KeyboardHook.VKeys.KEY_G)
            {
                return true;
            }
            if (key == RamGecTools.KeyboardHook.VKeys.KEY_I)
            {
                return true;
            }
            if (key == RamGecTools.KeyboardHook.VKeys.KEY_B)
            {
                return true;
            }
            if (key == RamGecTools.KeyboardHook.VKeys.KEY_D)
            {
                return true;
            }
            return false; 
        }
        void timer_Tick(object sender, EventArgs e)
        {
            this.ExpireSolution(true);
        }
        mikity.visualize.fullScreen full;
        System.Windows.Forms.Timer timer;
        System.Windows.Forms.ToolStripMenuItem __gm1, __vf1, __vf2, __vf3, __vf4;
        System.Windows.Forms.ToolStripMenuItem __bm1, __bm2, __bm3, __gm2;
        System.Windows.Forms.ToolStripMenuItem __ig1, __ig21, __ig22, __ig23, __ig3, __ig5;
        public override void AddedToDocument(Grasshopper.Kernel.GH_Document document)
        {
            base.AddedToDocument(document);
            //test();
            Rhino.RhinoDoc.ReplaceRhinoObject += RhinoDoc_ReplaceRhinoObject;
            timer = new System.Windows.Forms.Timer();
            timer.Tick += timer_Tick;
            timer.Enabled = false;
            timer.Interval = 1;

            // register evens
            keyboardHook.KeyDown = new RamGecTools.KeyboardHook.KeyboardHookCallback(keyboardHook_KeyDown);
            keyboardHook.KeyUp = new RamGecTools.KeyboardHook.KeyboardHookCallback(keyboardHook_KeyUp);
            keyboardHook._activate = new RamGecTools.KeyboardHook.activate(activate);
            keyboardHook._deactivate = new RamGecTools.KeyboardHook.deactivate(deactivate);
            keyboardHook.Install();

            mat = materialChoice.DHM;
            _material = new Minilla3D.Materials.harmonicMaterial();
            type = modelType.membrane;
            _isFixedBoundary = true;
            controlPanel = new mikity.visualize.FigureUI3();
            full = new mikity.visualize.fullScreen();
            _param.wA = 50;
            _param.wL = 50;
            _param.wC = 50;
            _param.wT = 0;
            _param.wT2 = 0;
            _param.Neo = 100;
            _param.Alpha = 50;
            _param.hmAlpha = 50;
            _param.dcmDist = 50;
            __rebuildControlPanel();
            controlPanel.Show();
            full.deactivate();
            full.Show();
            full.resetGo();
            full.offRF();
            full._selectMaterial = (s) => selectMaterial(s);

        }
        public override void RemovedFromDocument(Grasshopper.Kernel.GH_Document document)
        {
            base.RemovedFromDocument(document);
            keyboardHook.Uninstall();
            keyboardHook.KeyDown -= new RamGecTools.KeyboardHook.KeyboardHookCallback(keyboardHook_KeyDown);
            keyboardHook.KeyUp -= new RamGecTools.KeyboardHook.KeyboardHookCallback(keyboardHook_KeyUp);
            if (full != null)
            {
                full.Close();
                full = null;
            }
        }
        public override void DocumentContextChanged(Grasshopper.Kernel.GH_Document document, Grasshopper.Kernel.GH_DocumentContext context)
        {
            base.DocumentContextChanged(document, context);
            if (context == Grasshopper.Kernel.GH_DocumentContext.Unloaded)
            {

                keyboardHook.Uninstall();
                keyboardHook.KeyDown -= new RamGecTools.KeyboardHook.KeyboardHookCallback(keyboardHook_KeyDown);
                keyboardHook.KeyUp -= new RamGecTools.KeyboardHook.KeyboardHookCallback(keyboardHook_KeyUp);
                if (full != null)
                {
                    full.Close();
                    full = null;
                }
            }
        }
        public struct parameter
        {
            public int wL, wA, wC, Alpha, Neo, hmAlpha, dcmDist, wT, wT2;
        }
        delegate void updateParam();
        updateParam __UPDATE = () => { };
        private parameter _param;
        private double refArea = 0, area = 0;
        enum modelType
        {
            wire, membrane
        }
        enum materialChoice
        {
            DCM, L2G, CLM, SCP, DHM, SV, NH, Cotan, MKM
        }
        private materialChoice mat;
        private modelType type;

        private void __rebuildControlPanel()
        {
            controlPanel.clearSliders();
            switch (mat)
            {
                case materialChoice.DCM:
                    controlPanel.addSlider(0, 2, 100, _param.dcmDist, "dist");
                    controlPanel.listSlider[0].Converter = (v) => v * 0.02;
                    __UPDATE = () =>
                    {
                        if (mat == materialChoice.DCM)
                        {
                            Minilla3D.Materials.leastSquaresMaterial ls = _material as Minilla3D.Materials.leastSquaresMaterial;
                        }
                    };
                    break;
                case materialChoice.SCP:
                    controlPanel.addSlider(0, 2, 99, _param.hmAlpha, "Alpha");
                    controlPanel.listSlider[0].Converter = (v) => v * 0.02;
                    __UPDATE = () =>
                    {
                        if (mat == materialChoice.SCP)
                        {
                            Minilla3D.Materials.conformalMaterial cm = _material as Minilla3D.Materials.conformalMaterial;
                            cm.refArea = refArea * controlPanel.listSlider[0].value;
                            _param.hmAlpha = controlPanel.listSlider[0].originalValue;
                        }
                    }; break;
                case materialChoice.DHM:
                    __UPDATE = () => { };

                    break;
                case materialChoice.CLM:
                    controlPanel.addSlider(0, 2, 100, _param.wL, "wL");
                    controlPanel.addSlider(0, 2, 100, _param.wA, "wA");
                    controlPanel.addSlider(0, 2, 100, _param.wC, "wC");
                    controlPanel.addSlider(0, 2, 100, _param.wT, "t");
                    controlPanel.listSlider[0].Converter = (v) => v * 0.01;
                    controlPanel.listSlider[1].Converter = (v) => v * 0.01;
                    controlPanel.listSlider[2].Converter = (v) => v * 0.01;
                    controlPanel.listSlider[3].Converter = (v) => v / 50d + 1.0;
                    __UPDATE = () =>
                    {
                        if (mat == materialChoice.CLM)
                        {
                            Minilla3D.Materials.clarenzMaterial cl = _material as Minilla3D.Materials.clarenzMaterial;
                            cl.WL = controlPanel.listSlider[0].value;
                            cl.WA = controlPanel.listSlider[1].value;
                            cl.WC = controlPanel.listSlider[2].value;
                            cl.T = controlPanel.listSlider[3].value;
                            _param.wL = controlPanel.listSlider[0].originalValue;
                            _param.wA = controlPanel.listSlider[1].originalValue;
                            _param.wC = controlPanel.listSlider[2].originalValue;
                            _param.wT = controlPanel.listSlider[3].originalValue;
                        }
                    };
                    break;
                case materialChoice.MKM:
                    controlPanel.addSlider(0, 2, 100, _param.wT2, "t");
                    controlPanel.listSlider[0].Converter = (v) => v / 50d + 1.0;
                    __UPDATE = () =>
                    {
                        if (mat == materialChoice.MKM)
                        {
                            Minilla3D.Materials.mikityMaterial cl = _material as Minilla3D.Materials.mikityMaterial;
                            cl.T = controlPanel.listSlider[0].value;
                            _param.wT2 = controlPanel.listSlider[0].originalValue;
                        }
                    };
                    break;
                case materialChoice.NH:
                    controlPanel.addSlider(0, 2, 500, _param.Neo, "K");
                    controlPanel.listSlider[0].Converter = (v) => v * 0.1;
                    __UPDATE = () =>
                    {
                        if (mat == materialChoice.NH)
                        {
                            Minilla3D.Materials.neoHookeanMaterial nH = _material as Minilla3D.Materials.neoHookeanMaterial;
                            nH.mu1 = 0.1;
                            nH.K = controlPanel.listSlider[0].value;
                            _param.Neo = controlPanel.listSlider[0].originalValue;
                        }
                    };
                    break;
                case materialChoice.L2G:
                    __UPDATE = () => { };
                    break;
                case materialChoice.SV:
                    controlPanel.addSlider(1, 2, 99, _param.Alpha, "Alpha");
                    controlPanel.listSlider[0].Converter = (v) => v * 0.01;
                    __UPDATE = () =>
                    {
                        if (mat == materialChoice.SV)
                        {
                            Minilla3D.Materials.stVenantMaterial sv = _material as Minilla3D.Materials.stVenantMaterial;
                            sv.Alpha = controlPanel.listSlider[0].value;
                            _param.Alpha = controlPanel.listSlider[0].originalValue;
                        }
                    }; break;
                case materialChoice.Cotan:
                    __UPDATE = () => { };
                    break;

            }
        }
        enum initialGuess { prj, dhm1, dhm2, dhm3, dcm, scpl, scps };
        initialGuess InitialGuess = initialGuess.prj;
        private bool assignChecked(bool flag, System.Windows.Forms.ToolStripMenuItem m)
        {
            if (flag)
            {
                m.CheckState = System.Windows.Forms.CheckState.Checked;
            }
            else
            {
                m.CheckState = System.Windows.Forms.CheckState.Unchecked;
            }
            return flag;
        }
        public override void AppendAdditionalMenuItems(System.Windows.Forms.ToolStripDropDown menu)
        {
            base.AppendAdditionalMenuItems(menu);
            Menu_AppendSeparator(menu);
            __gm1 = Menu_AppendItem(menu, "Go?", Menu_GoClicked);
            __gm2 = Menu_AppendItem(menu, "Show Parameters", Menu_ShowParam);
            Menu_AppendSeparator(menu);
            System.Windows.Forms.ToolStripMenuItem tt = Menu_AppendItem(menu, "Initial guess");
            __ig1 = Menu_AppendItem(tt.DropDown, "Projection", Menu_InitialClicked);
            System.Windows.Forms.ToolStripMenuItem ttt = Menu_AppendItem(tt.DropDown, "DHM");
            __ig21 = Menu_AppendItem(ttt.DropDown, "DHM (Projection)", Menu_InitialClicked);
            __ig22 = Menu_AppendItem(ttt.DropDown, "DHM (Circle)", Menu_InitialClicked);
            __ig23 = Menu_AppendItem(ttt.DropDown, "DHM (Rectangle)", Menu_InitialClicked);
            __ig3 = Menu_AppendItem(tt.DropDown, "DCM", Menu_InitialClicked);
            //__ig4 = Menu_AppendItem(tt.DropDown, "SCPL", Menu_InitialClicked);
            __ig5 = Menu_AppendItem(tt.DropDown, "SCP", Menu_InitialClicked);
            tt = Menu_AppendItem(menu, "Optional features");
            __vf1 = Menu_AppendItem(tt.DropDown, "Show base vectors?", Menu_VisualClicked);
            __vf2 = Menu_AppendItem(tt.DropDown, "Show integrating points?", Menu_VisualClicked);
            __vf3 = Menu_AppendItem(tt.DropDown, "Show eigen vectors?", Menu_VisualClicked);
            __vf4 = Menu_AppendItem(tt.DropDown, "Show conformality?", Menu_VisualClicked);
            tt.DropDownItems.Add(__vf1);
            tt.DropDownItems.Add(__vf2);
            tt.DropDownItems.Add(__vf3);
            tt.DropDownItems.Add(__vf4);
            tt = Menu_AppendItem(menu, "Boundary Shapes");
            __bm1 = Menu_AppendItem(menu, "Boundary:Projection", Menu_BoundaryClicked);
            __bm2 = Menu_AppendItem(menu, "Boundary:Circle", Menu_BoundaryClicked);
            __bm3 = Menu_AppendItem(menu, "Boundary:Rectangle", Menu_BoundaryClicked);
            tt.DropDownItems.Add(__bm1);
            tt.DropDownItems.Add(__bm2);
            tt.DropDownItems.Add(__bm3);
            Menu_AppendSeparator(menu);
            assignChecked(InitialGuess == initialGuess.prj, __ig1);
            assignChecked(InitialGuess == initialGuess.dhm1, __ig21);
            assignChecked(InitialGuess == initialGuess.dhm2, __ig22);
            assignChecked(InitialGuess == initialGuess.dhm3, __ig23);
            assignChecked(InitialGuess == initialGuess.dcm, __ig3);
            assignChecked(InitialGuess == initialGuess.scps, __ig5);
            if (!assignChecked(_go, __gm1)) Clock = -1;

            assignChecked(controlPanel.Visibility == System.Windows.Visibility.Visible, __gm2);
            assignChecked(_base, __vf1);
            assignChecked(stress, __vf2);
            assignChecked(_intPoint, __vf3);
            assignChecked(_conformal, __vf4);

        }
        private string switches()
        {
            string output="";
            if (_drift1)
            {
                output += "Drift:1\n";
            }else if(_drift2)
            {
                output += "Drift:2\n";
            }else{
                output+="Drift:None\n";
            }

            return output;
        }
        private bool _go = false;
        private bool _base = false;
        private bool _conformal = false;
        private bool stress = false;
        private bool _intPoint = true;
        private bool _isFixedBoundary;
        private bool _drift1 = true, _drift2 = false;
        private bool _fixFlip = true, _RF = false;
        private int shift = 0;
        private Minilla3D.Materials.iMaterial _material = new Minilla3D.Materials.stVenantMaterial();
        private void Menu_ShowParam(Object sender, EventArgs e)
        {
            if (sender == __gm2)
            {
                if (controlPanel.Visibility == System.Windows.Visibility.Visible)
                {
                    controlPanel.Hide();
                }
                else
                {
                    controlPanel.Show();
                }
            }
        }

        private void Menu_BoundaryClicked(Object sender, EventArgs e)
        {
            System.Windows.Forms.ToolStripMenuItem _m = sender as System.Windows.Forms.ToolStripMenuItem;

            if (_m == __bm1)
            {
                projBoundary();
            }
            if (_m == __bm2)
            {
                circleBoundary();
            }
            if (_m == __bm3)
            {
                rectangleBoundary();
            }
        }
        private bool switchState(System.Windows.Forms.ToolStripMenuItem _m, ref bool flag)
        {
            if (_m.CheckState == System.Windows.Forms.CheckState.Checked)
            {
                _m.CheckState = System.Windows.Forms.CheckState.Unchecked;
                flag = false;
            }
            else if (_m.CheckState == System.Windows.Forms.CheckState.Unchecked)
            {
                _m.CheckState = System.Windows.Forms.CheckState.Checked;
                flag = true;
            }
            return flag;
        }
        private void Menu_VisualClicked(Object sender, EventArgs e)
        {
            System.Windows.Forms.ToolStripMenuItem _m = sender as System.Windows.Forms.ToolStripMenuItem;
            if (_m == __vf1) switchState(_m, ref _base);
            if (_m == __vf2) switchState(_m, ref _intPoint);
            if (_m == __vf3)
            {
                if (switchState(_m, ref stress)) { _conformal = false; }
            }
            if (_m == __vf4)
            {
                if (switchState(_m, ref _conformal)) { stress = false; }
            }
            reset = true;
            this.ExpireSolution(true);

        }

        private void Menu_GoClicked(Object sender, EventArgs e)
        {
            System.Windows.Forms.ToolStripMenuItem _m = sender as System.Windows.Forms.ToolStripMenuItem;
            if (_m == __gm1)
            {
                if (switchState(_m, ref _go))
                {
                    timer.Enabled = true;
                }
                else
                {
                    timer.Enabled = false;
                    Clock = -1;
                }
            }
        }

        public override bool Read(GH_IO.Serialization.GH_IReader reader)
        {
            _go = false;//default value;
            reader.TryGetBoolean("Go?", ref _go);
            _base = false;
            reader.TryGetBoolean("Show Base Vectors?", ref _base);
            stress = true;
            reader.TryGetBoolean("Show Eigen Vectors?", ref stress);
            _intPoint = true;
            reader.TryGetBoolean("Show Integrating Points?", ref _intPoint);
            return base.Read(reader);
        }
        public override bool Write(GH_IO.Serialization.GH_IWriter writer)
        {
            writer.SetBoolean("Go?", false);
            writer.SetBoolean("Show Base Vectors?", _base);
            writer.SetBoolean("Show Eigen Vectors?", stress);
            writer.SetBoolean("Show Integrating Points?", _intPoint);
            return base.Write(writer);
        }

        protected override void RegisterInputParams(Grasshopper.Kernel.GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Mesh/Surface", "m", "Mesh/Surface", Grasshopper.Kernel.GH_ParamAccess.item);
            pManager.AddIntegerParameter("D1", "D1", "D1", Grasshopper.Kernel.GH_ParamAccess.item, 0);
            pManager.AddIntegerParameter("D2", "D2", "D2", Grasshopper.Kernel.GH_ParamAccess.item, 0);
        }

        protected override void RegisterOutputParams(Grasshopper.Kernel.GH_Component.GH_OutputParamManager pManager)
        {
            //pManager.AddIntegerParameter("iteration", "out", "out", Grasshopper.Kernel.GH_ParamAccess.item);
            //pManager.AddNumberParameter("a/A", "out", "out", Grasshopper.Kernel.GH_ParamAccess.item);
            pManager.AddNumberParameter("maxC", "maxC", "maxC", Grasshopper.Kernel.GH_ParamAccess.item);
            pManager.AddNumberParameter("minC", "minC", "minC", Grasshopper.Kernel.GH_ParamAccess.item);
            pManager.AddNumberParameter("integral", "intC", "out4", Grasshopper.Kernel.GH_ParamAccess.item);
            //pManager.AddNumberParameter("area", "out5", "out5", Grasshopper.Kernel.GH_ParamAccess.item);
            //pManager.AddTextParameter("switches", "switches", "switches", Grasshopper.Kernel.GH_ParamAccess.item);
            //pManager.AddNumberParameter("damping", "damping", "damping", Grasshopper.Kernel.GH_ParamAccess.item);
        }
        public override void BakeGeometry(Rhino.RhinoDoc doc, Rhino.DocObjects.ObjectAttributes att, List<Guid> obj_ids)
        {
            if (this.BKGT != null)
            {
                this.BKGT(doc, att, obj_ids);
            }
        }
        public override void DrawViewportWires(Grasshopper.Kernel.IGH_PreviewArgs args)
        {
            if (this.DVPW != null)
            {
                this.DVPW(args);
            }
            base.DrawViewportWires(args);
        }
        List<Minilla3D.Elements.managedElement> elemList = new List<Minilla3D.Elements.managedElement>();
        List<Minilla3D.Elements.managedElement> elemList2 = new List<Minilla3D.Elements.managedElement>();
        int Clock = -1;
        Rhino.Geometry.Mesh m;
        Rhino.Geometry.Mesh triMesh = new Rhino.Geometry.Mesh();
        Rhino.Geometry.NurbsSurface ns;
        List<Rhino.Geometry.Point3d> iP = new List<Rhino.Geometry.Point3d>();
        List<Rhino.Geometry.Line> bV = new List<Rhino.Geometry.Line>();
        List<Rhino.Geometry.Line> eVT = new List<Rhino.Geometry.Line>();
        List<Rhino.Geometry.Line> reF = new List<Rhino.Geometry.Line>();
        List<Rhino.Geometry.Line> eVC = new List<Rhino.Geometry.Line>();
        List<Rhino.Geometry.Circle> cfm = new List<Rhino.Geometry.Circle>();
        int nParticles = 0;
        Minilla3D.Objects.generalSpring gS = new Minilla3D.Objects.generalSpring();
        Minilla3D.Objects.generalSpring cN = new Minilla3D.Objects.generalSpring();
        //double _wA = 1.0, _wL = 1.0, _wC = 1.0;
        int nElements = 0;
        double[,] pos, refPos, Re;
        DoubleArray vel, vel2,force, acc;
        List<double[]>[] star;
        SparseDoubleArray hessEd, hessA;
        SparseDoubleArray metric, invMetric;
        DrawViewPortWire DVPW = null;
        BakeGeometry BKGT = null;
        private void initialize(Rhino.Geometry.Mesh _m, Rhino.Geometry.Mesh _m2)
        {
            if (internalState == state.mesh)
            {
                elemList.Clear();
            }
            elemList2.Clear();
            triMesh.Vertices.Clear();
            triMesh.Faces.Clear();
            for (int i = 0; i < _m.Vertices.Count; i++)
            {
                triMesh.Vertices.Add(_m.Vertices[i]);
            }
            foreach (Rhino.Geometry.MeshFace F in _m.Faces)
            {
                if (F.IsQuad)
                {
                    if (internalState == state.mesh)
                    {
                        elemList.Add(new Minilla3D.Elements.I4D2(new int[4] { F.A, F.B, F.D, F.C }));
                    }
                    triMesh.Faces.AddFace(new Rhino.Geometry.MeshFace(F.A, F.B, F.C));
                    triMesh.Faces.AddFace(new Rhino.Geometry.MeshFace(F.A, F.C, F.D));
                }
                else if (F.IsTriangle)
                {
                    if (internalState == state.mesh)
                    {
                        elemList.Add(new Minilla3D.Elements.S3D2(new int[3] { F.A, F.B, F.C }));
                    }
                    faces.Add(new face(F.A, F.B, F.C));
                    triMesh.Faces.AddFace(new Rhino.Geometry.MeshFace(F.A, F.B, F.C));
                    elemList2.AddRange(faces[faces.Count - 1].edges);
                }
            }
            nElements = elemList.Count();
            nParticles = _m.Vertices.Count();
            gS.Clear();
            gS.AddRange(elemList);
            gS.initialize(nParticles);
            cN.Clear();
            foreach (face f in faces)
            {
                cN.AddRange(f.edges);
            }
            cN.initialize(nParticles);
            force = DoubleArray.Zeros(nParticles * 3);
            acc = DoubleArray.Zeros(nParticles * 3);
            Re = new double[nParticles, 3];
            hessEd = new SparseDoubleArray(nParticles * 2, nParticles * 2);
            hessA = new SparseDoubleArray(nParticles * 2, nParticles * 2);

            vel = DoubleArray.Zeros(nParticles * 3);
            vel2 = DoubleArray.Zeros(nParticles * 3);
            invMetric = SparseDoubleArray.Zeros(nParticles * 3, nParticles * 3);
            metric = SparseDoubleArray.Zeros(nParticles * 3, nParticles * 3);
            pos = new double[nParticles, 3];
            refPos = new double[nParticles, 3];
            star = new List<double[]>[nParticles];
            for (int i = 0; i < nParticles; i++)
            {
                star[i] = new List<double[]>();
            }
            for (int i = 0; i < nParticles; i++)
            {
                pos[i, 0] = _m2.Vertices[i].X;
                pos[i, 1] = _m2.Vertices[i].Y;
                pos[i, 2] = _m2.Vertices[i].Z;
                refPos[i, 0] = _m.Vertices[i].X;
                refPos[i, 1] = _m.Vertices[i].Y;
                refPos[i, 2] = _m.Vertices[i].Z;
            }
            getBoundary(_m2);
        }
        struct boundaryVertex
        {
            public int index;
            public double param;
        }



        

        private boundaryVertex[] __boundary;
        public void getBoundary(Rhino.Geometry.Mesh _m)
        {
            Rhino.Geometry.Polyline boundary = _m.GetNakedEdges()[0];

            double[] boundaryParam = new double[boundary.Count];
            double totalLength = 0;
            for (int i = 0; i < boundary.Count - 1; i++)
            {
                totalLength += (boundary[i + 1] - boundary[i]).Length;
            }
            double currentLength = 0;
            boundaryParam[0] = 0;
            for (int i = 1; i < boundary.Count; i++)
            {
                currentLength += (boundary[i] - boundary[i - 1]).Length;
                boundaryParam[i] = currentLength;
            }
            __boundary = new boundaryVertex[boundary.Count - 1];
            for (int i = 0; i < boundary.Count - 1; i++)
            {
                for (int j = 0; j < _m.Vertices.Count; j++)
                {
                    if (boundary[i] == _m.Vertices[j])
                    {
                        __boundary[i].index = j;
                        __boundary[i].param = boundaryParam[i];
                    }
                }
            }
        }
        bool isInitialized = false;
        private void computeInitialGuess(double[,] rX, double[,] x, bool flag)
        {
            if (isInitialized)
            {
                if (InitialGuess == initialGuess.prj)
                {
                    initialShapePrj(rX, x);
                }
                if (InitialGuess == initialGuess.dhm1)
                {
                    initialShapeDhm(rX, x, 1, flag);
                }
                if (InitialGuess == initialGuess.dhm2)
                {
                    initialShapeDhm(rX, x, 2, flag);
                }
                if (InitialGuess == initialGuess.dhm3)
                {
                    initialShapeDhm(rX, x, 3, flag);
                }
                if (InitialGuess == initialGuess.dcm)
                {
                    initialShapeDcm(rX, x);
                }
                if (InitialGuess == initialGuess.scpl)
                {
                    initialShapeScpl(rX, x);
                }
                if (InitialGuess == initialGuess.scps)
                {
                    initialShapeScps(rX, x);
                }
            }
        }
        private void Menu_InitialClicked(Object sender, EventArgs e)
        {
            System.Windows.Forms.ToolStripMenuItem _m = sender as System.Windows.Forms.ToolStripMenuItem;
            if (Clock == -1)
            {
                if (_m == __ig1)
                {
                    InitialGuess = initialGuess.prj;
                }
                if (_m == __ig21)
                {
                    InitialGuess = initialGuess.dhm1;
                }
                if (_m == __ig22)
                {
                    InitialGuess = initialGuess.dhm2;
                }
                if (_m == __ig23)
                {
                    InitialGuess = initialGuess.dhm3;
                }
                if (_m == __ig3)
                {
                    InitialGuess = initialGuess.dcm;
                }
                if (_m == __ig5)
                {
                    InitialGuess = initialGuess.scps;
                }
            }

            reset = true;
            ExpireSolution(true);
            //            computeInitialGuess(refPos, pos, true);
        }
        public void initialShapePrj(double[,] rX, double[,] x)
        {
            for (int i = 0; i < x.GetLength(0); i++)
            {
                x[i, 0] = rX[i, 0];
                x[i, 1] = rX[i, 1];
                x[i, 2] = 0;
            }
            _isFixedBoundary = false;
        }
        bool isBoundary(int n)
        {
            for (int i = 0; i < __boundary.Count(); i++)
            {
                if (__boundary[i].index == n)
                {
                    return true;
                }
            }
            return false;
        }
        double[] node = new double[3];
        double[,] node2 = new double[2, 3];
        double[] vals = new double[2];
        double  dt;
        List<Guid> fixedPointGuids = null;
        void deleteFixedPoints()
        {
            if (fixedPointGuids != null)
            {
                foreach (Guid g in fixedPointGuids)
                {
                    Rhino.RhinoDoc.ActiveDoc.Objects.Delete(g, true);
                }
            }
            fixedPointGuids = null;

        }
        void createFixedPoints()
        {
            if (fixedPointGuids != null)
            {
                deleteFixedPoints();
            }
            fixedPointGuids = new List<Guid>();
            for (int i = 0; i < __boundary.Count(); i++)
            {
                boundaryVertex bV = __boundary[i];
                fixedPointGuids.Add(Rhino.RhinoDoc.ActiveDoc.Objects.AddPoint(pos[bV.index, 0] + Shifting, pos[bV.index, 1], 0));
            }
        }
        void projBoundary()
        {
            Rhino.RhinoDoc.ReplaceRhinoObject -= RhinoDoc_ReplaceRhinoObject;
            //find the varycentric
            if (fixedPointGuids == null || fixedPointGuids.Count() != __boundary.Count())
            {
                createFixedPoints();
            }
            double x, y, z;
            for (int i = 0; i < __boundary.Count(); i++)
            {
                x = refPos[__boundary[i].index, 0] + Shifting;
                y = refPos[__boundary[i].index, 1];
                z = 0;
                Rhino.Geometry.Point3d P = new Rhino.Geometry.Point3d(x, y, z);
                Rhino.RhinoDoc.ActiveDoc.Objects.Replace(fixedPointGuids[i], P);
            }
            Rhino.RhinoDoc.ReplaceRhinoObject += RhinoDoc_ReplaceRhinoObject;
        }

        void circleBoundary()
        {
            Rhino.RhinoDoc.ReplaceRhinoObject -= RhinoDoc_ReplaceRhinoObject;
            //find the varycentric
            double cx = 0, cy = 0;
            for (int i = 0; i < nParticles; i++)
            {
                cx += refPos[i, 0];
                cy += refPos[i, 1];
            }
            cx /= nParticles;
            cy /= nParticles;
            cx += Shifting;
            //Calculate radius
            double totalLength = __boundary[__boundary.Count() - 1].param;
            double R = __boundary[__boundary.Count() - 1].param / (2 * Math.PI);
            if (fixedPointGuids == null || fixedPointGuids.Count() != __boundary.Count())
            {
                createFixedPoints();
            }
            for (int i = 0; i < __boundary.Count(); i++)
            {
                double x, y, z;
                double theta = 2 * Math.PI * ((double)i / __boundary.Count());
                x = cx + R * Math.Cos(theta);
                y = cy + R * Math.Sin(theta);
                z = 0;
                Rhino.Geometry.Point3d P = new Rhino.Geometry.Point3d(x, y, z);
                Rhino.RhinoDoc.ActiveDoc.Objects.Replace(fixedPointGuids[i], P);
            }
            Rhino.RhinoDoc.ReplaceRhinoObject += RhinoDoc_ReplaceRhinoObject;
        }

        void RhinoDoc_ReplaceRhinoObject(object sender, Rhino.DocObjects.RhinoReplaceObjectEventArgs e)
        {
            if (Clock == -1)
            {
                if (isInitialized && fixedPointGuids != null)
                {
                    for (int i = 0; i < __boundary.Count(); i++)
                    {
                        Guid gi = fixedPointGuids[i];
                        if (e.ObjectId.CompareTo(gi) == 0)
                        {
                            var PO = e.NewRhinoObject as Rhino.DocObjects.PointObject;
                            var P = PO.PointGeometry;
                            pos[__boundary[i].index, 0] = P.Location.X - Shifting;
                            pos[__boundary[i].index, 1] = P.Location.Y;
                            computeInitialGuess(refPos, pos, false);
                            __update();
                            this.ExpirePreview(true);
                            return;
                        }
                    }
                }
            }
        }
        void rectangleBoundary()
        {
            Rhino.RhinoDoc.ReplaceRhinoObject -= RhinoDoc_ReplaceRhinoObject;
            int _S = shift % (__boundary.Count());

            //find the varycentric
            double cx = 0, cy = 0;
            for (int i = 0; i < nParticles; i++)
            {
                cx += refPos[i, 0];
                cy += refPos[i, 1];
            }
            cx /= nParticles;
            cy /= nParticles;
            cx += Shifting;
            //Calculate radius
            double totalLength = __boundary[__boundary.Count() - 1].param;
            double R = __boundary[__boundary.Count() - 1].param / 4;
            if (fixedPointGuids == null || fixedPointGuids.Count() != __boundary.Count())
            {
                createFixedPoints();
            }
            for (int i = 0; i < (int)__boundary.Count() / 4; i++)
            {
                double x, y, z;
                int range = (int)__boundary.Count() / 4 - 0;
                x = cx - R / 2;
                y = cy - R / 2 + (double)i * R / range;
                z = 0;
                Rhino.Geometry.Point3d P = new Rhino.Geometry.Point3d(x, y, z);
                Rhino.RhinoDoc.ActiveDoc.Objects.Replace(fixedPointGuids[(i + shift) % (__boundary.Count())], P);
            }
            for (int i = (int)__boundary.Count() / 4; i < (int)__boundary.Count() / 2; i++)
            {
                double x, y, z;
                int _i = i - (int)__boundary.Count() / 4;
                int range = ((int)__boundary.Count() / 2) - ((int)__boundary.Count() / 4);
                x = cx - R / 2 + (double)_i * R / range;
                y = cy + R / 2;
                z = 0;
                Rhino.Geometry.Point3d P = new Rhino.Geometry.Point3d(x, y, z);
                Rhino.RhinoDoc.ActiveDoc.Objects.Replace(fixedPointGuids[(i + shift) % (__boundary.Count())], P);
            }
            for (int i = (int)__boundary.Count() / 2; i < (int)__boundary.Count() * 3 / 4; i++)
            {
                double x, y, z;
                int _i = i - (int)__boundary.Count() / 2;
                int range = ((int)__boundary.Count() * 3 / 4) - ((int)__boundary.Count() / 2);
                x = cx + R / 2;
                y = cy + R / 2 - (double)_i * R / range;
                z = 0;
                Rhino.Geometry.Point3d P = new Rhino.Geometry.Point3d(x, y, z);
                Rhino.RhinoDoc.ActiveDoc.Objects.Replace(fixedPointGuids[(i + shift) % (__boundary.Count())], P);
            }
            for (int i = (int)__boundary.Count() * 3 / 4; i < (int)__boundary.Count(); i++)
            {
                int _i = i - (int)__boundary.Count() * 3 / 4;
                double x, y, z;
                int range = (int)__boundary.Count() - ((int)__boundary.Count() * 3 / 4);
                x = cx + R / 2 - (double)_i * R / range;
                y = cy - R / 2;
                z = 0;
                Rhino.Geometry.Point3d P = new Rhino.Geometry.Point3d(x, y, z);
                Rhino.RhinoDoc.ActiveDoc.Objects.Replace(fixedPointGuids[(i + shift) % (__boundary.Count())], P);
            }
            Rhino.RhinoDoc.ReplaceRhinoObject += RhinoDoc_ReplaceRhinoObject;

        }
        void updateFixedPoints()
        {
            Rhino.RhinoDoc.ReplaceRhinoObject -= RhinoDoc_ReplaceRhinoObject;
            if (__boundary.Count() == fixedPointGuids.Count())
            {
                for (int i = 0; i < __boundary.Count(); i++)
                {
                    boundaryVertex bV = __boundary[i];
                    Rhino.DocObjects.PointObject obj = (Rhino.DocObjects.PointObject)Rhino.RhinoDoc.ActiveDoc.Objects.Find(fixedPointGuids[i]);
                    pos[bV.index, 0] = obj.PointGeometry.Location.X - Shifting;
                    pos[bV.index, 1] = obj.PointGeometry.Location.Y;
                    pos[bV.index, 2] = 0;
                    if (obj.PointGeometry.Location.Z != 0)
                    {
                        Rhino.RhinoDoc.ReplaceRhinoObject -= RhinoDoc_ReplaceRhinoObject;
                        Rhino.RhinoDoc.ActiveDoc.Objects.Transform(fixedPointGuids[i], Rhino.Geometry.Transform.Translation(0, 0, -obj.PointGeometry.Location.Z), true);
                        Rhino.RhinoDoc.ReplaceRhinoObject += RhinoDoc_ReplaceRhinoObject;
                    }
                }
            }
            Rhino.RhinoDoc.ReplaceRhinoObject += RhinoDoc_ReplaceRhinoObject;
        }
        System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
        public int P1 = 10, P2 = 20, D1 = 0, D2 = 0;
        public int uDim, vDim;
        public int nU, nV;
        Rhino.Geometry.Mesh initializeNurbs(Rhino.Geometry.NurbsSurface S)
        {
            double[] uKnot;
            double[] vKnot;

            nU = S.Points.CountU;
            nV = S.Points.CountV;
            int N = nU * nV;
            int uDim = S.OrderU;
            int vDim = S.OrderV;
            int uDdim = S.OrderU - 1;
            int vDdim = S.OrderV - 1;
            Rhino.Geometry.Mesh mesh = new Rhino.Geometry.Mesh();
            for (int j = 0; j < nV; j++)
            {
                for (int i = 0; i < nU; i++)
                {
                    Rhino.Geometry.ControlPoint cPoint = S.Points.GetControlPoint(i, j);
                    mesh.Vertices.Add(cPoint.Location);
                }
            }
            for (int j = 0; j < nV - 1; j++)
            {
                for (int i = 0; i < nU - 1; i++)
                {
                    int D1 = j * nU + i;
                    int D2 = j * nU + i + 1;
                    int D3 = (j + 1) * nU + i + 1;
                    int D4 = (j + 1) * nU + i;
                    mesh.Faces.AddFace(new Rhino.Geometry.MeshFace(D1, D2, D3, D4));
                }
            }

            uKnot = new double[nU - uDdim + 1 + uDdim * 2];
            vKnot = new double[nV - vDdim + 1 + vDdim * 2];
            for (int i = 0; i < uDdim; i++)
            {
                uKnot[i] = 0;
            }
            for (int i = 0; i < vDdim; i++)
            {
                vKnot[i] = 0;
            }
            for (int i = 0; i < nU - uDdim + 1; i++)
            {
                uKnot[i + uDdim] = i;
            }
            for (int i = 0; i < nV - vDdim + 1; i++)
            {
                vKnot[i + vDdim] = i;
            }
            for (int i = 0; i < uDdim; i++)
            {
                uKnot[i + nU + 1] = nU - uDdim;
            }
            for (int i = 0; i < vDdim; i++)
            {
                vKnot[i + nV + 1] = nV - vDdim;
            }
            elemList.Clear();
            for (int i = 1; i < nV - vDdim + 1; i++)
            {
                for (int j = 1; j < nU - uDdim + 1; j++)
                {
                    int[] index = new int[uDim * vDim];
                    for (int k = 0; k < vDim; k++)
                    {
                        for (int l = 0; l < uDim; l++)
                        {
                            index[k * uDim + l] = (i - 1 + k) * nU + j - 1 + l;
                        }
                    }
                    if (uDim == 2 & vDim == 2)
                    {
                        elemList.Add(new Minilla3D.Elements.N22D2(index, j, i, uKnot, vKnot));
                    }
                    if (uDim == 2 & vDim == 3)
                    {
                        elemList.Add(new Minilla3D.Elements.N23D2(index, j, i, uKnot, vKnot));
                    }
                    if (uDim == 3 & vDim == 2)
                    {
                        elemList.Add(new Minilla3D.Elements.N32D2(index, j, i, uKnot, vKnot));
                    }
                    if (uDim == 2 & vDim == 4)
                    {
                        elemList.Add(new Minilla3D.Elements.N24D2(index, j, i, uKnot, vKnot));
                    }
                    if (uDim == 4 & vDim == 2)
                    {
                        elemList.Add(new Minilla3D.Elements.N42D2(index, j, i, uKnot, vKnot));
                    }
                    if (uDim == 2 & vDim == 5)
                    {
                        elemList.Add(new Minilla3D.Elements.N25D2(index, j, i, uKnot, vKnot));
                    }
                    if (uDim == 5 & vDim == 2)
                    {
                        elemList.Add(new Minilla3D.Elements.N52D2(index, j, i, uKnot, vKnot));
                    }
                    if (uDim == 3 & vDim == 3)
                    {
                        elemList.Add(new Minilla3D.Elements.N33D2(index, j, i, uKnot, vKnot));
                    }
                    if (uDim == 3 & vDim == 4)
                    {
                        elemList.Add(new Minilla3D.Elements.N34D2(index, j, i, uKnot, vKnot));
                    }
                    if (uDim == 4 & vDim == 3)
                    {
                        elemList.Add(new Minilla3D.Elements.N43D2(index, j, i, uKnot, vKnot));
                    }
                    if (uDim == 3 & vDim == 5)
                    {
                        elemList.Add(new Minilla3D.Elements.N35D2(index, j, i, uKnot, vKnot));
                    }
                    if (uDim == 5 & vDim == 3)
                    {
                        elemList.Add(new Minilla3D.Elements.N53D2(index, j, i, uKnot, vKnot));
                    }
                    if (uDim == 4 & vDim == 4)
                    {
                        elemList.Add(new Minilla3D.Elements.N44D2(index, j, i, uKnot, vKnot));
                    }
                    if (uDim == 4 & vDim == 5)
                    {
                        elemList.Add(new Minilla3D.Elements.N45D2(index, j, i, uKnot, vKnot));
                    }
                    if (uDim == 5 & vDim == 4)
                    {
                        elemList.Add(new Minilla3D.Elements.N54D2(index, j, i, uKnot, vKnot));
                    }
                    if (uDim == 5 & vDim == 5)
                    {
                        elemList.Add(new Minilla3D.Elements.N55D2(index, j, i, uKnot, vKnot));
                    }

                }
            }

            return mesh;
        }
        void computeHessArea(SparseDoubleArray hess)
        {
            //hess.Clear();
            //hess.RemoveZeros();
            for (int i = 0; i < nParticles * 2; i++)
            {
                for (int j = 0; j < nParticles * 2; j++)
                {
                    hess[i, j] = 0;
                }
            }
            foreach (Rhino.Geometry.MeshFace t in triMesh.Faces)
            {
                double g = 0.5;
                if (t.IsTriangle)
                {
                    int A = t.A;
                    int B = t.B;
                    int C = t.C;
                    Rhino.Geometry.Vector3d AB = triMesh.Vertices[A] - triMesh.Vertices[B];
                    Rhino.Geometry.Vector3d BC = triMesh.Vertices[B] - triMesh.Vertices[C];
                    Rhino.Geometry.Vector3d N = Rhino.Geometry.Vector3d.CrossProduct(AB, BC);
                    if (N.Z > 0)
                    {
                        hess[A * 2, B * 2 + 1] += g;
                        hess[A * 2 + 1, B * 2] += -g;
                        hess[A * 2, C * 2 + 1] += -g;
                        hess[A * 2 + 1, C * 2] += g;

                        hess[B * 2, C * 2 + 1] += g;
                        hess[B * 2 + 1, C * 2] += -g;
                        hess[B * 2, A * 2 + 1] += -g;
                        hess[B * 2 + 1, A * 2] += g;

                        hess[C * 2, A * 2 + 1] += g;
                        hess[C * 2 + 1, A * 2] += -g;
                        hess[C * 2, B * 2 + 1] += -g;
                        hess[C * 2 + 1, B * 2] += g;
                    }
                    else
                    {
                        hess[A * 2, B * 2 + 1] += -g;
                        hess[A * 2 + 1, B * 2] += g;
                        hess[A * 2, C * 2 + 1] += g;
                        hess[A * 2 + 1, C * 2] += -g;

                        hess[B * 2, C * 2 + 1] += -g;
                        hess[B * 2 + 1, C * 2] += g;
                        hess[B * 2, A * 2 + 1] += g;
                        hess[B * 2 + 1, A * 2] += -g;

                        hess[C * 2, A * 2 + 1] += -g;
                        hess[C * 2 + 1, A * 2] += g;
                        hess[C * 2, B * 2 + 1] += g;
                        hess[C * 2 + 1, B * 2] += -g;
                    }
                }
            }
        }
        enum state { mesh, nurbs };
        state internalState;
        Rhino.Geometry.Mesh inputMesh = null;
        mikity.GeometryProcessing.MeshStructure meshStructure;
        private void firstActions()
        {
            P1 = nParticles - 1-D2;
            P2 = 0+D1;

            Minilla3D.Materials.harmonicMaterial hm = new Minilla3D.Materials.harmonicMaterial();
            gS.setMaterial(hm.getMaterial());
            gS.computeAll(refPos);
            gS.memoryMetric();
            gS.memoryVolume();
            gS.computeAll(refPos);
            area = gS.getTotalVolume();
            refArea = area;
            gS.setMaterial(_material.getMaterial());
            gS.computeAll(refPos);
            initialShapePrj(refPos, pos);
            gS.setMaterial(_material.getMaterial());
            gS.computeAll(pos);
            cN.computeAll(refPos);
            cN.memoryMetric();
            cN.memoryVolume();
            if (stress || _conformal) gS.computeEigenVectors();
            gS.computeHessEd(refPos);
            gS.getHessian(hessEd, 2);
            computeHessArea(hessA);

            isInitialized = true;
        }
        string dbg = "";
        int ctime=0, iteration=0;
        protected override void SolveInstance(Grasshopper.Kernel.IGH_DataAccess DA)
        {
            Rhino.Geometry.NurbsSurface inputNurbs = null;
            Object inputGeometry = null;
            dbg = "";
            double area = 0;
            if (full != null) dt = full.getDt(); else dt = 0.1;
            if (_go == false)
            {
                isInitialized = false;
                if (reset == false)
                {
                    if (!DA.GetData(0, ref inputGeometry)) { isInitialized = false; return; }
                    if (inputGeometry is Grasshopper.Kernel.Types.GH_Mesh)
                    {
                        inputMesh = (inputGeometry as Grasshopper.Kernel.Types.GH_Mesh).Value;
                        internalState = state.mesh;
                    }
                    else if (inputGeometry is Grasshopper.Kernel.Types.GH_Surface)
                    {
                        inputNurbs = (inputGeometry as Grasshopper.Kernel.Types.GH_Surface).Value.Surfaces[0].ToNurbsSurface();
                        Rhino.Geometry.ControlPoint cp = inputNurbs.Points.GetControlPoint(0, 0);
                        inputMesh = initializeNurbs(inputNurbs);
                        internalState = state.nurbs;
                    }
                    else
                    {
                        isInitialized = false;
                        return;
                    }
                    if (!DA.GetData(1, ref D1)) { D1 = 0; }
                    if (!DA.GetData(2, ref D2)) { D2 = 0; }
                    P1 = nParticles - 1 - D2;
                    P2 = 0 + D1;
                    m = inputMesh.DuplicateMesh();
                    if (internalState == state.nurbs)
                    {
                        ns = inputNurbs.Duplicate() as Rhino.Geometry.NurbsSurface;
                    }
                    initialize(m, m);
                    meshStructure=mikity.GeometryProcessing.MeshStructure.CreateFrom(m);
                    if (internalState == state.mesh)
                    {
                        this.BKGT = GetBKGT(m);
                        this.DVPW = GetDVPW(m);
                    }
                    else if (internalState == state.nurbs)
                    {
                        this.BKGT = GetBKGT(ns);
                        this.DVPW = GetDVPW(ns);
                    }
                    firstActions();
                }
                else
                {
                    reset = false;
                }
                isInitialized = true;
                if (_isFixedBoundary == true)
                {
                    if (fixedPointGuids == null)
                    {
                        createFixedPoints();
                    }
                    if (fixedPointGuids.Count() != __boundary.Count())
                    {
                        createFixedPoints();
                    }
                    updateFixedPoints();
                }
                else
                {
                    deleteFixedPoints();
                }
                computeInitialGuess(refPos, pos, true);


                if (stress || _conformal)
                {
                    gS.computeAll(pos);
                    gS.computeEigenVectors();
                }
                gS.computeVolume(pos);
                area = gS.getTotalVolume();
                __update();

                Clock = -1;
                dbg += "AreaRatio=" + (area / refArea).ToString()+"\n";
                dbg += "Area=" + area.ToString()+"\n";
                dbg += "refArea=" + refArea.ToString()+"\n";
                if (_conformal)
                {
                    dbg += "max(Mc)=" + maxC.ToString() + "\n";
                    dbg += "min(Mc)=" + minC.ToString() + "\n";
                }
                full.setDbgText(dbg);
            }
            else
            {
                if (!DA.GetData(1, ref D1)) { D1 = 0; }
                if (!DA.GetData(2, ref D2)) { D2 = 0; }
                P1 = nParticles - 1 - D2;
                P2 = 0 + D1;
                reset = false;
                sw.Reset();
                
                int S = 0;
                do
                {
                    S++;
                    iteration++;
                    sw.Start();
                    Clock++;
                    __UPDATE();
                    gS.computeVolume(pos);
                    area = gS.getTotalVolume();
                    if (_isFixedBoundary == true)
                    {
                        if (fixedPointGuids == null)
                        {
                            createFixedPoints();
                        }
                        if (fixedPointGuids.Count() != __boundary.Count())
                        {
                            createFixedPoints();
                        }
                        updateFixedPoints();
                    }
                    else
                    {
                        deleteFixedPoints();
                    }
                    if (type == modelType.wire)
                    {
                        foreach (face f in faces)
                        {
                            _wireMaterial(f, refPos, pos);
                        }
                        cN.computeAll(pos);
                        cN.getGrad(force);
                        gS.computeAll(pos);
                    }
                    else if (type == modelType.membrane)
                    {
                        if (mat == materialChoice.SCP)
                        {
                            Minilla3D.Materials.conformalMaterial cm = _material as Minilla3D.Materials.conformalMaterial;
                            cm.area = area;
                        }
                        gS.computeAll(pos);
                        gS.getGrad(force);
                        cN.computeAll(pos);
                        meshStructure.Update(pos);
                        if(_fixFlip)fixFlip(meshStructure);
                    }
                    threeTerm();
                    sw.Stop();
                } while (sw.ElapsedMilliseconds < 25);
                ctime += (int)sw.ElapsedMilliseconds;
                if (stress || _conformal) gS.computeEigenVectors();
                __update();
                //dbg += "repeatCycle=" + S.ToString() + "\n";
                //dbg += "AreaRatio=" + (area / refArea).ToString() + "\n";
                dbg += "itr:" + iteration.ToString() + "  tCt:" + ctime.ToString() + "ms\n";
                dbg += "Area=" + area.ToString() + "\n";
                //dbg += "refArea=" + refArea.ToString() + "\n";
                //if (_conformal)
                {
                    dbg += "max(Mc)=" + maxC.ToString() + "\n";
                    dbg += "min(Mc)=" + minC.ToString() + "\n";
                }
                dbg += "|F|=" + normW.ToString() + "\n";
                full.setDbgText(dbg);
                DA.SetData(0, maxC);
                DA.SetData(1, minC);
                DA.SetData(2, intC);
                
            }
        }
        double maxC = 0;
        double minC = 0;
        double intC = 0;
        double normW = 0;
        double[] V1 = new double[3];
        double[] V2 = new double[3];
        double[] g1 = new double[3];
        double[] g2 = new double[3];
        double[] V3 = new double[3];
        private void fixFlip(mikity.GeometryProcessing.MeshStructure model)
        {
            foreach (var v in model.innerVertices)
            {
                double val = 0;
                bool flag = false;
                foreach (var he in v.star)
                {
                    var ut = he.prev;
                    int P1 = he.next.P.N;
                    int P2 = ut.P.N;
                    int P3 = v.N;
                    V1[0] = he.next.P.x - v.x;
                    V1[1] = he.next.P.y - v.y;
                    V1[2] = 0;
                    V2[0] = ut.P.x - v.x;
                    V2[1] = ut.P.y - v.y;
                    V2[2] = 0;
                    V3[0] = 0;
                    V3[1] = 0;
                    V3[2] = V1[0] * V2[1] - V1[1] * V2[0];
                    if (val == 0) val = V3[2];
                    if (val * V3[2] < 0) { flag = true; break; }
                }
                if (flag)
                {
                    double x = 0;
                    double y = 0;
                    foreach (var he in v.star)
                    {
                        x += he.next.P.x;
                        y += he.next.P.y;
                    }
                    x /= v.star.Count;
                    y /= v.star.Count;
                    int P1 = v.N;
                    //vel[P1 * 3 + 0] = (x - pos[P1, 0]) / dt / 10d;
                    //vel[P1 * 3 + 1] = (y - pos[P1, 1]) / dt / 10d;
                    //pos[P1, 0] = vel[P1 * 3 + 0] * dt;
                    //pos[P1, 1] = vel[P1 * 3 + 1] * dt;
                    pos[P1, 0] = (x - pos[P1, 0]) / 10d + pos[P1, 0];
                    pos[P1, 1] = (y - pos[P1, 1]) / 10d + pos[P1, 1];
                    //v.x = pos[P1, 0];
                    //v.y = pos[P1, 1];
                }
            }
            if (!_isFixedBoundary)
            {
                foreach (var v in model.outerVertices)
                {
                    double val = 0;
                    bool flag = false;
                    foreach (var he in v.star)
                    {
                        var ut = he.prev;
                        int P1 = he.next.P.N;
                        int P2 = ut.P.N;
                        int P3 = v.N;
                        V1[0] = he.next.P.x - v.x;
                        V1[1] = he.next.P.y - v.y;
                        V1[2] = 0;
                        V2[0] = ut.P.x - v.x;
                        V2[1] = ut.P.y - v.y;
                        V2[2] = 0;
                        V3[0] = 0;
                        V3[1] = 0;
                        V3[2] = V1[0] * V2[1] - V1[1] * V2[0];
                        if (val == 0) val = V3[2];
                        if (val * V3[2] < 0) { flag = true; break; }
                    }
                    if (flag)
                    {
                        int P1 = v.N;
                        double x = (v.hf_begin.next.P.x + v.hf_end.P.x) / 2d;
                        double y = (v.hf_begin.next.P.y + v.hf_end.P.y) / 2d;
                        //vel[P1 * 3 + 0] = (x - pos[P1, 0]) / dt / 1d;
                        //vel[P1 * 3 + 1] = (y - pos[P1, 1]) / dt / 1d;
                        pos[P1, 0] = x;// (x - pos[P1, 0]) + pos[P1, 0];
                        pos[P1, 1] = y;// (y - pos[P1, 1]) + pos[P1, 1];
                        //pos[P1, 0] += vel[P1 * 3 + 0] * dt;
                        //pos[P1, 1] += vel[P1 * 3 + 1] * dt;
                        //v.x = pos[P1, 0];
                        //v.y = pos[P1, 1];
                    }
                }
            }
        }
        private void __update()
        {
            double S2 = 20;
            if (internalState == state.mesh)
            {
                for (int i = 0; i < nParticles; i++)
                {
                    m.Vertices[i] = new Rhino.Geometry.Point3f((float)pos[i, 0] + Shifting, (float)pos[i, 1], (float)pos[i, 2]);
                }
            }
            for (int i = 0; i < nParticles; i++)
            {
                triMesh.Vertices[i] = new Rhino.Geometry.Point3f((float)pos[i, 0] + Shifting, (float)pos[i, 1], (float)pos[i, 2]);
            }
            if (internalState == state.nurbs)
            {
                for (int i = 0; i < nV; i++)
                {
                    for (int j = 0; j < nU; j++)
                    {
                        int k = j + i * nU;
                        ns.Points.SetControlPoint(j, i, new Rhino.Geometry.ControlPoint(pos[k, 0] + Shifting, pos[k, 1], pos[k, 2]));
                    }
                }
            }
            iP.Clear();
            bV.Clear();
            eVT.Clear();
            eVC.Clear();
            reF.Clear();
            cfm.Clear();
            if (_isFixedBoundary || mat == materialChoice.DCM)
            {
                for (int i = 0; i < nParticles; i++)
                {
                    reF.Add(new Rhino.Geometry.Line(pos[i, 0] + Shifting, pos[i, 1], 0, pos[i, 0] + Re[i, 0] * S2 + Shifting, pos[i, 1] + Re[i, 1] * S2, 0));
                }
            }
            if (_base)
            {
                foreach (Minilla3D.Elements.managedElement e in elemList)
                {
                    for (int i = 0; i < e.nIntPoint; i++)
                    {
                        e.getGlobalCoord(node, i);
                        e.getBaseVectors(node2, i);
                        double S = 0.15;
                        bV.Add(new Rhino.Geometry.Line(node[0] + Shifting, node[1], node[2], node[0] + node2[0, 0] * S + Shifting, node[1] + node2[0, 1] * S, node[2] + node2[0, 2] * S));
                        bV.Add(new Rhino.Geometry.Line(node[0] + Shifting, node[1], node[2], node[0] + node2[1, 0] * S + Shifting, node[1] + node2[1, 1] * S, node[2] + node2[1, 2] * S));
                    }
                }

            }
            if (_intPoint)
            {
                if (type == modelType.membrane)
                {
                    foreach (Minilla3D.Elements.managedElement e in elemList)
                    {
                        for (int i = 0; i < e.nIntPoint; i++)
                        {
                            e.getGlobalCoord(node, i);
                            iP.Add(new Rhino.Geometry.Point3d(node[0] + Shifting, node[1], node[2]));
                        }
                    }
                }
            }
            if (stress)
            {
                foreach (Minilla3D.Elements.managedElement e in elemList)
                {
                    for (int i = 0; i < e.nIntPoint; i++)
                    {
                        e.getGlobalCoord(node, i);
                        e.getEigenVectors(node2, vals, i);
                        vals[0] = Math.Log(vals[0]);
                        vals[1] = Math.Log(vals[1]);
                        vals[0] *= 5;
                        vals[1] *= 5;
                        if (vals[0] > 0)
                        {
                            eVT.Add(new Rhino.Geometry.Line(node[0] + Shifting, node[1], node[2], node[0] + node2[0, 0] * vals[0] + Shifting, node[1] + node2[0, 1] * vals[0], node[2] + node2[0, 2] * vals[0]));
                            eVT.Add(new Rhino.Geometry.Line(node[0] + Shifting, node[1], node[2], node[0] - node2[0, 0] * vals[0] + Shifting, node[1] - node2[0, 1] * vals[0], node[2] - node2[0, 2] * vals[0]));
                        }
                        else
                        {
                            eVC.Add(new Rhino.Geometry.Line(node[0] + node2[0, 0] * vals[0] + Shifting, node[1] + node2[0, 1] * vals[0], node[2] + node2[0, 2] * vals[0], node[0] + Shifting, node[1], node[2]));
                            eVC.Add(new Rhino.Geometry.Line(node[0] - node2[0, 0] * vals[0] + Shifting, node[1] - node2[0, 1] * vals[0], node[2] - node2[0, 2] * vals[0], node[0] + Shifting, node[1], node[2]));
                        }
                        if (vals[1] > 0)
                        {
                            eVT.Add(new Rhino.Geometry.Line(node[0] + Shifting, node[1], node[2], node[0] + node2[1, 0] * vals[1] + Shifting, node[1] + node2[1, 1] * vals[1], node[2] + node2[1, 2] * vals[1]));
                            eVT.Add(new Rhino.Geometry.Line(node[0] + Shifting, node[1], node[2], node[0] - node2[1, 0] * vals[1] + Shifting, node[1] - node2[1, 1] * vals[1], node[2] - node2[1, 2] * vals[1]));
                        }
                        else
                        {
                            eVC.Add(new Rhino.Geometry.Line(node[0] + node2[1, 0] * vals[1] + Shifting, node[1] + node2[1, 1] * vals[1], node[2] + node2[1, 2] * vals[1], node[0] + Shifting, node[1], node[2]));
                            eVC.Add(new Rhino.Geometry.Line(node[0] - node2[1, 0] * vals[1] + Shifting, node[1] - node2[1, 1] * vals[1], node[2] - node2[1, 2] * vals[1], node[0] + Shifting, node[1], node[2]));
                        }
                    }
                }
            }
            //if (_conformal)
            {
                maxC = 0;
                minC = Double.MaxValue;
                intC = 0;
                double weight = 0;
                foreach (Minilla3D.Elements.managedElement e in elemList)
                {
                    for (int i = 0; i < e.nIntPoint; i++)
                    {
                        e.getGlobalCoord(node, i);
                        weight = e.getEigenVectors(node2, vals, i);
                        double criteria = 0;
                        //if (_conformal)
                        {
                            criteria = ((vals[1] - vals[0]) * (vals[1] - vals[0])) / (vals[1] * vals[0])    ;
                        }
                        if (criteria > maxC) maxC = criteria;
                        if (criteria < minC) minC = criteria;
                        intC += weight * criteria;
                        criteria = criteria * 20;
                        cfm.Add(new Rhino.Geometry.Circle(new Rhino.Geometry.Point3d(node[0] + Shifting, node[1], node[2]), criteria));
                    }
                }
            }
        }

        private void threeTerm()
        {
            
            for (int i = 0; i < nParticles; i++)
            {
                Re[i, 0] = 0;
                Re[i, 1] = 0;
                Re[i, 2] = 0;
            }
            if (_isFixedBoundary == true)
            {
                foreach (boundaryVertex b in __boundary)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        Re[b.index, j] = force[b.index * 3 + j];
                        force[b.index * 3 + j] = 0;
                    }
                }
            }
            if (mat == materialChoice.DCM)
            {
                Re[P1, 0] = force[P1 * 3 + 0];
                Re[P1, 1] = force[P1 * 3 + 1];
                Re[P2, 0] = force[P2 * 3 + 0];
                Re[P2, 1] = force[P2 * 3 + 1];
                force[P1 * 3 + 0] = 0;
                force[P1 * 3 + 1] = 0;
                force[P2 * 3 + 0] = 0;
                force[P2 * 3 + 1] = 0;
            }
            
            acc = force;
            for (int i = 0; i < nParticles; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    if (Double.IsNaN(acc[i*3+ j])) acc[i*3+ j] = 0;
                }
            }
            
            double norm1 = (vel * vel.T)[0,0];
            double norm2 = (vel * acc.T)[0,0];
            double norm3 = (acc * acc.T)[0,0];
            double norm = Math.Sqrt((acc * force.T)[0, 0]);
            normW = norm;
            double f = 0;
            if (norm1 * norm3 != 0)
            {
                f = -norm2 / Math.Sqrt(norm1 * norm3);
            }
            else
            {
                f = 1;
            }
            double damping1 = 0;
            damping1 = Drift0(f);
            if (norm < 1.0) norm = 1.0;
            vel = vel * damping1 - acc * (dt/norm);
            
            for (int i = 0; i < nParticles; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    pos[i, j] = pos[i, j] + vel[i * 3 + j] * dt;
                }
            }

        }
        public BakeGeometry GetBKGT(Rhino.Geometry.Mesh _m)
        {
            return new BakeGeometry((d, a, o) =>
            {
                Rhino.DocObjects.ObjectAttributes a2 = a.Duplicate();
                a2.LayerIndex = 2;
                if (this.type == modelType.membrane)
                {
                    o.Add(d.Objects.AddMesh(_m, a2));
                }
                else if (this.type == modelType.wire)
                {
                    o.Add(d.Objects.AddMesh(triMesh, a2));
                }
            });
        }
        public BakeGeometry GetBKGT(Rhino.Geometry.NurbsSurface _m)
        {
            return new BakeGeometry((d, a, o) =>
            {
                Rhino.DocObjects.ObjectAttributes a2 = a.Duplicate();
                a2.LayerIndex = 2;
                if (this.type == modelType.membrane)
                {
                    o.Add(d.Objects.AddSurface(_m, a2));
                }
                else if (this.type == modelType.wire)
                {
                    o.Add(d.Objects.AddMesh(triMesh, a2));
                }
            });
        }

        public DrawViewPortWire GetDVPW(Rhino.Geometry.Mesh _m)
        {
            return new DrawViewPortWire((args) =>
            {
                if (Hidden)
                {
                    return;
                }
                if (type != modelType.wire)
                {
                    args.Display.DrawMeshWires(_m, System.Drawing.Color.Black);
                }

                if (type != modelType.membrane)
                {
                    args.Display.DrawMeshWires(triMesh, System.Drawing.Color.Brown, 3);
                    args.Display.DrawPoints(triMesh.Vertices.ToPoint3dArray(), Rhino.Display.PointStyle.Simple, 2, System.Drawing.Color.White);
                }
                args.Display.DrawPoints(iP, Rhino.Display.PointStyle.ControlPoint, 1, System.Drawing.Color.Black);
                args.Display.DrawLines(bV, System.Drawing.Color.Red, 1);
                args.Display.DrawLines(eVT, System.Drawing.Color.Cyan);
                args.Display.DrawLines(eVC, System.Drawing.Color.Magenta);
                if ((mat == materialChoice.DCM && Clock >= 0) || (InitialGuess == initialGuess.dcm && Clock == -1))
                {
                    args.Display.DrawPoint(new Rhino.Geometry.Point3d(pos[P1, 0] + Shifting, pos[P1, 1], 0), Rhino.Display.PointStyle.X, 3, System.Drawing.Color.Red);
                    args.Display.DrawPoint(new Rhino.Geometry.Point3d(pos[P2, 0] + Shifting, pos[P2, 1], 0), Rhino.Display.PointStyle.X, 3, System.Drawing.Color.Red);
                }
                if ((mat == materialChoice.DCM || _isFixedBoundary)&&_RF)
                {
                    args.Display.DrawLines(reF, System.Drawing.Color.Red);
                    foreach (Rhino.Geometry.Line l in reF)
                    {
                        args.Display.DrawArrowHead(l.To, l.Direction, System.Drawing.Color.Red, 0.0, l.Length / 5d);
                    }
                }
                foreach (Rhino.Geometry.Line l in eVT)
                {
                    args.Display.DrawArrowHead(l.To, l.Direction, System.Drawing.Color.Cyan, 0.0, l.Length / 5d);
                }
                foreach (Rhino.Geometry.Line l in eVC)
                {
                    args.Display.DrawArrowHead(l.To, l.Direction, System.Drawing.Color.Magenta, 0.0, l.Length / 5d);
                }
                if (_conformal)
                {
                    foreach (Rhino.Geometry.Circle c in cfm)
                    {
                        args.Display.DrawCircle(c, System.Drawing.Color.Purple);
                    }
                }
                if (_isFixedBoundary == true)
                {
                    foreach (boundaryVertex b in __boundary)
                    {
                        args.Display.DrawPoint(_m.Vertices[b.index], Rhino.Display.PointStyle.X, 3, System.Drawing.Color.Red);
                    }
                }
                args.Display.DrawPoint(_m.Vertices[0], Rhino.Display.PointStyle.ActivePoint, 5, System.Drawing.Color.OrangeRed);
                args.Display.DrawPoint(inputMesh.Vertices[0], Rhino.Display.PointStyle.ActivePoint, 5, System.Drawing.Color.OrangeRed);
            });
        }
        public DrawViewPortWire GetDVPW(Rhino.Geometry.NurbsSurface _m)
        {
            return new DrawViewPortWire((args) =>
            {
                if (Hidden)
                {
                    return;
                }
                if (type != modelType.wire)
                {
                    args.Display.DrawSurface(_m, System.Drawing.Color.Red, 1);
                }

                if (type != modelType.membrane)
                {
                    args.Display.DrawMeshWires(triMesh, System.Drawing.Color.Brown, 3);
                    args.Display.DrawPoints(triMesh.Vertices.ToPoint3dArray(), Rhino.Display.PointStyle.Simple, 2, System.Drawing.Color.White);
                }
                args.Display.DrawPoints(iP, Rhino.Display.PointStyle.ControlPoint, 1, System.Drawing.Color.Black);
                args.Display.DrawLines(bV, System.Drawing.Color.Pink, 2);
                args.Display.DrawLines(eVT, System.Drawing.Color.Cyan);
                args.Display.DrawLines(eVC, System.Drawing.Color.Magenta);
                if ((mat == materialChoice.DCM&&Clock>=0)||(InitialGuess==initialGuess.dcm&&Clock==-1))
                {
                    args.Display.DrawPoint(new Rhino.Geometry.Point3d(pos[P1, 0] + Shifting, pos[P1, 1], 0), Rhino.Display.PointStyle.X, 3, System.Drawing.Color.Red);
                    args.Display.DrawPoint(new Rhino.Geometry.Point3d(pos[P2, 0] + Shifting, pos[P2, 1], 0), Rhino.Display.PointStyle.X, 3, System.Drawing.Color.Red);
                }
                if ((mat == materialChoice.DCM || _isFixedBoundary)&&_RF)
                {
                    args.Display.DrawLines(reF, System.Drawing.Color.Red);
                    foreach (Rhino.Geometry.Line l in reF)
                    {
                        args.Display.DrawArrowHead(l.To, l.Direction, System.Drawing.Color.Red, 0.0, l.Length / 5d);
                    }
                }
                foreach (Rhino.Geometry.Line l in eVT)
                {
                    args.Display.DrawArrowHead(l.To, l.Direction, System.Drawing.Color.Cyan, 0.0, l.Length / 5d);
                }
                foreach (Rhino.Geometry.Line l in eVC)
                {
                    args.Display.DrawArrowHead(l.To, l.Direction, System.Drawing.Color.Magenta, 0.0, l.Length / 5d);
                }
                if (_conformal)
                {
                    foreach (Rhino.Geometry.Circle c in cfm)
                    {
                        args.Display.DrawCircle(c, System.Drawing.Color.Purple);
                    }
                }

                if (_isFixedBoundary)
                {
                    foreach (boundaryVertex b in __boundary)
                    {
                        args.Display.DrawPoint(triMesh.Vertices[b.index], Rhino.Display.PointStyle.X, 3, System.Drawing.Color.Red);
                    }
                }
            });
        }
        public override Guid ComponentGuid
        {
            get { return new Guid("b3c112b7-59a1-40a5-9fd4-c63b8356f02c"); }
        }

    }*/
}
