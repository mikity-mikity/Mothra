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
        double[,] Force=null;
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
                
                //model.AddConstrs(bConds, bSenses, bRhs, null);
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
                for (int j = 0; j < nU * nV; j++)
                {
                    x[j, 2] = vars[j].Get(GRB.DoubleAttr.X);
                }
                //model.Set(GRB.DoubleAttr.LB, bVars, lb);
                //model.Set(GRB.DoubleAttr.UB, bVars, ub);
                //model.Optimize();
                model.Dispose();
                env.Dispose();

            }
            catch (GRBException e)
            {
                AddRuntimeMessage(Grasshopper.Kernel.GH_RuntimeMessageLevel.Error, "GUROBI ERROR!" + e.ErrorCode + ". " + e.Message);
                return;
            }
            x2Nurbs(x, airyNurbs);
            createAiryPoints(airyNurbs);
            myMasonry.setupNodesFromList(x);
            myMasonry.computeAiryFunction();
            myMasonry.getResidual(residual);
            System.Windows.Forms.MessageBox.Show(residual.Norm().ToString());
            Force = new double[nU * nV, 3];
            //myMasonry.computeEdgeForce(Force);

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
                force[i * 3 + 0, 0] = 0; // Force[i, 0] / 2.5;
                force[i * 3 + 1, 0] = 0; // Force[i, 1] / 2.5;
                force[i * 3 + 2, 0] = 0; // Force[i, 2] / 2.5;
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
            var F=hess* xx;
            for (int i = 0; i < nParticles; i++)
            {
                double Fx = F[i * 3 + 0, 0];
                double Fy = F[i * 3 + 1, 0];
                double Fz = F[i * 3 + 2, 0];

                if ((Fx * Fx + Fy * Fy + Fz * Fz) > 0.01)
                {
                    F[i * 3 + 0, 0] = Fx;
                    F[i * 3 + 1, 0] = Fy;
                    F[i * 3 + 2, 0] = Fz;
                    Force[i, 0] = Fx;
                    Force[i, 1] = Fy;
                    Force[i, 2] = Fz;
                }
                else
                {
                    F[i * 3 + 0, 0] = 0;
                    F[i * 3 + 1, 0] = 0;
                    F[i * 3 + 2, 0] = 0;
                    Force[i, 0] = 0;
                    Force[i, 1] = 0;
                    Force[i, 2] = 0;
                }
            }
            
            
            origX = DoubleArray.Zeros(nParticles * 3, 1);
            for (int i = 0; i < nParticles; i++)
            {
                origX[i * 3 + 0, 0] = x[i, 0];
                origX[i * 3 + 1, 0] = x[i, 1];
                origX[i * 3 + 2, 0] = x[i, 2];
            }
            shift = new List<int>();
            for (int i = 0; i < nParticles; i++)
            {
                shift.Add(i);
            }
            T1 = (nParticles - 4) * 3 - 1;
            T2 = nParticles * 3 - 1;
            C1 = 0;
            C2 = nParticles - 4;
            for (int i = 0; i < nParticles; i++)
            {
                if (i==0 || i==nU-1|| i==nU*nV-nU || i==nU*nV-1)
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
            shiftArray = new SparseDoubleArray(nParticles * 3, nParticles * 3);
            for (int i = 0; i < nParticles; i++)
            {
                shiftArray[i * 3, shift[i] * 3] = 1;
                shiftArray[i * 3 + 1, shift[i] * 3 + 1] = 1;
                shiftArray[i * 3 + 2, shift[i] * 3 + 2] = 1;
            }


            ED = shiftArray.T.Multiply(hess) as SparseDoubleArray;
            ED = ED.Multiply(shiftArray) as SparseDoubleArray;
            slice1 = new SparseDoubleArray(T1 + 1, T2 + 1);
            slice2 = new SparseDoubleArray(T2 + 1, T2 - T1);
            for (int i = 0; i < T1 + 1; i++)
            {
                slice1[i, i] = 1;
            }
            for (int i = 0; i < T2 - T1; i++)
            {
                slice2[i + T1 + 1, i] = 1;
            }
            DIB = (slice1.Multiply(ED) as SparseDoubleArray).Multiply(slice2) as SparseDoubleArray;
            DII = (slice1.Multiply(ED) as SparseDoubleArray).Multiply(slice1.T) as SparseDoubleArray;
            solver = new SparseLU(DII);
            origX = shiftArray.T * origX;
            fixX = origX.GetSlice(T1 + 1, T2, 0, 0);
            B = -DIB * fixX;
            
            force = (shiftArray.T * F).GetSlice(0, T1, 0, 0);
            B = B + force;
            dx = solver.Solve(B);

            ret = DoubleArray.Zeros(nParticles * 3, 1);
            for (int i = 0; i < T1 + 1; i++)
            {
                ret[i, 0] = dx[i, 0];
            }
            for (int i = T1 + 1; i <= T2; i++)
            {
                ret[i, 0] = fixX[i - T1 - 1, 0];
            }
            xx = shiftArray * ret;
            

            
            
            for (int i = 0; i < nParticles; i++)
            {   
                x[i, 0] = xx[i * 3, 0];
                x[i, 1] = xx[i * 3 + 1, 0];
                x[i, 2] = xx[i * 3 + 2, 0];
            }
            x2Nurbs(x, outputNurbs);
        }
        public override void BakeGeometry(Rhino.RhinoDoc doc, Rhino.DocObjects.ObjectAttributes att, List<Guid> obj_ids)
        {
            Rhino.DocObjects.ObjectAttributes a2 = att.Duplicate();
            a2.LayerIndex = 1;
            Guid id = doc.Objects.AddSurface(airyNurbs, a2);
            obj_ids.Add(id);
            Rhino.DocObjects.ObjectAttributes a3 = att.Duplicate();
            a2.LayerIndex = 2;
            Guid id2 = doc.Objects.AddSurface(outputNurbs, a2);
            obj_ids.Add(id2);
            base.BakeGeometry(doc, att, obj_ids);
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
            double S = 500;
            Nurbs2x(xyNurbs, x);
            myMasonry.setupNodesFromList(x);
            myMasonry.computeGlobalCoord();

            foreach (var e in myMasonry.elemList)
            {
                for (int i = 0; i < e.nIntPoint; i++)
                {
                    e.getEigenVectors(vec, val, i);
                    node = e.getIntPoint(i);
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
                for (int i = 0; i < e.nBIntPoint; i++)
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
            foreach (var e in myMasonry.edgeList)
            {
                for (int i = 0; i < e.nIntPoint; i++)
                {
                    var d = e.getIntPoint(i);
                    var d2 = new Rhino.Geometry.Point3d(d[0], d[1], d[2]);
                    outputP.Add(d2);
                    args.Display.DrawPoint(d2, Rhino.Display.PointStyle.X, 5, System.Drawing.Color.Purple);
                }
            }
            args.Display.DrawPoint(outputNurbs.Points.GetControlPoint(0, 0).Location, Rhino.Display.PointStyle.X, 2, System.Drawing.Color.Red);
            args.Display.DrawPoint(outputNurbs.Points.GetControlPoint(nU-1, 0).Location, Rhino.Display.PointStyle.X, 2, System.Drawing.Color.Red);
            args.Display.DrawPoint(outputNurbs.Points.GetControlPoint(0, nV-1).Location, Rhino.Display.PointStyle.X, 2, System.Drawing.Color.Red);
            args.Display.DrawPoint(outputNurbs.Points.GetControlPoint(nU-1, nV-1).Location, Rhino.Display.PointStyle.X, 2, System.Drawing.Color.Red);
            double S3 = 50;
            if (Force != null)
            {
                for (int i = 0; i < nU * nV; i++)
                {
                    double Fx = Force[i, 0] * S3;
                    double Fy = Force[i, 1] * S3;
                    double Fz = Force[i, 2] * S3;
                    if ((Fx * Fx + Fy * Fy + Fz * Fz) !=0)
                        args.Display.DrawArrow(new Line(x[i, 0], x[i, 1], x[i, 2], x[i, 0] + Fx, x[i, 1] + Fy, x[i, 2] + Fz), System.Drawing.Color.Orange);
                }
            }
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
            Rhino.Geometry.Polyline[] boundary = alternativeMesh.GetNakedEdges();
            
            foreach (var p in boundary)
            {
                for (int i = 0; i < p.Count - 1; i++)
                {
                    var f = alternativeMesh.Vertices.Select((e,index)=>new{e,index}).Where(e=>e.e==p[i]).Select(e=>e.index).First();
                    boundaryIndex.Add(f);
                }
            }
            /*boundaryIndex.Add(0);
            boundaryIndex.Add(nU-1);
            boundaryIndex.Add(nU * nV - 1-nU+1);
            boundaryIndex.Add(nU * nV - 1);
            */
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
                        int[] index2 = new int[uDim];
                        for (int l = 0; l < uDim; l++)
                        {
                            index2[l] = (j - 1) * nU + i - 1 + l;
                        }
                        myMasonry.edgeList.Add(new Minilla3D.Elements.nurbsCurve(uDim, index2, i, uKnot));
                    }
                    if (i == 1)
                    {
                        border = border | Minilla3D.Elements.nurbsElement.border.Left;
                        int[] index2 = new int[vDim];
                        for (int k = 0; k < vDim; k++)
                        {
                            index2[k] = (j - 1 + k) * nU + i - 1;
                        }
                        myMasonry.edgeList.Add(new Minilla3D.Elements.nurbsCurve(vDim, index2, j, vKnot));
                    }
                    if (j == nV - vDdim)
                    {
                        border = border | Minilla3D.Elements.nurbsElement.border.Bottom;
                        int[] index2 = new int[uDim];
                        for (int l = 0; l < uDim; l++)
                        {
                            index2[l] = (j - 1 + (vDim-1)) * nU + i - 1 + l;
                        }
                        myMasonry.edgeList.Add(new Minilla3D.Elements.nurbsCurve(uDim, index2, i, uKnot));
                    }
                    if (i == nU - uDdim)
                    {
                        border = border | Minilla3D.Elements.nurbsElement.border.Right;
                        int[] index2 = new int[vDim];
                        for (int k = 0; k < vDim; k++)
                        {
                            index2[k] = (j - 1 + k) * nU + i - 1 + uDim-1;
                            
                        }
                        myMasonry.edgeList.Add(new Minilla3D.Elements.nurbsCurve(vDim, index2, j, vKnot));
                    }
                    myMasonry.elemList.Add(new Minilla3D.Elements.nurbsElement(uDim, vDim, index, i, j, uKnot, vKnot, border));
                }
            }
        }
    }

}
