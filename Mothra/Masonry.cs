using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ShoNS.Array;
using Rhino.Geometry;
using System.IO;
using System.Reflection;
using Gurobi;
using Minilla3D.Elements;
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
            pManager.AddBooleanParameter("Boundary", "b", "Boundary", Grasshopper.Kernel.GH_ParamAccess.item);
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
        public Rhino.Geometry.Plane[] cuttingPlane=new Plane[4];
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
            myMasonry.computeGlobalCoord();
            int nConstraints2 = myMasonry.totalNumberOfIconst(boundary);
            int nConstraints = myMasonry.totalNumberOfBconst();
            ShoNS.Array.SparseDoubleArray jacobian = new SparseDoubleArray(nConstraints, nU * nV);
            ShoNS.Array.DoubleArray residual = new DoubleArray(nConstraints, 1);
            
            ShoNS.Array.SparseDoubleArray jacobianH = new SparseDoubleArray(nConstraints2, nU * nV);
            myMasonry.GetJacobianOfCurvature(jacobianH,boundary);
            
            myMasonry.getJacobian(jacobian);
            
            try
            {
                ////////////////////////////GUROBI///////////////////
                GRBEnv env = new GRBEnv("qcp.log");
                env.Set(GRB.IntParam.DualReductions, 0);
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
                        lb[count] = -GRB.INFINITY;
                        ub[count] = GRB.INFINITY;
                        type[count] = GRB.CONTINUOUS;
                        name[count] = "N[" + i.ToString("g") + "," + j.ToString("g") + "]";
                    }
                }
                lb[0] = 500;
                ub[0] = 500;
                lb[nU - 1] = 500;
                ub[nU - 1] = 500;
                lb[nU * nV - nU] = 500;
                ub[nU * nV - nU] = 500;
                lb[nU * nV - 1] = 500;
                ub[nU * nV - 1] = 500;

                //lb[(nU / 2) + (nV / 2) * nU-1] = 0;
                //ub[(nU / 2) + (nV / 2) * nU-1] = 0;
                lb[28] = 0;
                ub[28] = 0;
                lb[29] = 0;
                ub[29] = 0;
                lb[36] = 0;
                ub[36] = 0;
                lb[37] = 0;
                ub[37] = 0;
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
                    lb[j + 1] = 0;
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

                GRBQuadExpr obj = new GRBQuadExpr();

                model.Update();
                //Quadratic constraints  detH>0
                for (int i = 0; i < nConstraints2; i+=3)
                {
                    GRBQuadExpr Qlhs = new GRBQuadExpr(phis[i +2] * phis[i+ 2]);
                    GRBQuadExpr Qrhs = new GRBQuadExpr(phis[i +0] * phis[i + 1]);
                    model.AddQConstr(Qlhs, GRB.LESS_EQUAL, Qrhs, null);
                }
                
                if (!boundary)
                {
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
                    for (int i = 0; i < nConstraints; i++)
                    {
                        obj.AddTerm(1, bVars[i], bVars[i]);
                    }
                    model.SetObjective(obj);
                    model.Optimize();
                    switch (model.Get(GRB.IntAttr.Status))
                    {
                        case GRB.Status.OPTIMAL:
                            //System.Windows.Forms.MessageBox.Show("Congraturation");
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
                }
                else
                {
                    //Plane variables
                    GRBVar[][] planeVars = new GRBVar[4][];
                    for (int n = 0; n < 4; n++)
                    {
                        //z(x,y)=-D-Ax-By
                        type = new char[3];
                        name = new string[3];
                        lb = new double[3];
                        ub=new double[3];
                        for (int i = 0; i < 3; i++)
                        {
                            type[i] = GRB.CONTINUOUS;
                            lb[i] = -GRB.INFINITY;
                            ub[i] = GRB.INFINITY;
                        }
                        name[0] = "A" + n.ToString("g");
                        name[1] = "B" + n.ToString("g");
                        name[2] = "D" + n.ToString("g");
                        planeVars[n] = model.AddVars(lb,ub, null, type, name);
                    }
                    model.Update();
                    //GRBQuadExpr obj = new GRBQuadExpr();

                    //Coordinate variables of border integratingpoints
                    var borderVars = new List<GRBVar>[4];
                    var slackVars = new List<GRBVar>[4];
                    borderVars[0] = new List<GRBVar>();
                    borderVars[1] = new List<GRBVar>();
                    borderVars[2] = new List<GRBVar>();
                    borderVars[3] = new List<GRBVar>();
                    slackVars[0] = new List<GRBVar>();
                    slackVars[1] = new List<GRBVar>();
                    slackVars[2] = new List<GRBVar>();
                    slackVars[3] = new List<GRBVar>();
                    var tE = myMasonry.topEdge;
                    var lE = myMasonry.leftEdge;
                    var rE = myMasonry.rightEdge;
                    var bE = myMasonry.bottomEdge;
                    nConstraints=myMasonry.totalNumbrOfTopEdgeIntPoint();
                    jacobian = new SparseDoubleArray(nConstraints, nU * nV);
                    myMasonry.GetJacobianOfTopEdge(jacobian);
                    count=0;
                    for (int i = 0; i < tE.Count; i++)
                    {
                        var e = tE[i];
                        for (int j = 0; j < e.nIntPoint; j++,count++)
                        {
                            borderVars[0].Add(model.AddVar(-GRB.INFINITY, GRB.INFINITY, 0, GRB.CONTINUOUS, "topZ" + i.ToString() + "," + j.ToString()));
                            slackVars[0].Add(model.AddVar(-GRB.INFINITY,GRB.INFINITY,0,GRB.CONTINUOUS,"slack"));
                            model.Update();
                            var con = new GRBLinExpr();
                            var node=e.getIntPoint(j);
                            con.AddTerm(1, borderVars[0][count]);   //z
                            con.AddTerm(1, planeVars[0][2]);        //D
                            con.AddTerm(node[0], planeVars[0][0]);  //Ax
                            con.AddTerm(node[1], planeVars[0][1]);  //By
                            model.AddConstr(con, GRB.EQUAL, 0, null);
                            
                            var con2 = new GRBLinExpr();
                            con2.AddTerm(-1, slackVars[0][count]);
                            for (int k = 0; k < nU * nV; k++)
                            {
                                if (jacobian[count, k] != 0)
                                {
                                    con2.AddTerm(jacobian[count, k], vars[k]);
                                }
                            }
                            model.AddConstr(con2, GRB.EQUAL, 0, null);
                            obj.AddTerm(1, slackVars[0][count], slackVars[0][count]);
                            obj.AddTerm(-2, slackVars[0][count], borderVars[0][count]);
                            obj.AddTerm(1, borderVars[0][count], borderVars[0][count]);
                        }
                    }
                    
                    nConstraints = myMasonry.totalNumbrOfBottomEdgeIntPoint();
                    jacobian = new SparseDoubleArray(nConstraints, nU * nV);
                    myMasonry.GetJacobianOfBottomEdge(jacobian);
                    count = 0;
                    for (int i = 0; i < bE.Count; i++)
                    {
                        var e = bE[i];
                        for (int j = 0; j < e.nIntPoint; j++, count++)
                        {
                            borderVars[1].Add(model.AddVar(-GRB.INFINITY, GRB.INFINITY, 0, GRB.CONTINUOUS, "topZ" + i.ToString() + "," + j.ToString()));
                            slackVars[1].Add(model.AddVar(-GRB.INFINITY,GRB.INFINITY,0,GRB.CONTINUOUS,"slack"));
                            model.Update();
                            var con = new GRBLinExpr();
                            var node=e.getIntPoint(j);
                            con.AddTerm(1, borderVars[1][count]);   //z
                            con.AddTerm(1, planeVars[1][2]);        //D
                            con.AddTerm(node[0], planeVars[1][0]);  //Ax
                            con.AddTerm(node[1], planeVars[1][1]);  //By
                            model.AddConstr(con, GRB.EQUAL, 0, null);
                            
                            var con2 = new GRBLinExpr();
                            con2.AddTerm(-1, slackVars[1][count]);
                            for (int k = 0; k < nU * nV; k++)
                            {
                                if (jacobian[count, k] != 0)
                                {
                                    con2.AddTerm(jacobian[count, k], vars[k]);
                                }
                            }
                            model.AddConstr(con2, GRB.EQUAL, 0, null);
                            obj.AddTerm(1, slackVars[1][count], slackVars[1][count]);
                            obj.AddTerm(-2, slackVars[1][count], borderVars[1][count]);
                            obj.AddTerm(1, borderVars[1][count], borderVars[1][count]);
                        }
                    }
                    
                    nConstraints = myMasonry.totalNumbrOfLeftEdgeIntPoint();
                    jacobian = new SparseDoubleArray(nConstraints, nU * nV);
                    myMasonry.GetJacobianOfLeftEdge(jacobian);
                    count = 0;
                    for (int i = 0; i < lE.Count; i++)
                    {
                        var e = lE[i];
                        for (int j = 0; j < e.nIntPoint; j++, count++)
                        {
                            borderVars[2].Add(model.AddVar(-GRB.INFINITY, GRB.INFINITY, 0, GRB.CONTINUOUS, "topZ" + i.ToString() + "," + j.ToString()));
                            slackVars[2].Add(model.AddVar(-GRB.INFINITY,GRB.INFINITY,0,GRB.CONTINUOUS,"slack"));
                            model.Update();
                            var con = new GRBLinExpr();
                            var node=e.getIntPoint(j);
                            con.AddTerm(1, borderVars[2][count]);   //z
                            con.AddTerm(1, planeVars[2][2]);        //D
                            con.AddTerm(node[0], planeVars[2][0]);  //Ax
                            con.AddTerm(node[1], planeVars[2][1]);  //By
                            model.AddConstr(con, GRB.EQUAL, 0, null);
                            
                            var con2 = new GRBLinExpr();
                            con2.AddTerm(-1, slackVars[2][count]);
                            for (int k = 0; k < nU * nV; k++)
                            {
                                if (jacobian[count, k] != 0)
                                {
                                    con2.AddTerm(jacobian[count, k], vars[k]);
                                }
                            }
                            model.AddConstr(con2, GRB.EQUAL, 0, null);
                            obj.AddTerm(1, slackVars[2][count], slackVars[2][count]);
                            obj.AddTerm(-2, slackVars[2][count], borderVars[2][count]);
                            obj.AddTerm(1, borderVars[2][count], borderVars[2][count]);
                        }
                    }
                    
                    
                    nConstraints = myMasonry.totalNumbrOfRightEdgeIntPoint();
                    jacobian = new SparseDoubleArray(nConstraints, nU * nV);
                    myMasonry.GetJacobianOfRightEdge(jacobian);
                    count = 0;
                    for (int i = 0; i < rE.Count; i++)
                    {
                        var e = rE[i];
                        for (int j = 0; j < e.nIntPoint; j++, count++)
                        {
                            borderVars[3].Add(model.AddVar(-GRB.INFINITY, GRB.INFINITY, 0, GRB.CONTINUOUS, "topZ" + i.ToString() + "," + j.ToString()));
                            slackVars[3].Add(model.AddVar(-GRB.INFINITY,GRB.INFINITY,0,GRB.CONTINUOUS,"slack"));
                            model.Update();
                            var con = new GRBLinExpr();
                            var node=e.getIntPoint(j);
                            con.AddTerm(1, borderVars[3][count]);   //z
                            con.AddTerm(1, planeVars[3][2]);        //D
                            con.AddTerm(node[0], planeVars[3][0]);  //Ax
                            con.AddTerm(node[1], planeVars[3][1]);  //By
                            model.AddConstr(con, GRB.EQUAL, 0, null);
                            
                            var con2 = new GRBLinExpr();
                            con2.AddTerm(-1, slackVars[3][count]);
                            for (int k = 0; k < nU * nV; k++)
                            {
                                if (jacobian[count, k] != 0)
                                {
                                    con2.AddTerm(jacobian[count, k], vars[k]);
                                }
                            }
                            model.AddConstr(con2, GRB.EQUAL, 0, null);
                            obj.AddTerm(1, slackVars[3][count], slackVars[3][count]);
                            obj.AddTerm(-2, slackVars[3][count], borderVars[3][count]);
                            obj.AddTerm(1, borderVars[3][count], borderVars[3][count]);

                        }
                    }
                    
                    //Quadratic constraints  detH>0
                    
                    model.SetObjective(obj);
                    
                    model.Optimize();
                    switch (model.Get(GRB.IntAttr.Status))
                    {
                        case GRB.Status.OPTIMAL:
                            //
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
                    for (int i = 0; i < 4; i++)
                    {
                        double AA = planeVars[i][0].Get(GRB.DoubleAttr.X);
                        double BB = planeVars[i][1].Get(GRB.DoubleAttr.X);
                        double DD = planeVars[i][2].Get(GRB.DoubleAttr.X);
                        double Norm = Math.Sqrt(AA * AA + BB * BB + 1);
                        cuttingPlane[i] = new Plane(AA/Norm, BB/Norm, 1.0/Norm, DD/Norm);
                        switch (i)
                        {
                            case 0:
                                foreach (var e in myMasonry.topEdge)
                                {
                                    e.setPlane(AA, BB, 1, DD);
                                }
                                break;
                            case 1:
                                foreach (var e in myMasonry.bottomEdge)
                                {
                                    e.setPlane(AA, BB, 1, DD);
                                }
                                break;
                            case 2:
                                foreach (var e in myMasonry.leftEdge)
                                {
                                    e.setPlane(AA, BB, 1, DD);
                                }
                                break;
                            case 3:
                                foreach (var e in myMasonry.rightEdge)
                                {
                                    e.setPlane(AA, BB, 1, DD);
                                }
                                break;
                        }

                    }
                }
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
            myMasonry.precompute();
            myMasonry.computeAiryFunction();
            myMasonry.giveEdgeTension(0);
            if (boundary)
            {
                myMasonry.computeAngle();
            }
            Force = new double[nU * nV, 3];

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
            int C1 = 0,C2=0;
            if (boundary)
            {
                T1 = (nParticles - 4) * 3 - 1;
                T2 = nParticles * 3 - 1;
                C1 = 0;
                C2 = nParticles - 4;
                for (int i = 0; i < nParticles; i++)
                {
                    if (i == 0 || i == nU - 1 || i == nU * nV - nU || i == nU * nV - 1)
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

            }
            else
            {
                T1 = (nParticles - boundaryIndex.Count()) * 3 - 1;
                T2 = nParticles * 3 - 1;
                C1 = 0;
                C2 = nParticles - boundaryIndex.Count();
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
                force[i * 3 + 0, 0] = 0;
                force[i * 3 + 1, 0] = 0;
                force[i * 3 + 2, 0] = 0;
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

                if (/*!boundary&&*/(Fx * Fx + Fy * Fy + Fz * Fz) >0.01/* 0.01*/)
                {
                    F[i * 3 + 0, 0] = Fx;
                    F[i * 3 + 1, 0] = Fy;
                    F[i * 3 + 2, 0] = Fz+1;
                    Force[i, 0] = Fx;
                    Force[i, 1] = Fy;
                    Force[i, 2] = Fz+1;
                }
                else
                {
                    F[i * 3 + 0, 0] = 0;
                    F[i * 3 + 1, 0] = 0;
                    F[i * 3 + 2, 0] = 0+1;
                    Force[i, 0] = 0;
                    Force[i, 1] = 0;
                    Force[i, 2] = 0+1;
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
            for (int i = 0; i < 4; i++)
            {
                double[] f = null;
                if (cuttingPlane[i] != null)
                {
                    f = cuttingPlane[i].GetPlaneEquation();
                    //System.Windows.MessageBox.Show(f[0].ToString()+","+f[1].ToString()+","+f[2].ToString()+","+f[3].ToString());
                    var P1 = cuttingPlane[i].PointAt(-100, -100);
                    var P2 = cuttingPlane[i].PointAt(100, -100);
                    var P3 = cuttingPlane[i].PointAt(100, 100);
                    var P4 = cuttingPlane[i].PointAt(-100, 100);
                    args.Display.DrawLine(P1, P2, System.Drawing.Color.LightCoral, 2);
                    args.Display.DrawLine(P2, P3, System.Drawing.Color.LightCoral, 2);
                    args.Display.DrawLine(P3, P4, System.Drawing.Color.LightCoral, 2);
                    args.Display.DrawLine(P4, P1, System.Drawing.Color.LightCoral, 2);
                }
            }
            if (boundary)
            {
                var f = outputNurbs.IsoCurve(0, 0);
                var g = outputNurbs.IsoCurve(1, 0);
                var h = outputNurbs.IsoCurve(1, outputNurbs.Domain(0).Max);
                var c = outputNurbs.IsoCurve(0, outputNurbs.Domain(1).Max);
                args.Display.DrawCurve(f, System.Drawing.Color.Orange, 4);
                args.Display.DrawCurve(g, System.Drawing.Color.Orange, 4);
                args.Display.DrawCurve(h, System.Drawing.Color.Orange, 4);
                args.Display.DrawCurve(c, System.Drawing.Color.Orange, 4);
            }
            args.Display.DrawSurface(airyNurbs, System.Drawing.Color.Blue, 1);
            args.Display.DrawSurface(outputNurbs, System.Drawing.Color.Red, 1);
            args.Display.DrawSurface(xyNurbs, System.Drawing.Color.Green, 1);
            double[][] vec = new double[2][] { new double[3], new double[3] };
            double[] val = new double[2];
            double[] node=null;
            double S = 200;
            Nurbs2x(xyNurbs, x);
            myMasonry.setupNodesFromList(x);
            myMasonry.computeGlobalCoord();
            foreach (var e in myMasonry.edgeList)
            {
                for (int i = 0; i < e.nIntPoint; i++)
                {
                    double[] node2 = e.getIntPoint(i);
                    args.Display.Draw2dText(e.intP[i].tension.ToString("00.00"), System.Drawing.Color.White, new Rhino.Geometry.Point3d(node2[0], node2[1], node2[2]), true);
                }
            }
            
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
                if (!boundary)
                {

                    for (int i = 0; i < e.nBIntPoint; i++)
                    {
                        var d = e.getBIntPoint(i);
                        var d2 = new Rhino.Geometry.Point3d(d[0], d[1], d[2]);
                        args.Display.DrawPoint(d2, Rhino.Display.PointStyle.X, 2, System.Drawing.Color.Red);
                    }
                }                
            }
            if (boundary)
            {
                foreach (var e in myMasonry.topEdge)
                {
                    for (int i = 0; i < e.nIntPoint; i++)
                    {
                        var d = e.getIntPoint(i);
                        var d2 = new Rhino.Geometry.Point3d(d[0], d[1], d[2]);
                        outputP.Add(d2);
                        args.Display.DrawPoint(d2, Rhino.Display.PointStyle.X, 2, System.Drawing.Color.Green);
                    }
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
            if (boundary)
            {
                foreach (var e in myMasonry.edgeList)
                {
                    for (int i = 0; i < e.nIntPoint; i++)
                    {
                        var d = e.getIntPoint(i);
                        var d2 = new Rhino.Geometry.Point3d(d[0], d[1], d[2]);
                        outputP.Add(d2);
                        args.Display.DrawPoint(d2, Rhino.Display.PointStyle.X, 2, System.Drawing.Color.Green);
                    }
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
        bool boundary = false;
        protected override void SolveInstance(Grasshopper.Kernel.IGH_DataAccess DA)
        {
            Object inputGeometry = null;
            if (!DA.GetData(0, ref inputGeometry)) { return; }
            if (!DA.GetData(1, ref boundary)) { boundary = false; }
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
                    Minilla3D.Elements.nurbsElement.border _border=Minilla3D.Elements.nurbsElement.border.None;
                    if (j == 1)
                    {
                        _border = _border | Minilla3D.Elements.nurbsElement.border.Top;
                        {
                            int[] index2 = new int[uDim];
                            for (int l = 0; l < uDim; l++)
                            {
                                index2[l] = (j - 1) * nU + i - 1 + l;
                            }
                            var c = new Minilla3D.Elements.nurbsCurve(uDim, index2, i, uKnot);
                            myMasonry.edgeList.Add(c);
                            myMasonry.topEdge.Add(c);
                        }
                    }
                    if (i == 1)
                    {
                        _border = _border | Minilla3D.Elements.nurbsElement.border.Left;
                        {
                            int[] index2 = new int[vDim];
                            for (int k = 0; k < vDim; k++)
                            {
                                index2[k] = (j - 1 + k) * nU + i - 1;
                            }
                            var c = new Minilla3D.Elements.nurbsCurve(vDim, index2, j, vKnot);
                            myMasonry.edgeList.Add(c);
                            myMasonry.leftEdge.Add(c);
                        }
                    }
                    if (j == nV - vDdim)
                    {
                        _border = _border | Minilla3D.Elements.nurbsElement.border.Bottom;
                        {
                            int[] index2 = new int[uDim];
                            for (int l = 0; l < uDim; l++)
                            {
                                index2[l] = (j - 1 + (vDim - 1)) * nU + i - 1 + l;
                            }
                            var c = new Minilla3D.Elements.nurbsCurve(uDim, index2, i, uKnot);
                            myMasonry.edgeList.Add(c);
                            myMasonry.bottomEdge.Add(c);
                        }
                    }
                    if (i == nU - uDdim)
                    {
                        _border = _border | Minilla3D.Elements.nurbsElement.border.Right;
                        {
                            int[] index2 = new int[vDim];
                            for (int k = 0; k < vDim; k++)
                            {
                                index2[k] = (j - 1 + k) * nU + i - 1 + uDim - 1;

                            }
                            var c = new Minilla3D.Elements.nurbsCurve(vDim, index2, j, vKnot);
                            myMasonry.edgeList.Add(c);
                            myMasonry.rightEdge.Add(c);
                        }
                    }
                    myMasonry.elemList.Add(new Minilla3D.Elements.nurbsElement(uDim, vDim, index, i, j, uKnot, vKnot, _border));
                    switch (_border)
                    {
                        case nurbsElement.border.Left|nurbsElement.border.Top | nurbsElement.border.Right:
                            myMasonry.elemList.Last().stitch( myMasonry.leftEdge.Last(),myMasonry.topEdge.Last(), myMasonry.rightEdge.Last());
                            break;
                        case nurbsElement.border.Left | nurbsElement.border.Bottom | nurbsElement.border.Right:
                            myMasonry.elemList.Last().stitch(myMasonry.leftEdge.Last(), myMasonry.bottomEdge.Last(), myMasonry.rightEdge.Last());
                            break;
                        case nurbsElement.border.Left | nurbsElement.border.Top | nurbsElement.border.Bottom:
                            myMasonry.elemList.Last().stitch(myMasonry.leftEdge.Last(), myMasonry.topEdge.Last(), myMasonry.bottomEdge.Last());
                            break;
                        case nurbsElement.border.Right | nurbsElement.border.Top | nurbsElement.border.Bottom:
                            myMasonry.elemList.Last().stitch(myMasonry.rightEdge.Last(), myMasonry.topEdge.Last(), myMasonry.bottomEdge.Last());
                            break;
                        case nurbsElement.border.Left | nurbsElement.border.Right:
                            myMasonry.elemList.Last().stitch(myMasonry.leftEdge.Last(), myMasonry.rightEdge.Last());
                            break;
                        case nurbsElement.border.Top | nurbsElement.border.Bottom:
                            myMasonry.elemList.Last().stitch(myMasonry.topEdge.Last(), myMasonry.bottomEdge.Last());
                            break;
                        case nurbsElement.border.Left | nurbsElement.border.Top:
                            myMasonry.elemList.Last().stitch(myMasonry.leftEdge.Last(), myMasonry.topEdge.Last());
                            break;
                        case nurbsElement.border.Left | nurbsElement.border.Bottom:
                            myMasonry.elemList.Last().stitch(myMasonry.leftEdge.Last(), myMasonry.bottomEdge.Last());
                            break;
                        case nurbsElement.border.Right | nurbsElement.border.Top:
                            myMasonry.elemList.Last().stitch(myMasonry.rightEdge.Last(), myMasonry.topEdge.Last());
                            break;
                        case nurbsElement.border.Right | nurbsElement.border.Bottom:
                            myMasonry.elemList.Last().stitch(myMasonry.rightEdge.Last(), myMasonry.bottomEdge.Last());
                            break;
                        case nurbsElement.border.Left:
                            myMasonry.elemList.Last().stitch(myMasonry.leftEdge.Last());
                            break;
                        case nurbsElement.border.Right:
                            myMasonry.elemList.Last().stitch(myMasonry.rightEdge.Last());
                            break;
                        case nurbsElement.border.Top:
                            myMasonry.elemList.Last().stitch(myMasonry.topEdge.Last());
                            break;
                        case nurbsElement.border.Bottom:
                            myMasonry.elemList.Last().stitch(myMasonry.bottomEdge.Last());
                            break;
                    }
                }
            }
        }
    }

}
