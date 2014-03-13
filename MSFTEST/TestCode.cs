/* Copyright 2013, Gurobi Optimization, Inc. */

/* This example formulates and solves the following simple QCP model:

     maximize    x
     subject to  x + y + z = 1
                 x^2 + y^2 <= z^2 (second-order cone)
                 x^2 <= yz        (rotated second-order cone)
*/

using System;
using Gurobi;
using System.Windows.Forms;
using System.Drawing;
using System.Collections.Generic;
class qcp_cs
{
    struct cbl
    {
        public int P, Q;
    }
    void visualize()
    {
    }
    unsafe static void Main()
    {
        var form = new Form();
        form.Show();
        form.Width=500;
        form.Height=500;
        var pb=new PictureBox();
        form.Controls.Add(pb);
        pb.Left=5;
        pb.Top=5;
        pb.Width=490;
        pb.Height=490;
        var g = pb.CreateGraphics();
        try
        {
            int n = 11;
            int m = 11;
            int nParticles = n * m;
            double[] X=new double[nParticles], Y=new double[nParticles], Z=new double[nParticles];
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    X[j + i * n] = i;
                    Y[j + i * n] = j;
                    Z[j + i * n] = ((i - (m-1) / 2d) * (i - (m-1) / 2d) + (j - (n-1) / 2d) * (j - (n-1) / 2d))*0.3;
                }
            }
            int nCbls = (n - 1) * m + (m - 1) * n;
            cbl[] cbls=new cbl[nCbls];
            //Connectivity matrix
            double[,] C = new double[nCbls, nParticles];
            for (int i = 0; i < nCbls; i++)
            {
                for (int j = 0; j < nParticles; j++)
                {
                    C[i, j] = 0;
                }
            }
            int count=0;
            for (int j = 0; j < m; j++)
            {
                for (int i = 0; i < n - 1; i++)
                {
                    C[count, j * n + i] = 1;
                    C[count, j * n + i + 1] = -1;
                    cbls[count].P=j * n + i;
                    cbls[count].Q=j * n + i + 1;
                    count++;
                }
            }
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < m - 1; i++)
                {
                    C[count, i * n + j] = 1;
                    C[count, (i+1) * n + j] = -1;
                    cbls[count].P=i * n + j;
                    cbls[count].Q=(i+1) * n + j;
                    count++;
                }
            }
            Action visualize = () =>
                {
                    for (int i = 0; i < nCbls; i++)
                    {
                        float px = (float)X[cbls[i].P];
                        float py = (float)Y[cbls[i].P];
                        float pz = (float)Z[cbls[i].P];
                        float qx = (float)X[cbls[i].Q];
                        float qy = (float)Y[cbls[i].Q];
                        float qz = (float)Z[cbls[i].Q];
                        float s = (float)Math.Sqrt(3) / 2f;
                        float S = 10;
                        g.DrawLine(System.Drawing.Pens.Black, new PointF(S * s * (px - py) + 250, S * (0.5f * (px + py) + pz) + 50), new PointF(S * s * (qx - qy) + 250, S * (0.5f * (qx + qy) + qz) + 50));
                        //g.DrawLine(System.Drawing.Pens.Black, new PointF(S *  (px) + 250, S *  (pz) + 50), new PointF(S *  (qx) + 250, S * (qz) + 50));
                    }
                };
            //visualize();
            //fixed points
            int[] fp = new int[4] { 0, n - 1, n * m - 1 - n + 1, n * m - 1 };
            //Force Density
            double[,] Q = new double[nCbls, nCbls];
            for (int i = 0; i < nCbls; i++)
            {
                for (int j = 0; j < nCbls; j++)
                {
                    Q[i, j] = 0;
                }
            }
            for (int i = 0; i < nCbls; i++)
            {
                Q[i, i] = 1;
            }
            //Create Coefficient matrix
            double[,] D = new double[nParticles, nParticles];

            fixed (double* _ptr1 = &C[0, 0], _ptr2 = &Q[0, 0])
            {
                double* ptr1 = _ptr1;
                double* ptr2 = _ptr2;
                double* ptr3 = _ptr1;
                for (int i = 0; i < nParticles; i++)
                {
                    for (int j = 0; j < nParticles; j++)
                    {
                        ptr1 = _ptr1 + i;
                        ptr3 = _ptr1 + j;
                        double val = 0;
                        ptr2 = _ptr2;
                        for (int k = 0; k < nCbls; k++)
                        {
                            val += *ptr1 * *ptr2 * *ptr3;
                            ptr1 += nParticles;
                            ptr2 += 1 + nCbls;
                            ptr3 += nParticles;
                        }
                        D[i,j] = val;
                    }
                }
            }
            //Force
            double[] f = new double[nParticles];
            for (int i = 0; i < nParticles; i++)
            {
                f[i] = 0;
            }
            f[(n*m-1)/2] = 10;
            f[(n * m - 1) / 4] = 50;
            f[(n * m - 1) / 2 + (n * m - 1) / 4] = 20;

            GRBEnv env = new GRBEnv("qcp.log");
            env.Set(GRB.IntParam.Threads, 4);
            GRBModel model = new GRBModel(env);
            double[] lb = new double[nParticles];
            double[] ub = new double[nParticles];
            for (int i = 0; i < nParticles; i++)
            {
                lb[i] = 0;
                ub[i] = 5;
            }
            GRBVar[] vars = model.AddVars(lb, ub, null, null, null);
            model.Update();
            GRBQuadExpr obj = new GRBQuadExpr();
            for (int i = 0; i < nParticles; i++)
            {
                for (int j = 0; j < nParticles; j++)
                {
                    if (Array.IndexOf(fp, i) == -1 && Array.IndexOf(fp, j) == -1)
                    {
                        if (D[i, j] != 0)
                        {
                            obj.AddTerm(D[i, j], vars[i], vars[j]);
                        }
                    }
                }
            }
            for (int i = 0; i < nParticles; i++)
            {
                if (f[i] != 0)
                {
                    obj.AddTerm(-f[i], vars[i]);
                }
            }
            
            model.SetObjective(obj);
            double[] zcoord = new double[4] { 0,0,0,0};
            //Constraints
            for (int i = 0; i < fp.Length; i++)
            {
                GRBLinExpr expr = new GRBLinExpr();
                expr.AddTerm(1.0, vars[fp[i]]);
                model.AddConstr(expr, GRB.EQUAL, zcoord[i], null);
            }

            //Quadratic Constraints
            GRBQuadExpr cond = new GRBQuadExpr();
            for (int i = 0; i < nParticles; i++)
            {
                cond.AddTerm(1, vars[i], vars[i]);
            }
            model.AddQConstr(cond, GRB.LESS_EQUAL, 100, null);
            model.Tune();
            env.WriteParams("Gurobi.Params");
            //model.Write("Gurobi.env");
            model.Optimize();
            if (model.Get(GRB.IntAttr.Status) == GRB.Status.OPTIMAL)
            {
                for (int j = 0; j < nParticles; j++)
                    Z[j] = vars[j].Get(GRB.DoubleAttr.X);
            }
            visualize();
            
            Console.ReadKey();
            model.Dispose();
            env.Dispose();

        }
        catch (GRBException e)
        {
            Console.WriteLine("Error code: " + e.ErrorCode + ". " + e.Message);
        }
    }
}
