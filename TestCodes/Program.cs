using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ShoNS.Array;
namespace TestCodes
{
    class Program
    {
        static void Main(string[] args)
        {
            
            List<double> a = new List<double>();
            var d = a.Last();

            var v=new DoubleArray(5);
            Console.WriteLine(v.size0);
            Console.WriteLine(v.size1);
            v[0] = 10;
            v[1] = 20;
            v[2] = 30;

            System.Console.WriteLine((v*v.T).ToString());
            var M = new SparseDoubleArray(5, 5);
            M[0, 0] = 5;
            M[0, 0] += 5;
            M[4, 2] = 3;
            M[1, 1] = 5;
            M[2, 2] = 3;
            M[2, 2] = 4;
            var t = new DoubleArray(5, 2);
            t[0, 0] = 5;
            t[1, 0] = 5;
            t[2, 0] = 5;
            t[0, 1] = 3;
            t[1, 1] = 3;
            t[2, 1] = 3;
            System.Console.WriteLine(M.ToString());
            System.Console.WriteLine(M.size0.ToString());
            System.Console.WriteLine(M.size1.ToString());
            var f = t.T * M * t;
            System.Console.WriteLine(f);
            
            int[] _index = new int[] { 0, 1, 2, 5, 6, 7, 10, 11, 12 };
            double[] _nodes = new double[60];
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 5; j++)
                {
                    _nodes[(j + i * 5) * 3 + 0] = j+5;
                    _nodes[(j + i * 5) * 3 + 1] = i+5;
                    _nodes[(j + i * 5) * 3 + 2] = 0;
                }
            }
            double[] uKnot = new double[] { 0, 0, 0, 1, 2, 3, 3, 3 };
            double[] vKnot = new double[] { 0, 0, 0, 1, 1, 1 };
            Minilla3D.Elements.N33D2 n33d2 = new Minilla3D.Elements.N33D2(_index, 1, 1, uKnot, vKnot);
            n33d2.setupNodesFromList(_nodes);
            n33d2.computeGlobalCoord();
            double[] g = new double[3];
            double[,] b = new double[2,3];
            for (int i = 0; i < 16; i++)
            {
                n33d2.getGlobalCoord(g, i);
                n33d2.getBaseVectors(b, i);
                System.Console.WriteLine("{0}::{1},{2},{3}", i, g[0], g[1], g[2]);
                System.Console.WriteLine("base_0::{1},{2},{3}", i, b[0, 0], b[0, 1], b[0, 2]);
                System.Console.WriteLine("base_1::{1},{2},{3}", i, b[1, 0], b[1, 1], b[1, 2]);
            }
            System.Console.Read();
        }
    }
}
