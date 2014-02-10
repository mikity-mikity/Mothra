using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;
namespace mikity.GeometryProcessing
{
    class face
    {
        public face(int _N, params int[] indices)
        {
            corner = indices;
            N = _N;
        }
        public int N;
        public int[] corner;
        public halfedge firsthalfedge;
    }

    class halfedge
    {
        public vertex P;
        public face owner;
        public halfedge pair,next,prev;
        public bool isNaked
        {
            get
            {
                return pair == null ? true : false;
            }
        }
        public halfedge(vertex _P)
        {
            P = _P;
            if (_P.hf_begin == null) _P.hf_begin = this;
        }
    }
    class vertex
    {
        public int N;
        public double x, y, z;
        public List<halfedge> star=new List<halfedge>();
        public List<halfedge> onering=new List<halfedge>();
        public halfedge hf_begin;
        public halfedge hf_end;
        public bool isNaked
        {
            get
            {
                return hf_begin == hf_end ? false : true;
            }
        }
        public vertex(int _N)
        {
            N = _N;
        }
        public bool isInner
        {
            get
            {
                return onering[0] == onering[onering.Count - 1].next;
            }
        }
        public bool isBoundary
        {
            get
            {
                return onering[0] != onering[onering.Count - 1].next;
            }
        }
    }
    class MeshStructure
    {
        //to get boundary chain
        //boundaryStart->hf->next->P->hf->next->P->....
        public vertex boundaryStart;
        public List<vertex> vertices = new List<vertex>();
        public List<face> faces = new List<face>();
        public List<halfedge> halfedges = new List<halfedge>();
        public List<vertex> innerVertices = new List<vertex>();
        public List<vertex> outerVertices = new List<vertex>();
        private halfedge[,] __halfedgeTable;
        private List<face>[,] _faceTable;
        private orient[] __orientation;
        public int nVertices
        {
            get
            {
                return vertices.Count;
            }
        }
        public int nFaces
        {
            get
            {
                return faces.Count;
            }
        }
        private MeshStructure()
        {
            this.Clear();
        }
        private enum orient
        {
            unknown,positive,negative
        }
        private void Construct(Mesh val)
        {
            int _nVertices = val.Vertices.Count;
            int _nFaces = val.Faces.Count;
            
            __orientation = new orient[_nFaces];
            _faceTable = new List<face>[_nVertices, _nVertices];
            __halfedgeTable = new halfedge[_nVertices, _nVertices];
            
            for (int i = 0; i < __orientation.Count(); i++)
            {
                __orientation[i] = orient.unknown;
            }

            for (int i = 0; i < _nVertices; i++)
            {
                var _v = new vertex(i);
                vertices.Add(_v);
            }

            for (int i = 0; i < _nFaces; i++)
            {
                var f=val.Faces[i];
                var _f = f.IsQuad ? new face(i,f.A, f.B, f.C, f.D) : new face(i,f.A, f.B, f.C);
                faces.Add(_f);
                faceTableAdd(_f);
            }
            //Recursive
            halfEdgeAdd(faces[0]);
            //find pairs
            foreach (var h in halfedges)
            {
                int i = h.P.N;
                int j = h.next.P.N;
                if (__halfedgeTable[i, j] != null) throw new ArgumentOutOfRangeException(";)");
                __halfedgeTable[i, j] = h;
            }
            foreach (var h in halfedges)
            {
                int i = h.P.N;
                int j = h.next.P.N;
                //if boundary edge...
                if (__halfedgeTable[j, i] == null)
                {
                    h.pair = null;
                }
                else
                {
                    h.pair = __halfedgeTable[j, i];
                }
            }
            //post process to find boundary vertices
            foreach (var v in vertices)
            {
                var h=v.hf_begin;
                v.hf_end = h;
                do{
                    if (h.prev.isNaked)
                    {
                        v.hf_end=h.prev;
                        while(!h.isNaked){
                            h = h.pair.next;
                        }
                        v.hf_begin = h;
                        if (this.boundaryStart == null) this.boundaryStart = v;
                            break;
                    }
                    h = h.prev.pair;
                } while (h != v.hf_begin);
            }

            //post process to create stars
            foreach (var v in vertices)
            {
                var h = v.hf_begin;
                v.star.Clear();
                do
                {
                    v.star.Add(h);
                    if (h.prev.isNaked) break;
                    h = h.prev.pair;
                } while (h != v.hf_begin);
            }
            //post process to create onering
            foreach (var v in vertices)
            {
                var h = v.hf_begin;
                v.onering.Clear();
                do
                {
                    do
                    {
                        h = h.next;
                        v.onering.Add(h);
                    } while (h.next.next.P != v);
                    if (h.next.isNaked) break;
                    h = h.next.pair;
                } while (h != v.hf_begin);
            }
            //post process to split the vertices into inner and outer.
            innerVertices.Clear();
            outerVertices.Clear();
            foreach (var v in vertices)
            {
                if (v.hf_begin.isNaked) outerVertices.Add(v); else innerVertices.Add(v);
            }
        }
        private void halfEdgeAdd(face f)
        {
            var _o=orient.unknown;
            for (int i = 0; i < f.corner.Count(); i++)
            {
                int I = f.corner[i];
                int J = (i == f.corner.Count() - 1) ? f.corner[0] : f.corner[i + 1];
                if (_faceTable[I, J].Count == 2)
                {
                    if (_faceTable[I, J][0] == f)
                    {
                        if (__orientation[_faceTable[I, J][1].N] != orient.unknown)
                        {
                            _o = __orientation[_faceTable[I, J][1].N] == orient.positive ? orient.negative : orient.positive;
                        }
                    }
                    if (_faceTable[I, J][1] == f)
                    {
                        if (__orientation[_faceTable[I, J][0].N] != orient.unknown)
                        {
                            _o = __orientation[_faceTable[I, J][0].N] == orient.positive ? orient.negative : orient.positive;
                        }
                    }
                }
                else
                {
                    if (_faceTable[J, I] != null)
                    {
                        if (__orientation[_faceTable[J, I][0].N] != orient.unknown)
                        {
                            _o = __orientation[_faceTable[J, I][0].N];
                        }
                    }
                }                
            }
            __orientation[f.N] = _o == orient.unknown ? orient.positive : _o;
            //register a halfedge
            if (f.corner.Count() == 3 && __orientation[f.N] == orient.positive)
            {
                var he1 = new halfedge(vertices[f.corner[0]]);
                var he2 = new halfedge(vertices[f.corner[1]]);
                var he3 = new halfedge(vertices[f.corner[2]]);
                halfedges.Add(he1);
                halfedges.Add(he2);
                halfedges.Add(he3);
                he1.prev = he3; he1.next = he2; he1.owner = f;
                he2.prev = he1; he2.next = he3; he2.owner = f;
                he3.prev = he2; he3.next = he1; he3.owner = f;
                f.firsthalfedge = he1;
            }

            if (f.corner.Count() == 3 && __orientation[f.N] == orient.negative)
            {
                var he1 = new halfedge(vertices[f.corner[2]]);
                var he2 = new halfedge(vertices[f.corner[1]]);
                var he3 = new halfedge(vertices[f.corner[0]]);
                halfedges.Add(he1);
                halfedges.Add(he2);
                halfedges.Add(he3);
                he1.prev = he3; he1.next = he2; he1.owner = f;
                he2.prev = he1; he2.next = he3; he2.owner = f;
                he3.prev = he2; he3.next = he1; he3.owner = f;
                f.firsthalfedge = he1;
            }

            if (f.corner.Count() == 4 && __orientation[f.N] == orient.positive)
            {
                var he1 = new halfedge(vertices[f.corner[0]]);
                var he2 = new halfedge(vertices[f.corner[1]]);
                var he3 = new halfedge(vertices[f.corner[2]]);
                var he4 = new halfedge(vertices[f.corner[3]]);
                halfedges.Add(he1);
                halfedges.Add(he2);
                halfedges.Add(he3);
                halfedges.Add(he4);
                he1.prev = he4; he1.next = he2; he1.owner = f;
                he2.prev = he1; he2.next = he3; he2.owner = f;
                he3.prev = he2; he3.next = he4; he3.owner = f;
                he4.prev = he3; he4.next = he1; he4.owner = f;
                f.firsthalfedge = he1;
            }

            if (f.corner.Count() == 4 && __orientation[f.N] == orient.negative)
            {
                var he1 = new halfedge(vertices[f.corner[3]]);
                var he2 = new halfedge(vertices[f.corner[2]]);
                var he3 = new halfedge(vertices[f.corner[1]]);
                var he4 = new halfedge(vertices[f.corner[0]]);
                halfedges.Add(he1);
                halfedges.Add(he2);
                halfedges.Add(he3);
                halfedges.Add(he4);
                he1.prev = he4; he1.next = he2; he1.owner = f;
                he2.prev = he1; he2.next = he3; he2.owner = f;
                he3.prev = he2; he3.next = he4; he3.owner = f;
                he4.prev = he3; he4.next = he1; he4.owner = f;
                f.firsthalfedge = he1;
            }

            //list up neighbors that are not oriented
            for (int i = 0; i < f.corner.Count(); i++)
            {
                int I = f.corner[i];
                int J = (i == f.corner.Count() - 1) ? f.corner[0] : f.corner[i + 1];
                if (_faceTable[I, J].Count == 2)
                {
                    if (_faceTable[I, J][0] == f)
                    {
                        if (__orientation[_faceTable[I, J][1].N] == orient.unknown)
                        {
                            halfEdgeAdd(_faceTable[I, J][1]);
                        }
                    }
                    if (_faceTable[I, J][1] == f)
                    {
                        if (__orientation[_faceTable[I, J][0].N] == orient.unknown)
                        {
                            halfEdgeAdd(_faceTable[I, J][0]);
                        }
                    }
                }
                else
                {
                    if (_faceTable[J, I] != null)
                    {
                        if (__orientation[_faceTable[J, I][0].N] == orient.unknown)
                        {
                            halfEdgeAdd(_faceTable[J, I][0]);
                        }
                    }
                }
            }
        }
        private void faceTableAdd(int i, int j, face f)
        {
            if (_faceTable[i, j] == null)
            {
                _faceTable[i, j] = new List<face>();
            }
            _faceTable[i, j].Add(f);
        }
        private void faceTableAdd(face f)
        {
            for (int i = 0; i < f.corner.Count(); i++)
            {
                int I = f.corner[i];
                int J = (i == f.corner.Count() - 1) ? f.corner[0] : f.corner[i + 1];
                faceTableAdd(I, J, f);
            }
        }
        public static MeshStructure CreateFrom(Mesh val)
        {
            var ret=new MeshStructure();
            ret.Construct(val);
            return ret;
        }
        public void Clear()
        {
            vertices.Clear();
            faces.Clear();
            halfedges.Clear();
            innerVertices.Clear();
            outerVertices.Clear();
            boundaryStart = null;
        }
        public void Update(double[,] X)
        {
            foreach (var v in vertices)
            {
                v.x = X[v.N, 0];
                v.y = X[v.N, 1];
                v.z = X[v.N, 2];
            }
        }
    }
}
