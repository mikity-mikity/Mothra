using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ShoNS.Array;
namespace mikity.ghComponents
{
    public delegate void DrawViewPortWire(Grasshopper.Kernel.IGH_PreviewArgs args);
    public delegate void UpdateGeometry(double x, double y, double z);
    public delegate void BakeGeometry(Rhino.RhinoDoc doc, Rhino.DocObjects.ObjectAttributes att, List<Guid> obj_ids);
    
}