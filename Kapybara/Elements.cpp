#include "stdafx.h"
#include"Elements.h"
#include"NurbsSurface.h"
#include"ManagedElements.h"
#include"IntegratingPoint.h"
namespace Minilla3D {
	namespace Elements{
#define __nIntPoint (_uDim+1)*(_vDim+1)
#define __nNode (_uDim)*(_vDim)
#define __elemDim 2

		template<int _uDim,int _vDim> double nurbsSurfaceElement<_uDim,_vDim>::_weight[__nIntPoint];
		template<int _uDim,int _vDim> double nurbsSurfaceElement<_uDim,_vDim>::_localCoord[__nIntPoint][__elemDim];
		template<int _uDim,int _vDim> double nurbsSurfaceElement<_uDim,_vDim>::_cu[__elemDim][_uDim+1];
		template<int _uDim,int _vDim> double nurbsSurfaceElement<_uDim,_vDim>::_pu[__elemDim][_uDim+1];
#undef __nIntPoint
#undef __nNode
#undef __elemDim
	}
}
