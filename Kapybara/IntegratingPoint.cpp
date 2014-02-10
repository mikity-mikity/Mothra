#include "stdafx.h"
#include"IntegratingPoint.h"
namespace Minilla3D {
	namespace Elements{
		void ___integratingPoint::initialize(int _elemDim,double* _weight,double* _metric,double* _refMetric,double* _invMetric,double* _refInvMetric,double* _Cauchy,double* _SPK,double *_dv,double *_refDv)
		{
			elemDim=_elemDim;
			weight=_weight;
			metric=_metric;
			refMetric=_refMetric;
			invMetric=_invMetric;
			refInvMetric=_refInvMetric;
			Cauchy=_Cauchy;
			SPK=_SPK;
			dv=_dv;
			refDv=_refDv;
		}
	}
}