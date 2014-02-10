#include "stdafx.h"
#include"Kapybara.h"
#include"IntegratingPoint.h"
#include"Material.h"

namespace Minilla3D
{
	namespace Materials
	{
		defaultMaterial::defaultMaterial()
		{
			m_mat=gcnew material(this, &defaultMaterial::Compute);
		}
		void  defaultMaterial::Compute(Elements::___integratingPoint& iP,double Volume,double refVolume)
		{
			for(int i=0;i<iP.elemDim*iP.elemDim;i++)
			{
				iP.Cauchy[i]=0.5*iP.invMetric[i];
			}
			double J=(*iP.dv)/(*iP.refDv);
			int nd=iP.elemDim*iP.elemDim;
			if(J<=0)
			{
				for(int i=0;i<nd;i++)
				{
					iP.SPK[i]=0;
					iP.Cauchy[i]=0;
				}
			}else
			{
				for(int i=0;i<nd;i++)
				{
					iP.SPK[i]=0.5*J*iP.invMetric[i];
				}
			}
		}
		material^ defaultMaterial::getMaterial()
		{
			return m_mat;
		}

		///////////////////////////
		stVenantMaterial::stVenantMaterial()
		{
			Young=1.0;
			Poisson=0.2;
			Alpha=0.5;
			m_mat=gcnew material(this, &stVenantMaterial::Compute);
		}
		void  stVenantMaterial::Compute(Elements::___integratingPoint& iP,double Volume,double refVolume)
		{
			double J=(*iP.dv)/(*iP.refDv);
			if(J<=0)
			{
				int nd2=iP.elemDim*iP.elemDim;
				for(int i=0;i<nd2;i++)
				{
					iP.SPK[i]=0;
					iP.Cauchy[i]=0;
				}
			}else{
				double _m=m_Mu*2.0;
				double _l=m_Lambda;
				double _a=m_Alpha;
				int& nD=iP.elemDim;
				if(nD==2)
				{

					for(int i=0;i<2;i++)
					{
						for(int j=0;j<2;j++)
						{
							double val=0;
							for(int k=0;k<2;k++)
							{
								for(int l=0;l<2;l++)
								{
									val+=0.25*(iP.metric[l+k*2]-iP.refMetric[l+k*2])*(_l*iP.refInvMetric[l+k*2]*iP.refInvMetric[j+i*2]+_m*iP.refInvMetric[i+k*2]*iP.refInvMetric[j+l*2]);
								}
							}
							iP.SPK[j+i*2]=val;
						}
					}
					//Eigen value decomposition
					//0. Matrix
					double a=iP.SPK[0]*iP.refMetric[0]+iP.SPK[1]*iP.refMetric[2];
					double c=iP.SPK[0]*iP.refMetric[1]+iP.SPK[1]*iP.refMetric[3];//transpose
					double b=iP.SPK[2]*iP.refMetric[0]+iP.SPK[3]*iP.refMetric[2];//transpose
					double d=iP.SPK[2]*iP.refMetric[1]+iP.SPK[3]*iP.refMetric[3];
					//1. EigenValues
					double f=a+d;double s=std::sqrt(f*f-4*(a*d-b*c));
					double l1=(f+s)/2;
					double l2=(f-s)/2;

					//2. EigenVectors
					double P11=0,P21=0,P12=0,P22=0;
					double det=0;
					if(c!=0)
					{
						P11=l1-d;
						P21=c;

						P12=l2-d;
						P22=c;
					}else if(b!=0)
					{
						P11=b;
						P21=l1-a;

						P12=b;
						P22=l2-a;
					}else
					{
						P11=1;
						P21=0;

						P12=0;
						P22=1;
					}
					det=P11*P22-P21*P12;
				
				
					double p11=0,p12=0,p21=0,p22=0,A11=0,A12=0,A21=0,A22=0;
					double F11=0,F12=0,F21=0,F22=0;
					if(det!=0)
					{
						//3. inverse of eigenvectors
						p11=P22/det;
						p12=-P12/det;
						p21=-P21/det;
						p22=P11/det;
						//4. new eigenvalue
						if(l1<=0)l1=l1*Alpha*2.;else l1=l1*(1-Alpha)*2.;
						if(l2<=0)l2=l2*Alpha*2.;else l2=l2*(1-Alpha)*2.;
						A11=l1;
						A12=0;
						A21=0;
						A22=l2;
						//5. new Transform Matrix“YŽš‚Íã¨‰º‚Ì‡ƒxƒNƒgƒ‹‚Í‰E‚©‚ç‚©‚¯‚éII
						F11=P11*A11*p11+P11*A12*p21+P12*A21*p11+P12*A22*p21;
						F12=P11*A11*p12+P11*A12*p22+P12*A21*p12+P12*A22*p22;
						F21=P21*A11*p11+P21*A12*p21+P22*A21*p11+P22*A22*p21;
						F22=P21*A11*p12+P21*A12*p22+P22*A21*p12+P22*A22*p22;
					}else
					{
						F11=a;
						F12=b;
						F21=c;
						F22=d;
					}
				
					//6. ‚Q‚Â–Ú‚Ì“YŽš‚ðã‚°‚é (F12‚ÆF21‚ª“ü‚ê‘Ö‚í‚Á‚Ä‚¢‚é‚Ì‚Í“]’u
					iP.SPK[0]=F11*iP.refInvMetric[0]+F21*iP.refInvMetric[2];
					iP.SPK[1]=F11*iP.refInvMetric[1]+F21*iP.refInvMetric[3];
					iP.SPK[2]=F12*iP.refInvMetric[0]+F22*iP.refInvMetric[2];
					iP.SPK[3]=F12*iP.refInvMetric[1]+F22*iP.refInvMetric[3];

					iP.Cauchy[0]=iP.SPK[0]/J;
					iP.Cauchy[1]=iP.SPK[1]/J;
					iP.Cauchy[2]=iP.SPK[2]/J;
					iP.Cauchy[3]=iP.SPK[3]/J;
				}else
				{
					//SPK first
					for(int i=0;i<nD;i++)
					{
						for(int j=0;j<nD;j++)
						{
							double val=0;
							for(int k=0;k<nD;k++)
							{
								for(int l=0;l<nD;l++)
								{
									val+=0.25*(iP.metric[l+k*nD]-iP.refMetric[l+k*nD])*(_l*iP.refInvMetric[l+k*nD]*iP.refInvMetric[j+i*nD]+_m*iP.refInvMetric[i+k*nD]*iP.refInvMetric[j+l*nD]);
								}
							}
							iP.SPK[j+i*nD]=val;
							iP.Cauchy[j+i*nD]=val/J;
						}
					}
				}
			}
		}
		material^ stVenantMaterial::getMaterial()
		{
			return m_mat;
		}

		neoHookeanMaterial::neoHookeanMaterial()
		{
			m_mat=gcnew material(this, &neoHookeanMaterial::Compute);
			mu1=0.1;
			K=10;
		}
		void neoHookeanMaterial::Compute(Elements::___integratingPoint& iP,double Volume,double refVolume)
		{
			double J = (*iP.dv)/(*iP.refDv);
			if(J<=0)
			{
				int nd=iP.elemDim*iP.elemDim;
				for(int i=0;i<nd;i++)
				{
					iP.SPK[i]=0;
					iP.Cauchy[i]=0;
				}
			}else
			{
				int& nD=iP.elemDim;
				double K = m_K;
				double u1 = m_mu1;
				double KK = K / 2. * (J - 1.);

				double JJ=Math::Pow(J,-2.0/nD);
				u1*=JJ;
				KK*=J;
				double I1 = 0;
				//First invariant
				for (int i = 0; i < nD; i++)
				{
					for (int j = 0; j < nD; j++)
					{
						I1 += iP.metric[j+i*nD] * iP.refInvMetric[j+i*nD];
					}
				}
				I1/=(double)nD;
				for (int i = 0; i < nD; i++)
				{
					for (int j = 0; j < nD; j++)
					{
						iP.SPK[j+i*nD] = (iP.refInvMetric[j+i*nD] - I1  * iP.invMetric[j+i*nD]) * u1
							+ KK * iP.invMetric[j+i*nD];
						iP.Cauchy[j+i*nD]=iP.SPK[j+i*nD]/J;
					}
				}
			}

		}
		material^ neoHookeanMaterial::getMaterial()
		{
			return m_mat;
		}
		clarenzMaterial::clarenzMaterial()
		{
			m_mat=gcnew material(this, &clarenzMaterial::Compute);
//			W1=3./4.;
//			W2=1./8.;
//			W3=1./8.;
			WL=1./8.;
			WA=6./8.;
			WC=1./8.;
			T=1.0;
		}
		
		void clarenzMaterial::Compute(Elements::___integratingPoint& iP,double Volume,double refVolume)
		{
			double J = (*iP.dv)/(*iP.refDv);
			if(J<=0)
			{
				int nd=iP.elemDim*iP.elemDim;
				for(int i=0;i<nd;i++)
				{
					iP.SPK[i]=0;
					iP.Cauchy[i]=0;
				}
			}else
			{
				int& nD=iP.elemDim;
				double w1 = m_wL;
				double w2 = m_wA;
				double w3 = m_wC;
				double _t=m_t;
				//double beta=0.0;//
			
				//beta=w1+w2;
				double I1 = 0;
				//First invariant
				for (int i = 0; i < nD; i++)
				{
					for (int j = 0; j < nD; j++)
					{
						I1 += iP.metric[j+i*nD] * iP.refInvMetric[j+i*nD];
					}
				}
				double I2=J*J;
				double coeff=_t*System::Math::Pow(I1*I1/I2-4,_t-1.0);
				w3*=coeff;
				int ee=nD*nD;
				for (int i = 0; i < ee; i++)
				{
					iP.SPK[i]=w1*2*(iP.refInvMetric[i]-iP.invMetric[i]/I2)+w2*2*iP.invMetric[i]*(I2-1/I2)+
						w3*coeff*(4*iP.refInvMetric[i]*I1/I2-2*iP.invMetric[i]*I1*I1/I2);
					iP.Cauchy[i]=iP.SPK[i]/J;
				}
			}
		}
		material^ clarenzMaterial::getMaterial()
		{
			return m_mat;
		}

		mikityMaterial::mikityMaterial()
		{
			m_mat=gcnew material(this, &mikityMaterial::Compute);
			T=1.0;
		}
		
		void mikityMaterial::Compute(Elements::___integratingPoint& iP,double Volume,double refVolume)
		{
			double J = (*iP.dv)/(*iP.refDv);
			if(J<=0)
			{
				int nd=iP.elemDim*iP.elemDim;
				for(int i=0;i<nd;i++)
				{
					iP.SPK[i]=0;
					iP.Cauchy[i]=0;
				}
			}else
			{
				int& nD=iP.elemDim;
				double _t=m_t;
			
				double I1 = 0;
				//First invariant
				for (int i = 0; i < nD; i++)
				{
					for (int j = 0; j < nD; j++)
					{
						I1 += iP.metric[j+i*nD] * iP.refInvMetric[j+i*nD];
					}
				}
				double I2=J*J;

				double coeff=_t*System::Math::Pow(1-4*I2/(I1*I1),_t-1.0);
				double p1=-8*I2/(I1*I1);
				double p2=16*I2/(I1*I1*I1);
				for (int i = 0; i < nD; i++)
				{
					for (int j = 0; j < nD; j++)
					{

						iP.SPK[j+i*nD] =coeff*(iP.invMetric[j+i*nD]*p1
							+iP.refInvMetric[j+i*nD]*p2);
						iP.Cauchy[j+i*nD]=iP.SPK[j+i*nD]/J;
					}
				}
			}
		}
		material^ mikityMaterial::getMaterial()
		{
			return m_mat;
		}

		sanderMaterial::sanderMaterial()
		{
			m_mat=gcnew material(this, &sanderMaterial::Compute);
			mu1=0.1;
		}
		void sanderMaterial::Compute(Elements::___integratingPoint& iP,double Volume,double refVolume)
		{
			int& nD=iP.elemDim;
			double u1=m_mu1;
			double J=(*iP.dv)/(*iP.refDv);
			if(J<=0)
			{
				int nd=iP.elemDim*iP.elemDim;
				for(int i=0;i<nd;i++)
				{
					iP.SPK[i]=0;
					iP.Cauchy[i]=0;
				}
			}else
			{
				for (int i = 0; i < nD; i++)
				{
					for (int j = 0; j < nD; j++)
					{
						double val=0;
						for(int k=0;k<nD;k++)
						{
							for(int l=0;l<nD;l++)
							{
								val+=iP.refMetric[l+k*nD]*iP.invMetric[k+i*nD]*iP.invMetric[l+j*nD]*u1;
							}
						}
						iP.SPK[j+i*nD] = -val;
						iP.Cauchy[j+i*nD]=iP.SPK[j+i*nD]/J;
					}
				}
			}
		}
		material^ sanderMaterial::getMaterial()
		{
			return m_mat;
		}
		harmonicMaterial::harmonicMaterial()
		{
			m_mat=gcnew material(this, &harmonicMaterial::Compute);
			mu1=0.1;
		}
		void harmonicMaterial::Compute(Elements::___integratingPoint& iP,double Volume,double refVolume)
		{
			int& nD=iP.elemDim;
			double u1=m_mu1;
			double J = (*iP.dv)/(*iP.refDv);
			
			if(J<=0)
			{
				int nd=iP.elemDim*iP.elemDim;
				for(int i=0;i<nd;i++)
				{
					iP.SPK[i]=0;
					iP.Cauchy[i]=0;
				}
			}else
			{
				for (int i = 0; i < nD; i++)
				{
					for (int j = 0; j < nD; j++)
					{
						iP.SPK[j+i*nD] = iP.refInvMetric[j+i*nD]*u1;
						iP.Cauchy[j+i*nD]=iP.SPK[j+i*nD]/J;
					}
				}
			}
		}
		material^ harmonicMaterial::getMaterial()
		{
			return m_mat;
		}

		conformalMaterial::conformalMaterial()
		{
			m_mat=gcnew material(this, &conformalMaterial::Compute);
			area=50;
			refArea=50;
		}
		void conformalMaterial::Compute(Elements::___integratingPoint& iP,double Volume,double refVolume)
		{
			int& nD=iP.elemDim;
			int ee=nD*nD;
			double J = (*iP.dv)/(*iP.refDv);
			if(J<=0)
			{
				int nd=iP.elemDim*iP.elemDim;
				for(int i=0;i<nd;i++)
				{
					iP.SPK[i]=0;
					iP.Cauchy[i]=0;
				}
			}else
			{
				for (int i = 0; i < ee; i++)
				{
					iP.SPK[i]=iP.refInvMetric[i]+(area/refArea-2.0)*iP.invMetric[i]*J;
					iP.Cauchy[i] = iP.SPK[i]/J;
				}
			}
		}
		material^ conformalMaterial::getMaterial()
		{
			return m_mat;
		}

		mooneyRivlinMaterial::mooneyRivlinMaterial()
		{
			m_mat=gcnew material(this, &mooneyRivlinMaterial::Compute);
			u1=0.1;
			u2=0.1;
			K=10;
		}
		void  mooneyRivlinMaterial::Compute(Elements::___integratingPoint& iP,double Volume,double refVolume)
		{
			double J = (*iP.dv)/(*iP.refDv);
			if(J<=0)
			{
				int nd=iP.elemDim*iP.elemDim;
				for(int i=0;i<nd;i++)
				{
					iP.SPK[i]=0;
					iP.Cauchy[i]=0;
				}
			}else
			{
				int& nD=iP.elemDim;
				double K = m_K;
				double u1 = m_mu1;
				double u2 = m_mu2;
				double I1 = 0;

				for (int i = 0; i < nD; i++)
				{
					for (int j = 0; j < nD; j++)
					{
						I1 += iP.metric[j+i*nD] * iP.refInvMetric[j+i*nD];
					}
				}
				double f = 0;
				for (int k = 0; k < nD; k++)
				{
					for (int l = 0; l < nD; l++)
					{
						for (int u = 0; u < nD; u++)
						{
							for (int s = 0; s < nD; s++)
							{
								f += iP.refInvMetric[l+k*nD] * iP.refInvMetric[s+u*nD] * iP.metric[u+l*nD] * iP.metric[s+k*nD];
							}
						}
					}
				}
				double I1_3=I1/3.;
				double I1_I1_3=I1*I1/3.;
				double f_3=f/3.;
				double K_2=K/2.*J;
				for (int i = 0; i < nD; i++)
				{
					for (int j = 0; j < nD; j++)
					{
						double t = 0;
						for (int k = 0; k < nD; k++)
						{
							for (int l = 0; l < nD; l++)
							{
								t += iP.refInvMetric[i+k*nD] * iP.refInvMetric[l+j*nD] * iP.metric[l+k*nD];
							}
						}
						iP.SPK[j+i*nD] = (iP.refInvMetric[j+i*nD] - I1_3 * iP.invMetric[j+i*nD]) * u1
							+ (I1 * iP.refInvMetric[j+i*nD] - I1_I1_3 * iP.invMetric[j+i*nD] - t + f_3 * iP.invMetric[j+i*nD]) * u2
								+ K_2 * (J - 1.) * iP.invMetric[j+i*nD];
						iP.Cauchy[j+i*nD]=iP.SPK[j+i*nD]/J;
					}
				}
			}
		}
		material^ mooneyRivlinMaterial::getMaterial()
		{
			return m_mat;
		}

		formFindingMaterial::formFindingMaterial()
		{
			m_mat=gcnew material(this, &formFindingMaterial::Compute);
			Weight=1.0;
			Power=2.0;
		}
		void  formFindingMaterial::Compute(Elements::___integratingPoint& iP,double Volume,double refVolume)
		{
			double J=(*iP.dv)/(*iP.refDv);
            if(J<=0)
			{
				int nd=iP.elemDim*iP.elemDim;
				for(int i=0;i<nd;i++)
				{
					iP.SPK[i]=0;
					iP.Cauchy[i]=0;
				}
			}else
			{
				double n = m_Power * m_Weight * std::pow(Volume, m_Power - 1);
				for(int i=0;i<iP.elemDim*iP.elemDim;i++)
				{
					iP.Cauchy[i]=n*iP.invMetric[i];
				}
				if(iP.refDv!=0)
				{
					for(int i=0;i<iP.elemDim*iP.elemDim;i++)
					{
						iP.SPK[i]=n*iP.invMetric[i]*J;
					}
				}
			}
        }

		material^ formFindingMaterial::getMaterial()
		{
			return m_mat;
		}

		slbeMaterial::slbeMaterial()
		{
			m_mat=gcnew material(this, &slbeMaterial::Compute);
			Weight=1.0;
		}
		void  slbeMaterial::Compute(Elements::___integratingPoint& iP,double Volume,double refVolume)
		{
			double J=(*iP.dv)/(*iP.refDv);
			if(J<=0)
			{
				int nd=iP.elemDim*iP.elemDim;
				for(int i=0;i<nd;i++)
				{
					iP.SPK[i]=0;
					iP.Cauchy[i]=0;
				}
			}else
			{
				double f=(Volume-refVolume)/refVolume;
				double n = 2.*m_Weight *f;
				for(int i=0;i<iP.elemDim*iP.elemDim;i++)
				{
					iP.Cauchy[i]=n*iP.invMetric[i];
				}
				if(iP.refDv!=0)
				{
					for(int i=0;i<iP.elemDim*iP.elemDim;i++)
					{
						iP.SPK[i]=iP.Cauchy[i]*J;
					}
				}
			}
        }

		material^ slbeMaterial::getMaterial()
		{
			return m_mat;
		}

		leastSquaresMaterial::leastSquaresMaterial()
		{
			m_mat=gcnew material(this, &leastSquaresMaterial::Compute);
		}
		void  leastSquaresMaterial::Compute(Elements::___integratingPoint& iP,double Volume,double refVolume)
		{
			double J=(*iP.dv)/(*iP.refDv);
			if(J<=0)
			{
				int nd=iP.elemDim*iP.elemDim;
				for(int i=0;i<nd;i++)
				{
					iP.SPK[i]=0;
					iP.Cauchy[i]=0;
				}
			}else
			{
				for(int i=0;i<iP.elemDim*iP.elemDim;i++)
				{
					iP.Cauchy[i]=(iP.refInvMetric[i]/J-iP.invMetric[i]);
					iP.SPK[i] = iP.Cauchy[i]*J;
				}
			}
		}
		material^ leastSquaresMaterial::getMaterial()
		{
			return m_mat;
		}
	}
}