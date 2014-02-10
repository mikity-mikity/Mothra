#pragma once
#include"Kapybara.h"
#include "IntegratingPoint.h"
namespace Minilla3D
{
	namespace Materials{
		public delegate void material(Minilla3D::Elements::___integratingPoint& iP,double Volume,double refVolume);

		public interface class iMaterial
		{
		public:
			virtual material^ getMaterial();
		};

		public ref class defaultMaterial:public iMaterial
		{
		private:
			material^ m_mat;
			void Compute(Elements::___integratingPoint& iP,double Volume,double refVolume);
		public:
			defaultMaterial();
			virtual material^ getMaterial();
		};
		public ref class stVenantMaterial:public iMaterial
		{
		private:
			double m_Young,m_Poisson,m_Mu,m_Lambda,m_Alpha;
			material^ m_mat;
			void Compute(Elements::___integratingPoint& iP,double Volume,double refVolume);
		public:
			property double Young
			{
				double get()
				{
					return m_Young;
				}
				void set(double value)
				{
					m_Young=value;
					m_Mu=m_Young/(2.*(1.+m_Poisson));
					m_Lambda=m_Young*m_Poisson/((1+m_Poisson)*(1-2*m_Poisson));//•½–Ê‰ž—Í‚Ì‚Æ‚«‚Í(1-2*m_Poisson)‚ð(1-m_Poisson)‚É‚·‚é
				}
			}
			property double Poisson
			{
				double get()
				{
					return m_Poisson;
				}
				void set(double value)
				{
					m_Poisson=value;
					m_Mu=m_Young/(2.*(1.+m_Poisson));
					m_Lambda=m_Young*m_Poisson/((1+m_Poisson)*(1-2*m_Poisson));
				}
			}
			property double Alpha
			{
				double get()
				{
					return m_Alpha;
				}
				void set(double value)
				{
					if(value<0)value=0;
					if(value>1)value=1;
					m_Alpha=value;
				}
			}
			stVenantMaterial();
			virtual material ^getMaterial();
		};
		public ref class neoHookeanMaterial:public iMaterial
		{
		private:
			double m_mu1,m_K,m_Alpha;
		public:
			property double mu1
			{
				double get()
				{
					return m_mu1;
				}
				void set(double value)
				{
					m_mu1 = value;
				}
			}
			property double K
			{
				double get()
				{
					return m_K;
				}
				void set(double value)
				{
					m_K = value;
				}
			}
		    property double Alpha
			{
				double get()
				{
					return m_Alpha;
				}
				void set(double value)
				{
					if(value<0)value=0;
					if(value>1)value=1;
					m_Alpha=value;
				}
			}
		private:
			material^ m_mat;
			void Compute(Elements::___integratingPoint& iP,double Volume,double refVolume);
		public:
			neoHookeanMaterial();
			virtual material ^getMaterial();
		};
		public ref class clarenzMaterial:public iMaterial
		{
		private:
			double m_wL,m_wA,m_wC,m_t;
		public:
			property double WL
			{
				double get()
				{
					return m_wL;
				}
				void set(double value)
				{
					m_wL = value;
				}
			}
			property double WA
			{
				double get()
				{
					return m_wA;
				}
				void set(double value)
				{
					m_wA = value;
				}
			}
			property double WC
			{
				double get()
				{
					return m_wC;
				}
				void set(double value)
				{
					m_wC = value;
				}
			}
			property double T
			{
				double get()
				{
					return m_t;
				}
				void set(double value)
				{
					m_t = value;
				}
			}

		private:
			material^ m_mat;
			void Compute(Elements::___integratingPoint& iP,double Volume,double refVolume);
		public:
			clarenzMaterial();
			virtual material ^getMaterial();
		};
		public ref class mikityMaterial:public iMaterial
		{
		private:
			double m_t;
		public:
			property double T
			{
				double get()
				{
					return m_t;
				}
				void set(double value)
				{
					m_t = value;
				}
			}

		private:
			material^ m_mat;
			void Compute(Elements::___integratingPoint& iP,double Volume,double refVolume);
		public:
			mikityMaterial();
			virtual material ^getMaterial();
		};
		public ref class sanderMaterial:public iMaterial
		{
		private:
			double m_mu1;
		public:
			property double mu1
			{
				double get()
				{
					return m_mu1;
				}
				void set(double value)
				{
					m_mu1 = value;
				}
			}

		private:
			material^ m_mat;
			void Compute(Elements::___integratingPoint& iP,double Volume,double refVolume);
		public:
	        sanderMaterial();
			virtual material ^getMaterial();
		};
		public ref class harmonicMaterial:public iMaterial
		{
		private:
			double m_mu1;
		public:
			property double mu1
			{
				double get()
				{
					return m_mu1;
				}
				void set(double value)
				{
					m_mu1 = value;
				}
			}

		private:
			material^ m_mat;
			void Compute(Elements::___integratingPoint& iP,double Volume,double refVolume);
		public:
	        harmonicMaterial();
			virtual material ^getMaterial();
		};

		public ref class conformalMaterial:public iMaterial
		{
		private:
			double m_area;
			double m_refArea;
		public:
			property double area
			{
				double get()
				{
					return m_area;
				}
				void set(double value)
				{
					m_area = value;
				}
			}
			property double refArea
			{
				double get()
				{
					return m_refArea;
				}
				void set(double value)
				{
					m_refArea = value;
				}
			}

		private:
			material^ m_mat;
			void Compute(Elements::___integratingPoint& iP,double Volume,double refVolume);
		public:
	        conformalMaterial();
			virtual material ^getMaterial();
		};





		public ref class mooneyRivlinMaterial:public iMaterial
		{
		private:
			double m_mu1, m_mu2, m_K;
		public:
			property double u1
			{
				void set(double value)
				{
					m_mu1 = value;
				}
				double get()
				{
					return m_mu1;
				}
			}
			property double u2
			{
				void set(double value)
				{
					m_mu2 = value;
				}
				double get()
				{
					return m_mu2;
				}
			}
			property double K
			{
				void set(double value)
				{
					m_K = value;
				}
				double get()
				{
					return m_K;
				}
			}

		private:
			material^ m_mat;
			void Compute(Elements::___integratingPoint& iP,double Volume,double refVolume);
		public:
			mooneyRivlinMaterial();
			virtual material ^getMaterial();
		};
		public ref class formFindingMaterial:public iMaterial
		{
		private:
			double m_Weight;
			double m_Power;
		public:
			property double Weight
			{
				void set(double value)
				{
					m_Weight = value;
				}
				double get()
				{
					return m_Weight;
				}
			}
			property double Power
			{
				void set(double value)
				{
					m_Power = value;
				}
				double get()
				{
					return m_Power;
				}
			}

		private:
			void Compute(Elements::___integratingPoint& iP,double Volume,double refVolume);
			material^ m_mat;
		public:
			formFindingMaterial();
			virtual material ^getMaterial();
		};
		public ref class slbeMaterial:public iMaterial
		{
		private:
			double m_Weight;
		public:
			property double Weight
			{
				void set(double value)
				{
					m_Weight = value;
				}
				double get()
				{
					return m_Weight;
				}
			}
		private:
			void Compute(Elements::___integratingPoint& iP,double Volume,double refVolume);
			material^ m_mat;
		public:
			slbeMaterial();
			virtual material ^getMaterial();
		};

		public ref class leastSquaresMaterial:public iMaterial
		{
		private:
			void Compute(Elements::___integratingPoint& iP,double Volume,double refVolume);
			material^ m_mat;
		public:
			leastSquaresMaterial();
			virtual material ^getMaterial();
		};
	}
}