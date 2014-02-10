#pragma once
#include<iostream>
#include"Elements.h"
#include"NurbsSurface.h"
#include"IntegratingPoint.h"
#include"Kapybara.h"
using namespace System;
namespace Minilla3D {
	namespace Elements{
	public ref class managedElement
		{
		protected:
			iElement* _elem;
		public:
			int nIntPoint;
			int nNode;
			int elemDim;
		private:
			Minilla3D::Materials::material ^constitutive;
		public:
			void initialize(int _nNode,int _elemDim,int _nIntPoint,Minilla3D::Elements::iElement *e)
			{
				_elem=e;
				nNode=_nNode;
				elemDim=_elemDim;
				nIntPoint=_nIntPoint;
				Materials::defaultMaterial^ dM=gcnew Materials::defaultMaterial();
			    constitutive=dM->getMaterial();
			}
		public :
		public:
			managedElement()
			{
			}
			managedElement(int _nNode,int _elemDim,int _nIntPoint,Minilla3D::Elements::iElement *e)
			{
				initialize(_nNode,_elemDim,_nIntPoint,e);
			}
			managedElement(array<int,1> ^_index,int _nNode,int _elemDim,int _nIntPoint,Minilla3D::Elements::iElement *e)
			{
				initialize(_nNode,_elemDim,_nIntPoint,e);
				this->setupIndex(_index);
			}
			property int default[ int ]
			{
			  int get( int index )
			  {
				  return _elem->getIndex(index);
			  }
			}
			property double default[ int ,int]
			{
			  double get( int i,int j)
			  {
				  return _elem->getNode(i,j);
			  }
			}
			void setMaterial(Minilla3D::Materials::material ^_m)
			{
				constitutive=_m;
			}
			void setRefVolume(double v)
			{
				_elem->setRefVolume(v);
			}
			double getVolume()
			{
				return _elem->Volume;
			}
			void getGlobalCoord(array<double,1> ^out,int num)
			{
				pin_ptr<double> ptr=&out[0];
				_elem->copyGlobalCoord(ptr,num);
				ptr=nullptr;
			}
			void getBaseVectors(array<double,2> ^out,int num)
			{
				pin_ptr<double> ptr=&out[0,0];
				_elem->copyBaseVectors(ptr,num);
				ptr=nullptr;
			}
			double getEigenVectors(array<double,2> ^out,array<double,1> ^out2,int num)
			{
				pin_ptr<double> ptr=&out[0,0];
				pin_ptr<double> ptr2=&out2[0];
				return _elem->copyEigenVectors(ptr,ptr2,num);
				ptr=nullptr;
				ptr2=nullptr;
			}
			void getGradient(array<double,1> ^out)
			{
				pin_ptr<double> ptr=&out[0];
				_elem->copyGradient(ptr,nNode*__DIM);
				ptr=nullptr;
			}
			void getStress(array<double,2> ^out,int n)
			{
				pin_ptr<double> ptr=&out[0,0];
				_elem->copyStress(ptr,n);
				ptr=nullptr;
			}
			void computeGlobalCoord()
			{
				_elem->computeGlobalCoord();
			}
			void computeBaseVectors()
			{
				_elem->computeBaseVectors();
			}
			void computeEigenVectors()
			{
				_elem->computeEigenVectors();
			}
			void computeMetric()
			{
				_elem->computeMetric();
			}
			void memoryMetric()
			{
				_elem->memoryMetric();
			}
			void computeVolume()
			{
				_elem->computeVolume();
			}
			void memoryVolume()
			{
				_elem->memoryVolume();
			}
			void computeStress()
			{
				for(int i=0;i<nIntPoint;i++)
				{
					this->constitutive(_elem->iP[i],_elem->Volume,_elem->refVolume);
				}
			}
			void computeGradient()
			{
				_elem->computeGradient();
			}
			void computeHessEd()
			{
				_elem->computeHessEd();
			}
			void setupNodes(array<double,1> ^x)
			{
				pin_ptr<double> _ptr1=&x[0];
				if(x->Length!=nNode*__DIM)
				{
					throw gcnew ArgumentOutOfRangeException();
				}
				_elem->setupNodes(_ptr1);
				_ptr1=nullptr;
			}
			void setupNodesFromList(array<double,1> ^x)
			{
				pin_ptr<double> _ptr1=&x[0];
				_elem->setupNodesFromList(_ptr1);
				_ptr1=nullptr;
			}
			void setupNodesFromList(double *_ptr1)
			{
				_elem->setupNodesFromList(_ptr1);
			}
			void setupIndex(array<int,1> ^f)
			{
				if(f->Length!=nNode)
				{
					throw gcnew System::ArgumentException();
				}else
				{
					pin_ptr<int> _ptr1=&f[0];
					_elem->setupIndex(_ptr1);
					_ptr1=nullptr;
				}
			}
			void mergeGradient(double* ptr)
			{
				_elem->mergeGradient(ptr);
			}
			void mergeHessian(ShoNS::Array::SparseDoubleArray ^hess,int dim)
			{
				_elem->mergeHessian(hess,dim);
			}
			void createStar(array<System::Collections::Generic::List<array<double,1>^>^,1> ^f)
			{
				_elem->collectStar(f);
			}

		};

		public ref class N22D2:public managedElement
		{
			static int ___uDim=2;
			static int ___vDim=2;
		public:
			N22D2(int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<2,2>(uNum,vNum,uk,vk));
				uk=nullptr;
				vk=nullptr;
			}
			N22D2(array<int,1>^_index,int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<2,2>(uNum,vNum,uk,vk));
				this->setupIndex(_index);
				uk=nullptr;
				vk=nullptr;
			}
			~N22D2()
			{
				this->!N22D2();
			}
			!N22D2()
			{
				delete(_elem);
			}
		};
		public ref class N23D2:public managedElement
		{
			static int ___uDim=2;
			static int ___vDim=3;
		public:
			N23D2(int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<2,3>(uNum,vNum,uk,vk));
				uk=nullptr;
				vk=nullptr;
			}
			N23D2(array<int,1>^_index,int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<2,3>(uNum,vNum,uk,vk));
				this->setupIndex(_index);
				uk=nullptr;
				vk=nullptr;
			}
			~N23D2()
			{
				this->!N23D2();
			}
			!N23D2()
			{
				delete(_elem);
			}
		};
		public ref class N32D2:public managedElement
		{
			static int ___uDim=3;
			static int ___vDim=2;
		public:
			N32D2(int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<3,2>(uNum,vNum,uk,vk));
				uk=nullptr;
				vk=nullptr;
			}
			N32D2(array<int,1>^_index,int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<3,2>(uNum,vNum,uk,vk));
				this->setupIndex(_index);
				uk=nullptr;
				vk=nullptr;
			}
			~N32D2()
			{
				this->!N32D2();
			}
			!N32D2()
			{
				delete(_elem);
			}
		};
		public ref class N24D2:public managedElement
		{
			static int ___uDim=2;
			static int ___vDim=4;
		public:
			N24D2(int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<2,4>(uNum,vNum,uk,vk));
				uk=nullptr;
				vk=nullptr;
			}
			N24D2(array<int,1>^_index,int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<2,4>(uNum,vNum,uk,vk));
				this->setupIndex(_index);
				uk=nullptr;
				vk=nullptr;
			}
			~N24D2()
			{
				this->!N24D2();
			}
			!N24D2()
			{
				delete(_elem);
			}
		};
		public ref class N42D2:public managedElement
		{
			static int ___uDim=4;
			static int ___vDim=2;
		public:
			N42D2(int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<4,2>(uNum,vNum,uk,vk));
				uk=nullptr;
				vk=nullptr;
			}
			N42D2(array<int,1>^_index,int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<4,2>(uNum,vNum,uk,vk));
				this->setupIndex(_index);
				uk=nullptr;
				vk=nullptr;
			}
			~N42D2()
			{
				this->!N42D2();
			}
			!N42D2()
			{
				delete(_elem);
			}
		};
		public ref class N25D2:public managedElement
		{
			static int ___uDim=2;
			static int ___vDim=5;
		public:
			N25D2(int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<2,5>(uNum,vNum,uk,vk));
				uk=nullptr;
				vk=nullptr;
			}
			N25D2(array<int,1>^_index,int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<2,5>(uNum,vNum,uk,vk));
				this->setupIndex(_index);
				uk=nullptr;
				vk=nullptr;
			}
			~N25D2()
			{
				this->!N25D2();
			}
			!N25D2()
			{
				delete(_elem);
			}
		};
		public ref class N52D2:public managedElement
		{
			static int ___uDim=5;
			static int ___vDim=2;
		public:
			N52D2(int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<5,2>(uNum,vNum,uk,vk));
				uk=nullptr;
				vk=nullptr;
			}
			N52D2(array<int,1>^_index,int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<5,2>(uNum,vNum,uk,vk));
				this->setupIndex(_index);
				uk=nullptr;
				vk=nullptr;
			}
			~N52D2()
			{
				this->!N52D2();
			}
			!N52D2()
			{
				delete(_elem);
			}
		};
		public ref class N33D2:public managedElement
		{
			static int ___uDim=3;
			static int ___vDim=3;
		public:
			N33D2(int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<3,3>(uNum,vNum,uk,vk));
				uk=nullptr;
				vk=nullptr;
			}
			N33D2(array<int,1>^_index,int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<3,3>(uNum,vNum,uk,vk));
				this->setupIndex(_index);
				uk=nullptr;
				vk=nullptr;
			}
			~N33D2()
			{
				this->!N33D2();
			}
			!N33D2()
			{
				delete(_elem);
			}
		};
		public ref class N34D2:public managedElement
		{
			static int ___uDim=3;
			static int ___vDim=4;
		public:
			N34D2(int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<3,4>(uNum,vNum,uk,vk));
				uk=nullptr;
				vk=nullptr;
			}
			N34D2(array<int,1>^_index,int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<3,4>(uNum,vNum,uk,vk));
				this->setupIndex(_index);
				uk=nullptr;
				vk=nullptr;
			}
			~N34D2()
			{
				this->!N34D2();
			}
			!N34D2()
			{
				delete(_elem);
			}
		};
		public ref class N43D2:public managedElement
		{
			static int ___uDim=4;
			static int ___vDim=3;
		public:
			N43D2(int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<4,3>(uNum,vNum,uk,vk));
				uk=nullptr;
				vk=nullptr;
			}
			N43D2(array<int,1>^_index,int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<4,3>(uNum,vNum,uk,vk));
				this->setupIndex(_index);
				uk=nullptr;
				vk=nullptr;
			}
			~N43D2()
			{
				this->!N43D2();
			}
			!N43D2()
			{
				delete(_elem);
			}
		};
		public ref class N35D2:public managedElement
		{
			static int ___uDim=3;
			static int ___vDim=5;
		public:
			N35D2(int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<3,5>(uNum,vNum,uk,vk));
				uk=nullptr;
				vk=nullptr;
			}
			N35D2(array<int,1>^_index,int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<3,5>(uNum,vNum,uk,vk));
				this->setupIndex(_index);
				uk=nullptr;
				vk=nullptr;
			}
			~N35D2()
			{
				this->!N35D2();
			}
			!N35D2()
			{
				delete(_elem);
			}
		};
		public ref class N53D2:public managedElement
		{
			static int ___uDim=5;
			static int ___vDim=3;
		public:
			N53D2(int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<5,3>(uNum,vNum,uk,vk));
				uk=nullptr;
				vk=nullptr;
			}
			N53D2(array<int,1>^_index,int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<5,3>(uNum,vNum,uk,vk));
				this->setupIndex(_index);
				uk=nullptr;
				vk=nullptr;
			}
			~N53D2()
			{
				this->!N53D2();
			}
			!N53D2()
			{
				delete(_elem);
			}
		};
		public ref class N44D2:public managedElement
		{
			static int ___uDim=4;
			static int ___vDim=4;
		public:
			N44D2(int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<4,4>(uNum,vNum,uk,vk));
				uk=nullptr;
				vk=nullptr;
			}
			N44D2(array<int,1>^_index,int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<4,4>(uNum,vNum,uk,vk));
				this->setupIndex(_index);
				uk=nullptr;
				vk=nullptr;
			}
			~N44D2()
			{
				this->!N44D2();
			}
			!N44D2()
			{
				delete(_elem);
			}
		};
		public ref class N45D2:public managedElement
		{
			static int ___uDim=4;
			static int ___vDim=5;
		public:
			N45D2(int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<4,5>(uNum,vNum,uk,vk));
				uk=nullptr;
				vk=nullptr;
			}
			N45D2(array<int,1>^_index,int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<4,5>(uNum,vNum,uk,vk));
				this->setupIndex(_index);
				uk=nullptr;
				vk=nullptr;
			}
			~N45D2()
			{
				this->!N45D2();
			}
			!N45D2()
			{
				delete(_elem);
			}
		};
		public ref class N54D2:public managedElement
		{
			static int ___uDim=5;
			static int ___vDim=4;
		public:
			N54D2(int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<5,4>(uNum,vNum,uk,vk));
				uk=nullptr;
				vk=nullptr;
			}
			N54D2(array<int,1>^_index,int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<5,4>(uNum,vNum,uk,vk));
				this->setupIndex(_index);
				uk=nullptr;
				vk=nullptr;
			}
			~N54D2()
			{
				this->!N54D2();
			}
			!N54D2()
			{
				delete(_elem);
			}
		};
		public ref class N55D2:public managedElement
		{
			static int ___uDim=5;
			static int ___vDim=5;
		public:
			N55D2(int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<5,5>(uNum,vNum,uk,vk));
				uk=nullptr;
				vk=nullptr;
			}
			N55D2(array<int,1>^_index,int uNum,int vNum,array<double,1>^uKnot,array<double,1>^vKnot)
			{
				pin_ptr<double> uk=&uKnot[0];
				pin_ptr<double>	vk=&vKnot[0];
				this->initialize(___uDim*___vDim,2,(___uDim+1)*(___vDim+1),new nurbsSurfaceElement<5,5>(uNum,vNum,uk,vk));
				this->setupIndex(_index);
				uk=nullptr;
				vk=nullptr;
			}
			~N55D2()
			{
				this->!N55D2();
			}
			!N55D2()
			{
				delete(_elem);
			}
		};
	}
}