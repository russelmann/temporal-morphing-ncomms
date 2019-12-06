//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#pragma once

#include <PhySim/CommonIncludes.h>
#include <PhySim/PhySimInterface.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class Node;

	class Model_Particles;

	class EnergyElement : public IEnergyElement
	{
	protected:

		Model_Particles *m_pModel;
		vector<DoFSet*> m_vDoFs;
		Material* m_pMaterial;

		double m_energy;
		VectorXd m_vgradient;			// Full element gradient
		MatrixXd m_mHessian;			// Full element Hessian
		MatrixXp m_mHessianPFree;		// Pointers to global free Hessian
		MatrixXp m_mHessianPFull;		// Pointers to global full Hessian

		Real m_intVolume;

	public:

		EnergyElement(Model_Particles* pModel, Material* pMaterial);

		virtual ~EnergyElement(void);

		inline virtual size_t GetStencilSize() { return this->m_vDoFs.size(); }
		inline virtual const vector<DoFSet*>& GetDoFStencil() const { return this->m_vDoFs; }
		inline virtual void SetDoFStencil(const vector<DoFSet*> vns) { this->m_vDoFs = vns; }

		virtual void Init() {}

		virtual void ComputeAndStore_Energy() {}
		virtual void ComputeAndStore_Gradient() {}
		virtual void ComputeAndStore_Hessian() {}

		virtual void AssembleGlobal_Gradient(VectorXd& vtotalGradient, bool full = false);
		virtual void AssembleGlobal_Hessian(VectorTd& vtotalHessian, bool full = false);

		virtual void AllocateGlobal_Hessian(const CoefMap& mp, bool full = false);
		virtual void AssembleGlobal_FastPreallocatedHessian(bool full = false);

		inline virtual const Real& GetIntegrationVolume() const { return this->m_intVolume; }
		inline virtual void SetIntegrationVolume(const Real& iv) { this->m_intVolume = iv; }

		inline Real GetElementEnergy() const { return this->m_energy; }
		inline const VectorXd& GetElementGradient() const { return this->m_vgradient; }
		inline const MatrixXd& GetElementHessian() const { return this->m_mHessian; }

		inline virtual Material* GetMaterial() { return this->m_pMaterial; }
		inline virtual const Material* GetMaterial() const { return this->m_pMaterial; }
		inline virtual void SetMaterial(Material* pMat) { this->m_pMaterial = pMat; }

	};
}
