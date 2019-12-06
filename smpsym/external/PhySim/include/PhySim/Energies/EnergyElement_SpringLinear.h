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

#include <PhySim/Energies/EnergyElement.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class EnergyElement_SpringLinear : public EnergyElement
	{

	protected:

		Edge* m_pEdge;
		Real m_restLength;

	public:
		EnergyElement_SpringLinear(Model_Particles* pModel, Edge* pEdge, Material* pMaterial);
		virtual ~EnergyElement_SpringLinear(void);

		virtual void Init();

		virtual Real ComputeRestLength();
		inline virtual Real GetRestLength() const { return this->m_restLength; }
		inline virtual void SetRestLength(Real rl) { this->m_restLength = rl; }

		// This element reimplements gradient and HessianFull computation
		// to take advantage of the symmetries of the gradient and the
		// HessianFull to optimize assembly.

		virtual void AssembleGlobal_Gradient(VectorXd& vtotalGradient, bool full = false);
		virtual void AssembleGlobal_Hessian(VectorTd & vtotalHessian, bool full = false);
		virtual void AssembleGlobal_FastPreallocatedHessian(bool full = false);

		virtual void ComputeAndStore_Energy();
		virtual void ComputeAndStore_Gradient();
		virtual void ComputeAndStore_Hessian();

	};
}

