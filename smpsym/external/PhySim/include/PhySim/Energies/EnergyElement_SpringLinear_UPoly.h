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

#include <PhySim/Energies/EnergyElement_SpringLinear.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class EnergyElement_SpringLinear_UPoly : public EnergyElement_SpringLinear
	{
	protected:
		VectorXd m_vcoeff;

	public:
		EnergyElement_SpringLinear_UPoly(Model_Particles* pModel, Edge* pEdge, Material* pMaterial);
		virtual ~EnergyElement_SpringLinear_UPoly(void);

		inline virtual const VectorXd& GetPolyCoeff() const { return this->m_vcoeff; }
		inline virtual void SetPolyCoeff(const VectorXd& vc) { this->m_vcoeff = vc; }

		virtual void AssembleGlobal_Gradient(VectorXd& vtotalGradient, bool full = false) { EnergyElement::AssembleGlobal_Gradient(vtotalGradient, full); }
		virtual void AssembleGlobal_Hessian(VectorTd& vtotalHessian, bool full = false) { EnergyElement::AssembleGlobal_Hessian(vtotalHessian, full); }
		virtual void AssembleGlobal_FastPreallocatedHessian(bool full = false) { EnergyElement::AssembleGlobal_FastPreallocatedHessian(full); }

		virtual void ComputeAndStore_Energy();
		virtual void ComputeAndStore_Gradient();
		virtual void ComputeAndStore_Hessian();

	};
}

