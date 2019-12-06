#pragma once
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

	class EnergyElement_SpringLinear_MPoly : public EnergyElement_SpringLinear
	{

	protected:

		Edge* m_pEdge;
		Real m_restLength;
		MatrixXd m_mcoeff;

	public:
		EnergyElement_SpringLinear_MPoly(Model_Particles* pModel, Edge* pEdge, Material* pMaterial);
		virtual ~EnergyElement_SpringLinear_MPoly(void);

		virtual void Init();

		virtual Real ComputeRestLength();

		inline virtual Real GetRestLength() const { return this->m_restLength; }
		inline virtual void SetRestLength(Real rl) { this->m_restLength = rl; }

		inline virtual const MatrixXd& GetPolyCoeff() const { return this->m_mcoeff; }
		inline virtual void SetPolyCoeff(const MatrixXd& vc) { this->m_mcoeff = vc; }

		virtual void AssembleGlobal_Gradient(VectorXd& vtotalGradient, bool full = false) { EnergyElement::AssembleGlobal_Gradient(vtotalGradient, full); }
		virtual void AssembleGlobal_Hessian(VectorTd & vtotalHessian, bool full = false) { EnergyElement::AssembleGlobal_Hessian(vtotalHessian, full); }
		virtual void AssembleGlobal_FastPreallocatedHessian(bool full = false) { EnergyElement::AssembleGlobal_FastPreallocatedHessian(full); }

		virtual void ComputeAndStore_Energy();
		virtual void ComputeAndStore_Gradient();
		virtual void ComputeAndStore_Hessian();

	};
}

