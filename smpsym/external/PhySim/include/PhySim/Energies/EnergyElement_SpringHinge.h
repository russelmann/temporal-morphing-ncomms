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

	class EnergyElement_SpringHinge : public EnergyElement
	{

	protected:

		Edge* m_pEdge0;
		Edge* m_pEdge1;
		Real m_restAngle;

	public:
		EnergyElement_SpringHinge(Model_Particles* pModel, Edge* pEdge0, Edge* pEdge1, Material* pMaterial);
		virtual ~EnergyElement_SpringHinge(void);

		virtual void Init();

		virtual Real ComputeRestAngle();
		inline virtual Real GetRestAngle() const { return this->m_restAngle; }
		inline virtual void SetRestAngle(Real rl) { this->m_restAngle = rl; }

		virtual void ComputeAndStore_Energy();
		virtual void ComputeAndStore_Gradient();
		virtual void ComputeAndStore_Hessian();

	};
}

