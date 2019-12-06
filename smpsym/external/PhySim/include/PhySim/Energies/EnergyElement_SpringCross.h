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

	class EnergyElement_SpringCross : public EnergyElement
	{
	protected:

		Edge* m_pEdge0;
		Edge* m_pEdge1;
		Real m_restStrain;

	public:
		EnergyElement_SpringCross(Model_Particles* pModel, Edge* pEdge0, Edge* pEdge1, Material* pMaterial);
		virtual ~EnergyElement_SpringCross(void);

		virtual void Init();

		virtual Real ComputeStrain(Space s = Space::DEF);
		inline virtual Real GetRestStrain() const { return this->m_restStrain; }
		inline virtual void GetRestStrain(Real rl) { this->m_restStrain = rl; }

		virtual void ComputeAndStore_Energy();
		virtual void ComputeAndStore_Gradient();
		virtual void ComputeAndStore_Hessian();

	};
}

