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

	class EnergyElement_Force : public EnergyElement
	{

	protected:
		vector<VectorXd>		m_vforces;

	public:
		EnergyElement_Force(Model_Particles* pModel, const vector<DoFSet*>& vdofs);
		virtual ~EnergyElement_Force(void);

		virtual void Init();

		virtual const vector<VectorXd>& GetForces() const { return this->m_vforces; }
		virtual void SetForces(const vector<VectorXd>& vf) { this->m_vforces = vf; }

		virtual void ComputeAndStore_Energy();
		virtual void ComputeAndStore_Gradient();
	};

}


