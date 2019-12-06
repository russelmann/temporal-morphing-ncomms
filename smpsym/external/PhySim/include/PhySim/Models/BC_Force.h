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

#include <PhySim/Models/BCondition.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class BC_Force : public BCondition
	{
	protected:
		EnergyElement_Force*	m_pEleForce;

	public:
		BC_Force(Model* pModel, const BCSetup& bcSetup);
		virtual ~BC_Force(void);

		virtual void Init() override;
		virtual void Apply() override;

	};
}
