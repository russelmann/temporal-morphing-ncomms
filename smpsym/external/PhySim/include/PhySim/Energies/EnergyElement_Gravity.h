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

	class EnergyElement_Gravity : public EnergyElement
	{

	protected:
		Vector3d m_vgravity;

	public:
		EnergyElement_Gravity(Model_Particles* pModel);
		virtual ~EnergyElement_Gravity(void);

		virtual void Init();

		virtual const Vector3d& GetGravityAcceleration() const { return this->m_vgravity; }
		virtual void SetGravityAcceleration(const Vector3d& vg) { this->m_vgravity = vg; }

		virtual void ComputeAndStore_Energy();
		virtual void ComputeAndStore_Gradient();
	};

}


