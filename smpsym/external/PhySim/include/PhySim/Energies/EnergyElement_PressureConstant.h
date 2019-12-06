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

	class EnergyElement_PressureConstant : public EnergyElement
	{

	protected:

		Face_Tri* m_pFace;

	public:
		EnergyElement_PressureConstant(Model_Particles* pModel, Face_Tri* pFace, Material* pMaterial);
		virtual ~EnergyElement_PressureConstant(void);

		virtual void Init();

		virtual void ComputeAndStore_Energy();
		virtual void ComputeAndStore_Gradient();
		virtual void ComputeAndStore_Hessian();

	};
}