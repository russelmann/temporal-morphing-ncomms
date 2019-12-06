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

#include <PhySim/Geometry/Node.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class Node;

	class EnergyElement_FaceNodeColl_Quad : public EnergyElement
	{

	public:
		EnergyElement_FaceNodeColl_Quad(Model_Particles* pModel, const vector<Node*>& m_vnodes, Material* pMaterial);

		virtual ~EnergyElement_FaceNodeColl_Quad(void);

		virtual void ComputeAndStore_Energy();
		virtual void ComputeAndStore_Gradient();
		virtual void ComputeAndStore_Hessian();

	};
}