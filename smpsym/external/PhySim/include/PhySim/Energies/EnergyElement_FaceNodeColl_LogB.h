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

	class EnergyElement_FaceNodeColl_LogB : public EnergyElement
	{

	public:
		EnergyElement_FaceNodeColl_LogB(Model_Particles* pModel, const vector<Node*>& m_vnodes, Material* pMaterial) : EnergyElement(pModel, pMaterial)
		{
			this->m_vDoFs.resize(4);
			this->m_vDoFs[0] = m_vnodes[0]->DoF();
			this->m_vDoFs[1] = m_vnodes[1]->DoF();
			this->m_vDoFs[2] = m_vnodes[2]->DoF();
			this->m_vDoFs[3] = m_vnodes[3]->DoF();
		}

		virtual ~EnergyElement_FaceNodeColl_LogB(void)
		{
			// Nothing to do here...
		}

		virtual void ComputeAndStore_Energy();
		virtual void ComputeAndStore_Gradient();
		virtual void ComputeAndStore_Hessian();

	};
}