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

	class EnergyElement_FaceNodeColl_Step : public EnergyElement
	{
	
	public:
		EnergyElement_FaceNodeColl_Step(Model_Particles* pModel, const vector<Node*>& m_vnodes, Material* pMaterial) : EnergyElement(pModel, pMaterial)
		{
			this->m_vDoFs.resize(4);
			this->m_vDoFs[0] = m_vnodes[0]->DoF();
			this->m_vDoFs[1] = m_vnodes[1]->DoF();
			this->m_vDoFs[2] = m_vnodes[2]->DoF();
			this->m_vDoFs[3] = m_vnodes[3]->DoF();
		}

		virtual ~EnergyElement_FaceNodeColl_Step(void)
		{
			// Nothing to do here...
		}

		virtual void ComputeAndStore_Energy()
		{
			const Vector3d& x0 = this->m_vDoFs[0]->GetPosition_x();
			const Vector3d& x1 = this->m_vDoFs[1]->GetPosition_x();
			const Vector3d& x2 = this->m_vDoFs[2]->GetPosition_x();
			const Vector3d& x3 = this->m_vDoFs[3]->GetPosition_x();
			Vector3d e0 = (x1 - x0).normalized();
			Vector3d e1 = (x2 - x0).normalized();
			Vector3d e3 = x3 - x0;
			Real T = (*this->m_pMaterial)[Material::Property::CollT];

			this->m_energy = (e0.cross(e1).dot(e3) > T) ? 0 : HUGE_VAL;
		}

	};
}

