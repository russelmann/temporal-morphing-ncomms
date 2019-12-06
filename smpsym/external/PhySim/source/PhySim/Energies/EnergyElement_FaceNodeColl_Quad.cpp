//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Energies/EnergyElement_FaceNodeColl_Quad.h>

#include <PhySim/Geometry/Node.h>
#include <PhySim/Geometry/Edge.h>
#include <PhySim/Models/Model_Particles.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	EnergyElement_FaceNodeColl_Quad::EnergyElement_FaceNodeColl_Quad(Model_Particles* pModel, const vector<Node*>& m_vnodes, Material* pMaterial) : EnergyElement(pModel, pMaterial)
	{
		this->m_vDoFs.resize(4);
		this->m_vDoFs[0] = m_vnodes[0]->DoF();
		this->m_vDoFs[1] = m_vnodes[1]->DoF();
		this->m_vDoFs[2] = m_vnodes[2]->DoF();
		this->m_vDoFs[3] = m_vnodes[3]->DoF();

		this->m_vgradient.resize(12);
		this->m_mHessian.resize(12, 12);
		this->m_mHessianPFree.resize(12, 12);
		this->m_mHessianPFull.resize(12, 12);

		this->Init();
	}

	EnergyElement_FaceNodeColl_Quad::~EnergyElement_FaceNodeColl_Quad(void)
	{
		// Nothing to do here...
	}

	void EnergyElement_FaceNodeColl_Quad::ComputeAndStore_Energy()
	{
		const Vector3d& x0 = this->m_vDoFs[0]->GetPosition_x();
		const Vector3d& x1 = this->m_vDoFs[1]->GetPosition_x();
		const Vector3d& x2 = this->m_vDoFs[2]->GetPosition_x();
		const Vector3d& x3 = this->m_vDoFs[3]->GetPosition_x();
		Vector3d e0 = (x1 - x0).normalized();
		Vector3d e1 = (x2 - x0).normalized();
		Vector3d e3 = x3 - x0;
		Real D = e0.cross(e1).dot(e3);
		Real T = (*this->m_pMaterial)[Material::Property::CollT];
		Real K = (*this->m_pMaterial)[Material::Property::CollK];
		if (D > T)
		{
			this->m_energy = 0;
			return; // Oneside
		}

#include "../codegen/FaceNodeCollQuad_Energy.mcg"

		this->m_energy = t47;
	}

	void EnergyElement_FaceNodeColl_Quad::ComputeAndStore_Gradient()
	{
		const Vector3d& x0 = this->m_vDoFs[0]->GetPosition_x();
		const Vector3d& x1 = this->m_vDoFs[1]->GetPosition_x();
		const Vector3d& x2 = this->m_vDoFs[2]->GetPosition_x();
		const Vector3d& x3 = this->m_vDoFs[3]->GetPosition_x();
		Vector3d e0 = (x1 - x0).normalized();
		Vector3d e1 = (x2 - x0).normalized();
		Vector3d e3 = x3 - x0;
		Real D = e0.cross(e1).dot(e3);
		Real T = (*this->m_pMaterial)[Material::Property::CollT];
		Real K = (*this->m_pMaterial)[Material::Property::CollK];
		if (D > T)
		{
			this->m_vgradient.setZero(12);
			return; // Oneside (return)
		}

		Real vgx[12];

#include "../codegen/FaceNodeCollQuad_Gradient.mcg"

		for (int i = 0; i < 12; ++i)
			m_vgradient(i) = vgx[i];
	}

	void EnergyElement_FaceNodeColl_Quad::ComputeAndStore_Hessian()
	{
		const Vector3d& x0 = this->m_vDoFs[0]->GetPosition_x();
		const Vector3d& x1 = this->m_vDoFs[1]->GetPosition_x();
		const Vector3d& x2 = this->m_vDoFs[2]->GetPosition_x();
		const Vector3d& x3 = this->m_vDoFs[3]->GetPosition_x();
		Vector3d e0 = (x1 - x0).normalized();
		Vector3d e1 = (x2 - x0).normalized();
		Vector3d e3 = x3 - x0;
		Real D = e0.cross(e1).dot(e3);
		Real T = (*this->m_pMaterial)[Material::Property::CollT];
		Real K = (*this->m_pMaterial)[Material::Property::CollK];
		if (D > T)
		{
			this->m_mHessian.setZero(12,12);
			return; // Oneside (return)
		}

		Real mHx[12][12];

#include "../codegen/FaceNodeCollQuad_Hessian.mcg"

		for (int i = 0; i < 12; ++i)
			for (int j = 0; j < 12; ++j)
				m_mHessian(i,j) = mHx[i][j];
	}

}