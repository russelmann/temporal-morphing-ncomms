//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Energies/EnergyElement_SpringHinge.h>

#include <PhySim/Geometry/Node.h>
#include <PhySim/Geometry/Edge.h>
#include <PhySim/Models/Model_Particles.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	EnergyElement_SpringHinge::EnergyElement_SpringHinge(Model_Particles* pModel, Edge* pEdge0, Edge* pEdge1, Material* pMaterial) : EnergyElement(pModel, pMaterial)
	{
		this->m_pEdge0 = pEdge0;
		this->m_pEdge1 = pEdge1;

		this->m_vDoFs.resize(3);
		this->m_vDoFs[0] = pEdge0->GetOrigin()->DoF();
		this->m_vDoFs[1] = pEdge0->GetHead()->DoF();
		this->m_vDoFs[2] = pEdge1->GetHead()->DoF();

		this->m_vgradient.resize(9);
		this->m_mHessian.resize(9, 9);
		this->m_mHessianPFree.resize(9, 9);
		this->m_mHessianPFull.resize(9, 9);

		this->Init();
	}

	EnergyElement_SpringHinge::~EnergyElement_SpringHinge()
	{
		// Nothing to do...
	}

	void EnergyElement_SpringHinge::Init()
	{
		this->m_restAngle = this->ComputeRestAngle();
		this->m_intVolume = 0.5*(this->m_pEdge0->ComputeVolume(Space::MAT) + 
								 this->m_pEdge1->ComputeVolume(Space::MAT));
	}

	Real EnergyElement_SpringHinge::ComputeRestAngle()
	{
		const Vector3d& a0 = this->m_vDoFs[0]->GetPosition_0();
		const Vector3d& b0 = this->m_vDoFs[1]->GetPosition_0();
		const Vector3d& c0 = this->m_vDoFs[2]->GetPosition_0();
		Vector3d n0 = (b0 - a0).normalized();
		Vector3d n1 = (c0 - b0).normalized();
		return acos(n0.dot(n1));
	}

	void EnergyElement_SpringHinge::ComputeAndStore_Energy()
	{
		const Vector3d& x1 = this->m_vDoFs[0]->GetPosition_x();
		const Vector3d& x2 = this->m_vDoFs[1]->GetPosition_x();
		const Vector3d& x3 = this->m_vDoFs[2]->GetPosition_x();

		Real K = (*this->m_pMaterial)[Material::Property::BendingK];
		Real L0 = m_intVolume;
		Real A0 = m_restAngle;

		{
#include "../codegen/SpringHingeEnergy.mcg"
		
			this->m_energy = t37;
		}
	}

	void EnergyElement_SpringHinge::ComputeAndStore_Gradient()
	{
		const Vector3d& x1 = this->m_vDoFs[0]->GetPosition_x();
		const Vector3d& x2 = this->m_vDoFs[1]->GetPosition_x();
		const Vector3d& x3 = this->m_vDoFs[2]->GetPosition_x();

		Real K = (*this->m_pMaterial)[Material::Property::BendingK];
		Real L0 = m_intVolume;
		Real A0 = m_restAngle;

		Real vgx[9];

		{
#include "../codegen/SpringHingeGradient.mcg"
		}

		for (int i = 0; i < 9; ++i)
			m_vgradient(i) = vgx[i];
	}

	void EnergyElement_SpringHinge::ComputeAndStore_Hessian()
	{
		const Vector3d& x1 = this->m_vDoFs[0]->GetPosition_x();
		const Vector3d& x2 = this->m_vDoFs[1]->GetPosition_x();
		const Vector3d& x3 = this->m_vDoFs[2]->GetPosition_x();

		Real K = (*this->m_pMaterial)[Material::Property::BendingK];
		Real L0 = m_intVolume;
		Real A0 = m_restAngle;

		Real mHx[9][9];

		{
#include "../codegen/SpringHingeHessian.mcg"
		}

		for (int i = 0; i < 9; ++i)
			for (int j = 0; j < 9; ++j)
				m_mHessian(i,j) = mHx[i][j];
	}

}