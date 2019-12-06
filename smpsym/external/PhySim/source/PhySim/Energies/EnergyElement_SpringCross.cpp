//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Energies/EnergyElement_SpringCross.h>

#include <PhySim/Geometry/Node.h>
#include <PhySim/Geometry/Edge.h>
#include <PhySim/Models/Model_Particles.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	EnergyElement_SpringCross::EnergyElement_SpringCross(Model_Particles* pModel, Edge* pEdge0, Edge* pEdge1, Material* pMaterial) : EnergyElement(pModel, pMaterial)
	{
		this->m_pEdge0 = pEdge0;
		this->m_pEdge1 = pEdge1;

		this->m_vDoFs.resize(4);
		this->m_vDoFs[0] = pEdge0->GetOrigin()->DoF();
		this->m_vDoFs[1] = pEdge0->GetHead()->DoF();
		this->m_vDoFs[2] = pEdge1->GetOrigin()->DoF();
		this->m_vDoFs[3] = pEdge1->GetHead()->DoF();

		this->m_vgradient.resize(12);
		this->m_mHessian.resize(12, 12);
		this->m_mHessianPFree.resize(12, 12);
		this->m_mHessianPFull.resize(12, 12);

		this->Init();
	}

	EnergyElement_SpringCross::~EnergyElement_SpringCross()
	{
		// Nothing to do...
	}

	void EnergyElement_SpringCross::Init()
	{
		this->m_restStrain = this->ComputeStrain();
		this->m_intVolume = (this->m_pEdge0->ComputeVolume(Space::MAT) +
						     this->m_pEdge1->ComputeVolume(Space::MAT));
	}

	Real EnergyElement_SpringCross::ComputeStrain(Space s)
	{
		const Vector3d& a0 = this->m_vDoFs[0]->GetPosition(s);
		const Vector3d& b0 = this->m_vDoFs[1]->GetPosition(s);
		const Vector3d& c0 = this->m_vDoFs[2]->GetPosition(s);
		const Vector3d& d0 = this->m_vDoFs[3]->GetPosition(s);
		Real L0 = (b0 - a0).norm();
		Real L1 = (d0 - c0).norm();
		return L0 - L1;
	}

	void EnergyElement_SpringCross::ComputeAndStore_Energy()
	{
		const Vector3d& x1 = this->m_vDoFs[0]->GetPosition_x();
		const Vector3d& x2 = this->m_vDoFs[1]->GetPosition_x();
		const Vector3d& x3 = this->m_vDoFs[2]->GetPosition_x();
		const Vector3d& x4 = this->m_vDoFs[3]->GetPosition_x();

		Real K = (*this->m_pMaterial)[Material::Property::ShearK];
		Real L0 = m_intVolume;
		Real R0 = m_restStrain;

		{
#include "../codegen/SpringCrossEnergy.mcg"

			this->m_energy = t21;
		}
	}

	void EnergyElement_SpringCross::ComputeAndStore_Gradient()
	{
		const Vector3d& x1 = this->m_vDoFs[0]->GetPosition_x();
		const Vector3d& x2 = this->m_vDoFs[1]->GetPosition_x();
		const Vector3d& x3 = this->m_vDoFs[2]->GetPosition_x();
		const Vector3d& x4 = this->m_vDoFs[3]->GetPosition_x();

		Real K = (*this->m_pMaterial)[Material::Property::ShearK];
		Real L0 = m_intVolume;
		Real R0 = m_restStrain;

		Real vgx[12];

		{
#include "../codegen/SpringCrossGradient.mcg"
		}

		for (int i = 0; i < 12; ++i)
			m_vgradient(i) = vgx[i];
	}

	void EnergyElement_SpringCross::ComputeAndStore_Hessian()
	{
		const Vector3d& x1 = this->m_vDoFs[0]->GetPosition_x();
		const Vector3d& x2 = this->m_vDoFs[1]->GetPosition_x();
		const Vector3d& x3 = this->m_vDoFs[2]->GetPosition_x();
		const Vector3d& x4 = this->m_vDoFs[3]->GetPosition_x();

		Real K = (*this->m_pMaterial)[Material::Property::ShearK];
		Real L0 = m_intVolume;
		Real R0 = m_restStrain;

		Real mHx[12][12];

		{
#include "../codegen/SpringCrossHessian.mcg"
		}

		for (int i = 0; i < 12; ++i)
			for (int j = 0; j < 12; ++j)
				m_mHessian(i, j) = mHx[i][j];
	}

}