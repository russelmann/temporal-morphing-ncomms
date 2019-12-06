//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Energies/EnergyElement_SpringLinear_MPoly.h>

#include <PhySim/Geometry/Node.h>
#include <PhySim/Geometry/Edge.h>
#include <PhySim/Models/Model_Particles.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	EnergyElement_SpringLinear_MPoly::EnergyElement_SpringLinear_MPoly(Model_Particles* pModel, Edge* pEdge, Material* pMaterial) : EnergyElement_SpringLinear(pModel, pEdge, pMaterial)
	{
		this->m_pEdge = pEdge;

		this->m_vDoFs.resize(2);
		this->m_vDoFs[0] = pEdge->GetOrigin()->DoF();
		this->m_vDoFs[1] = pEdge->GetHead()->DoF();

		this->m_vgradient.resize(6);
		this->m_mHessian.resize(6, 6);
		this->m_mHessianPFree.resize(6, 6);
		this->m_mHessianPFull.resize(6, 6);

		this->Init();
	}

	EnergyElement_SpringLinear_MPoly::~EnergyElement_SpringLinear_MPoly()
	{
		// Nothing to do...
	}

	void EnergyElement_SpringLinear_MPoly::Init()
	{
		this->m_intVolume = this->m_restLength = this->ComputeRestLength();
	}

	Real EnergyElement_SpringLinear_MPoly::ComputeRestLength()
	{
		return this->m_pEdge->ComputeVolume();
	}

	void EnergyElement_SpringLinear_MPoly::ComputeAndStore_Energy()
	{
		const Vector3d& x1 = this->m_vDoFs[0]->GetPosition_x();
		const Vector3d& x2 = this->m_vDoFs[1]->GetPosition_x();
		const double& L0 = this->m_restLength;
		const double* mC = this->m_mcoeff.data();

		// TODO

		// this->m_energy;
	}

	void EnergyElement_SpringLinear_MPoly::ComputeAndStore_Gradient()
	{
		const Vector3d& x1 = this->m_vDoFs[0]->GetPosition_x();
		const Vector3d& x2 = this->m_vDoFs[1]->GetPosition_x();
		const double& L0 = this->m_restLength;
		const double* vC = this->m_mcoeff.data();

		// TODO

		// this->m_vgradient.block<3, 1>(0, 0) = gb;
	}

	void EnergyElement_SpringLinear_MPoly::ComputeAndStore_Hessian()
	{
		const Vector3d& x1 = this->m_vDoFs[0]->GetPosition_x();
		const Vector3d& x2 = this->m_vDoFs[1]->GetPosition_x();
		const double& L0 = this->m_restLength;
		const double* vC = this->m_mcoeff.data();

		// TODO

		// this->m_mHessian.block<3, 3>(0, 0) = k*(L0inv*ububT + s*((mI - ububT)*Lxinv));
	}

}