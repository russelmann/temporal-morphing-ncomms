//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Energies/EnergyElement_SpringLinear_UPoly.h>

#include <PhySim/Geometry/Node.h>
#include <PhySim/Geometry/Edge.h>
#include <PhySim/Models/Model_Particles.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	EnergyElement_SpringLinear_UPoly::EnergyElement_SpringLinear_UPoly(Model_Particles* pModel, Edge* pEdge, Material* pMaterial) : EnergyElement_SpringLinear(pModel, pEdge, pMaterial)
	{
		this->m_pEdge = pEdge;

		this->m_vDoFs.resize(2);
		this->m_vDoFs[0] = pEdge->GetOrigin()->DoF();
		this->m_vDoFs[1] = pEdge->GetHead()->DoF();

		this->m_vgradient.resize(6);
		this->m_mHessian.resize(6, 6);
		this->m_mHessianPFree.resize(6,6);
		this->m_mHessianPFull.resize(6,6);

		this->m_vcoeff.setZero(9);

		this->Init();
	}

	EnergyElement_SpringLinear_UPoly::~EnergyElement_SpringLinear_UPoly()
	{
		// Nothing to do...
	}

	void EnergyElement_SpringLinear_UPoly::ComputeAndStore_Energy()
	{
		const Vector3d& x1 = this->m_vDoFs[0]->GetPosition_x();
		const Vector3d& x2 = this->m_vDoFs[1]->GetPosition_x();
		const double& L0 = this->m_restLength;
		const double* vC = this->m_vcoeff.data();

		{
#include "../codegen/SpringLinearUPolyEnergy.mcg"

			this->m_energy = t25;
		}
	}

	void EnergyElement_SpringLinear_UPoly::ComputeAndStore_Gradient()
	{
		const Vector3d& x1 = this->m_vDoFs[0]->GetPosition_x();
		const Vector3d& x2 = this->m_vDoFs[1]->GetPosition_x();
		const double& L0 = this->m_restLength;
		const double* vC = this->m_vcoeff.data();

		{
			double vgx[6];

#include "../codegen/SpringLinearUPolyGradient.mcg"

			for (int i = 0; i < 6; ++i)
				this->m_vgradient(i) = vgx[i];
		}

		//Vector3d e = x2 - x1;
		//double L = e.norm();
		//Vector3d u = e / L;

		//double s = L0 - L;

		//VectorXd vx(9);
		//for (int i = 0; i < 9; ++i)
		//	vx(i) = pow(s, i);
		//VectorXd vC = m_vcoeff;

		//Real dUds = 0;
		//for (int i = 1; i < 9; ++i)
		//	dUds += vC(i)*i*pow(s, i-1);

		//this->m_vgradient.block(0, 0, 3, 1) = dUds*u;
		//this->m_vgradient.block(3, 0, 3, 1) = -dUds*u;
	}

	void EnergyElement_SpringLinear_UPoly::ComputeAndStore_Hessian()
	{
		const Vector3d& x1 = this->m_vDoFs[0]->GetPosition_x();
		const Vector3d& x2 = this->m_vDoFs[1]->GetPosition_x();
		const double& L0 = this->m_restLength;
		const double* vC = this->m_vcoeff.data();

		{
			double mHx[6][6];

#include "../codegen/SpringLinearUPolyHessian.mcg"

			for (int i = 0; i < 6; ++i)
				for (int j = 0; j < 6; ++j)
					this->m_mHessian(i,j) = mHx[i][j];
		}
	}

}