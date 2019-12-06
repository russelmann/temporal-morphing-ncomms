//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Energies/EnergyElement_PressureConstant.h>

#include <PhySim/Geometry/Node.h>
#include <PhySim/Geometry/Face_Tri.h>
#include <PhySim/Models/Model_Particles.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	EnergyElement_PressureConstant::EnergyElement_PressureConstant(Model_Particles* pModel, Face_Tri* pFace, Material* pMaterial) : EnergyElement(pModel, pMaterial)
	{
		this->m_pFace = pFace;

		this->m_vDoFs.resize(3);
		this->m_vDoFs[0] = pFace->Nodes()[0]->DoF();
		this->m_vDoFs[1] = pFace->Nodes()[1]->DoF();
		this->m_vDoFs[2] = pFace->Nodes()[2]->DoF();

		this->m_vgradient.resize(9);
		this->m_mHessian.resize(9, 9);
		this->m_mHessianPFree.resize(9, 9);
		this->m_mHessianPFull.resize(9, 9);

		this->Init();
	}

	EnergyElement_PressureConstant::~EnergyElement_PressureConstant()
	{
		// Nothing to do here...
	}

	void EnergyElement_PressureConstant::Init()
	{
		// Nothing to do here...
	}

	void EnergyElement_PressureConstant::ComputeAndStore_Energy()
	{
		Vector3d x0 = m_pFace->Nodes()[0]->DoF()->GetPosition_x();
		Vector3d x1 = m_pFace->Nodes()[1]->DoF()->GetPosition_x();
		Vector3d x2 = m_pFace->Nodes()[2]->DoF()->GetPosition_x();

		double P = -(*this->m_pMaterial)[Material::Property::Pressure];

		// Compute energy the following energy E = -p*V
		// such that the forces, f = -dEdx = p * A * n
		// Volume is computed inside the same way it's
		// computed in Face::computeVolume, but for a
		// single face...

		{
#include "../codegen/PressureConstant_Energy.mcg"

			this->m_energy = t12;
		}
	}

	void EnergyElement_PressureConstant::ComputeAndStore_Gradient()
	{
		Vector3d x0 = m_pFace->Nodes()[0]->DoF()->GetPosition_x();
		Vector3d x1 = m_pFace->Nodes()[1]->DoF()->GetPosition_x();
		Vector3d x2 = m_pFace->Nodes()[2]->DoF()->GetPosition_x();

		double P = -(*this->m_pMaterial)[Material::Property::Pressure];

		{
			Real vgx[9];

#include "../codegen/PressureConstant_Gradient.mcg"

			for (int i = 0; i < 9; ++i)
				m_vgradient(i) = vgx[i];
		}
	}

	void EnergyElement_PressureConstant::ComputeAndStore_Hessian()
	{
		Vector3d x0 = m_pFace->Nodes()[0]->DoF()->GetPosition_x();
		Vector3d x1 = m_pFace->Nodes()[1]->DoF()->GetPosition_x();
		Vector3d x2 = m_pFace->Nodes()[2]->DoF()->GetPosition_x();

		double P = -(*this->m_pMaterial)[Material::Property::Pressure];

		{
			Real mHx[9][9];

#include "../codegen/PressureConstant_Hessian.mcg"

			for (int i = 0; i < 9; ++i)
				for (int j = 0; j < 9; ++j)
					m_mHessian(i, j) = mHx[i][j];
		}
	}

}