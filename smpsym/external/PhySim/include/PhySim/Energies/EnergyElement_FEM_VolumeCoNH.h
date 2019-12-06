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

#include <PhySim/Geometry/Cell.h>

#include <PhySim/Energies/EnergyElement_FEM.h>

#include <PhySim/Energies/Material.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class EnergyElement_FEM_VolumeCoNH : public EnergyElement_FEM
	{

	public:
		EnergyElement_FEM_VolumeCoNH(Model_Particles* pModel, Cell* pCell, Material* pMaterial, int numQ = 1) : EnergyElement_FEM(pModel, pCell, pMaterial, numQ)
		{
			// Nothing to do here...
		}

		virtual ~EnergyElement_FEM_VolumeCoNH(void)
		{
			// Nothing to do here...
		}

		virtual void ComputeEnergyForF(const VectorXd& Fv, Real& U) const
		{
			Real lame1 = (*this->m_pMaterial)[Material::Property::Lame1];
			Real lame2 = (*this->m_pMaterial)[Material::Property::Lame2];

#include "../codegen/FEMVol_CoNH_Energy.mcg"

			U = t31;
		}

		virtual void ComputeGradientForF(const VectorXd& Fv, VectorXd& vg) const
		{
			Real lame1 = (*this->m_pMaterial)[Material::Property::Lame1];
			Real lame2 = (*this->m_pMaterial)[Material::Property::Lame2];

			vg.resize(9);
			Real vgx[9];

#include "../codegen/FEMVol_CoNH_Gradient.mcg"

			for (int i = 0; i < 9; ++i)
				vg(i) = vgx[i]; // Copy
		}

		virtual void ComputeHessianForF(const VectorXd& Fv, MatrixXd& mH) const
		{
			Real lame1 = (*this->m_pMaterial)[Material::Property::Lame1];
			Real lame2 = (*this->m_pMaterial)[Material::Property::Lame2];

			mH.resize(9,9);
			Real mHx[9][9];

#include "../codegen/FEMVol_CoNH_Hessian.mcg"

			for (int i = 0; i < 9; ++i)
			for (int j = 0; j < 9; ++j)
				mH(i,j) = mHx[i][j]; // Copy
		}

	};
}

