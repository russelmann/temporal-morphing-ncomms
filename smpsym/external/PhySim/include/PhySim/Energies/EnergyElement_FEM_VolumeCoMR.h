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

	class EnergyElement_FEM_VolumeCoMR : public EnergyElement_FEM
	{

	public:
		EnergyElement_FEM_VolumeCoMR(Model_Particles* pModel, Cell* pCell, Material* pMaterial, int numQ = 1) : EnergyElement_FEM(pModel, pCell, pMaterial, numQ)
		{
			// Nothing to do here...
		}

		virtual ~EnergyElement_FEM_VolumeCoMR(void)
		{
			// Nothing to do here...
		}

		virtual void ComputeEnergyForF(const VectorXd& Fv, Real& U) const
		{
			Real bulk = (*this->m_pMaterial)[Material::Property::Bulk];
			Real mr01 = (*this->m_pMaterial)[Material::Property::Mooney01];
			Real mr10 = (*this->m_pMaterial)[Material::Property::Mooney10];

#include "../codegen/FEMVol_CoMR_Energy.mcg"

			U = t61;
		}

		virtual void ComputeGradientForF(const VectorXd& Fv, VectorXd& vg) const
		{
			Real bulk = (*this->m_pMaterial)[Material::Property::Bulk];
			Real mr01 = (*this->m_pMaterial)[Material::Property::Mooney01];
			Real mr10 = (*this->m_pMaterial)[Material::Property::Mooney10];

			vg.resize(9);
			Real vgx[9];

#include "../codegen/FEMVol_CoMR_Gradient.mcg"

			for (int i = 0; i < 9; ++i)
				vg(i) = vgx[i]; // Copy
		}

		virtual void ComputeHessianForF(const VectorXd& Fv, MatrixXd& mH) const
		{
			Real bulk = (*this->m_pMaterial)[Material::Property::Bulk];
			Real mr01 = (*this->m_pMaterial)[Material::Property::Mooney01];
			Real mr10 = (*this->m_pMaterial)[Material::Property::Mooney10];

			mH.resize(9, 9);
			Real mHx[9][9];

#include "../codegen/FEMVol_CoMR_Hessian.mcg"

			for (int i = 0; i < 9; ++i)
				for (int j = 0; j < 9; ++j)
					mH(i, j) = mHx[i][j]; // Copy
		}

	};
}

