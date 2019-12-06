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

#include <PhySim/Geometry/Face.h>

#include <PhySim/Energies/EnergyElement_FEM.h>

#include <PhySim/Energies/Material.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class EnergyElement_FEM_MembraneInNH : public EnergyElement_FEM
	{

	public:
		EnergyElement_FEM_MembraneInNH(Model_Particles* pModel, Face* pFace, Material* pMaterial, int numQ = 1) : EnergyElement_FEM(pModel, pFace, pMaterial, numQ)
		{
			// Nothing to do here...
		}

		virtual ~EnergyElement_FEM_MembraneInNH(void)
		{
			// Nothing to do here...
		}

		virtual void ComputeEnergyForF(const VectorXd& Fv, Real& U) const
		{
			Real lame1 = (*this->m_pMaterial)[Material::Property::Lame1];
			Real lame2 = (*this->m_pMaterial)[Material::Property::Lame2];

#include "../codegen/FEMMem_InNH_Energy.mcg"

			U = t26;
		}

		virtual void ComputeGradientForF(const VectorXd& Fv, VectorXd& vg) const
		{
			Real lame1 = (*this->m_pMaterial)[Material::Property::Lame1];
			Real lame2 = (*this->m_pMaterial)[Material::Property::Lame2];

			vg.resize(6);
			Real vgx[6];

#include "../codegen/FEMMem_InNH_Gradient.mcg"

			for (int i = 0; i < 6; ++i)
				vg(i) = vgx[i]; // Copy
		}

		virtual void ComputeHessianForF(const VectorXd& Fv, MatrixXd& mH) const
		{
			Real lame1 = (*this->m_pMaterial)[Material::Property::Lame1];
			Real lame2 = (*this->m_pMaterial)[Material::Property::Lame2];

			mH.resize(6, 6);
			Real mHx[6][6];

#include "../codegen/FEMMem_InNH_Hessian.mcg"

			for (int i = 0; i < 6; ++i)
				for (int j = 0; j < 6; ++j)
					mH(i, j) = mHx[i][j]; // Copy
		}

	};
}

