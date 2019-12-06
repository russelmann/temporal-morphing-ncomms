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

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class Polytope;

	class EnergyElement_FEM: public EnergyElement
	{

	protected:

		Polytope* m_pPoly;
		
		int m_numQ;
		vector<double> m_vQwei;
		vector<VectorXd> m_vQpos;

		MatrixXd m_mrestStrain;

		vector<MatrixXd> m_vmDFDx; // Stores F derivative w.r.t. x for each quadrature point
		vector<MatrixXd> m_vmBH0i; // Stores B*(N0*B)^-1 at rest for each quadrature point

		int m_order;

	public:
		EnergyElement_FEM(Model_Particles* pModel, Polytope* pPoly, Material* pMaterial, int numQ = 1);
		virtual ~EnergyElement_FEM(void);

		virtual void Init();

		virtual void ComputeAndStore_Energy();
		virtual void ComputeAndStore_Gradient();
		virtual void ComputeAndStore_Hessian();

		virtual void ComputeEnergyForF(const VectorXd& vF, Real& U) const = 0;
		virtual void ComputeGradientForF(const VectorXd& vF, VectorXd& vg) const = 0;
		virtual void ComputeHessianForF(const VectorXd& vF, MatrixXd& mH) const = 0;

		virtual MatrixXd ComputeDeformationGradient() const;

		virtual MatrixXd ComputeRestStrain();
		inline virtual MatrixXd GetRestStrain() const { return this->m_mrestStrain; }
		inline virtual void SetRestStrain(MatrixXd rl) { this->m_mrestStrain = rl; }

	};
}

