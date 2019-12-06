//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#pragma once

#include <PhySim/Solvers/OptimProblem.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class OptimProblem_BasicDynamic : public OptimProblem
	{
	protected:
		IModel* m_pModel;
		Real m_timeStep;

	public:

		OptimProblem_BasicDynamic(IModel* pM, Real dt);
		virtual ~OptimProblem_BasicDynamic();

		virtual void Init(IModel* pM, Real dt);

		virtual void GetVariables(VectorXd& vx) const;
		virtual void SetVariables(const VectorXd& vx);

		virtual bool GetEnergy(Real& e);
		virtual bool GetGradient(VectorXd& vg);
		virtual bool GetHessian(MatrixSd& mH);
	};
}