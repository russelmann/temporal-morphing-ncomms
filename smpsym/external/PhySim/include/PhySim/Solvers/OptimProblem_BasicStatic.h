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

	class OptimProblem_BasicStatic : public OptimProblem
	{
	protected:
		IModel*				m_pM;
		IModel::StateP		m_pS;


	public:

		OptimProblem_BasicStatic(IModel* pM);
		virtual ~OptimProblem_BasicStatic();

		virtual void Init(IModel* pM);

		virtual void GetVariables(VectorXd& vx) const override;
		virtual void SetVariables(const VectorXd& vx) override;

		virtual bool GetEnergy(Real& e) override;
		virtual bool GetGradient(VectorXd& vg) override;
		virtual bool GetHessian(MatrixSd& mH) override;

		virtual bool StartStepCallback() override;
		virtual bool PrePerformStepCallback() override;

		virtual bool IsFullyConstrained() override
		{
			return this->m_pM->BoundaryConditionsLoaded();
		}

	};
}