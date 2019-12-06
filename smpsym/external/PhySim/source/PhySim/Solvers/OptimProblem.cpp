//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Solvers/OptimProblem.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	void OptimProblem::GetUpperBound(VectorXd& vub) const
	{
		vub.setConstant(this->m_N, HUGE_VAL);
	}

	void OptimProblem::GetLowerBound(VectorXd& vlb) const
	{
		vlb.setConstant(this->m_N, -HUGE_VAL);
	}

	bool OptimProblem::GetEnergy(Real& e)
	{
		e = HUGE_VAL;

		return false;
	}

	bool OptimProblem::GetGradient(VectorXd& vg)
	{
		vg.setConstant(m_N, 0);

		return false;
	}

	bool OptimProblem::GetHessian(MatrixSd& mH)
	{
		mH = MatrixSd(m_N, m_N);
		
		return false;
	}

	bool OptimProblem::GetConstraints(VectorXd& vc)
	{
		vc.setConstant(m_N, 0);

		return false;
	}

	bool OptimProblem::GetJacobian(MatrixSd& mJ)
	{
		mJ = MatrixSd(m_M, m_N);

		return false;
	}

}
