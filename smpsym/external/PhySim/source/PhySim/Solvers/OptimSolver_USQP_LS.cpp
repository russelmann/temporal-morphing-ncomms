//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Solvers/OptimSolver_USQP_LS.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	SolveResult OptimSolver_USQP_LS::SolveStep()
	{
		const Options& options = this->GetOptions();

		if (!this->m_isInit)
			return SolveResult::Failure;

		options.m_pProblem->StartStepCallback();

		double E0 = 0;
		double G0 = 0;
		VectorXd dx;
		VectorXd vx;
		VectorXd vg;
		options.m_pProblem->GetVariables(vx);

		options.m_pProblem->GetEnergy(E0);
		options.m_pProblem->GetGradient(vg);
		G0 = vg.norm();

		logSimu("\n--\nStarting problem %s STEP. Energy: %f. Gradient: %f", options.m_pProblem->GetName().c_str(), E0, G0);

		// Check if the gradient is near zero

		if (G0 < options.m_maxError)
			return SolveResult::Success;

		if (options.profileTime) this->m_timerSolve.start();

		// Compute the problem HessianFull

		dx = this->m_dxPrev;
		this->ComputeStep(vx, vg, dx);
		this->m_dxPrev = dx;

		// Limit the size of dx if necessary

		double stepLength = dx.norm();
		if (stepLength > options.m_maxStepSize)
			dx = (dx / stepLength)*options.m_maxStepSize;

		// Line-search in dx. Explanation: the energy G is only a quadratic approximation to E,
		// and hence dx does not necessary minimize E. However, it does for sufficiently small
		// step in the dx direction. Therefore, we perform a line-search in dx until finding
		// a step-length that reduces the potential energy E.

		int itL = 0;
		VectorXd vxL = vx;
		VectorXd dxL = dx;
		double EL = E0;
		bool improved = false;
		for (itL = 0; itL < options.m_lineSearchIters; ++itL)
		{
			// Perform step!

			stepLength = dxL.norm();

			options.m_pProblem->PrePerformStepCallback();
			options.m_pProblem->SetVariables(vx + dxL);
			options.m_pProblem->PosPerformStepCallback();

			// Recompute optimized energy

			options.m_pProblem->GetEnergy(EL);

			// Improvement?

			if (EL < E0)
			{
				improved = true;
				break; // Found!
			}

			// It didn't improve so reduce the step

			dxL = dxL*options.m_lineSearchFactor;
		}

		this->m_dxPrev = dxL;

		if (options.profileTime) this->m_timerSolve.stopStoreLog();

		options.m_pProblem->EndStepCallback();

		// Change improvement

		if (improved)
		{
			logSimu("\n[SUCCESS] Improved. Energy: %f. Bis: %i. Step: %f\n--", EL, itL, stepLength);

			return SolveResult::Success;
		}
		else
		{
			logSimu("\n[FAILURE] Non-Descent. Energy: %f. Bis: %i. Step: %f\n--", EL, itL, stepLength);

			// Roll back the changes!

			options.m_pProblem->SetVariables(vx);

			return SolveResult::NonDesc;
		}
	}

}
