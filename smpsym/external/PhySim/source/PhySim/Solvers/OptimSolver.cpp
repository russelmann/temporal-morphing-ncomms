//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Solvers/OptimSolver.h>

#include <PhySim/Solvers/LinearSolver_EigenLDLT.h>
#include <PhySim/Solvers/LinearSolver_EigenCG.h>
#include <PhySim/Solvers/LinearSolver_BiCGSTAB.h>
#include <PhySim/Solvers/LinearSolver_CholmodLDLT.h>
#include <PhySim/Solvers/LinearSolver_SSparseSPQR.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	OptimSolver::OptimSolver()
	{
		this->m_isInit = false;
		
		this->m_pOptions.reset(new Options());
	}

	OptimSolver::~OptimSolver()
	{
		// Nothing to do here...
	}

	void OptimSolver::Init()
	{
		// Initialize custom timer

		this->m_timerSolve = CustomTimer(10, "OPT_SOL", "");

		// Initialize step/solve calls

		this->m_fullCallback = NULL;
		this->m_stepCallback = NULL;

		this->m_isInit = true;
	}

	OptimSolver::Options& OptimSolver::GetOptions()
	{
		return *this->m_pOptions.get();
	}

	Real OptimSolver::ComputeOptimality()
	{
		VectorXd vg;
		this->GetOptions().m_pProblem->GetGradient(vg);
		return vg.norm();
	}

	Real OptimSolver::ComputeFeasibility()
	{
		// TODO
		return 0;
	}

	bool OptimSolver::IsOptimal()
	{
		return this->ComputeOptimality() < GetOptions().m_maxError;
	}

	bool OptimSolver::IsFeasible()
	{
		// TODO
		return false;
	}

	SolveResult OptimSolver::SolveFull()
	{
		if (!this->m_isInit)
			return SolveResult::Failure;

		this->m_forceSolveStop = false;

		this->m_pOptions->m_pProblem->StartSolveCallback();

		VectorXd vg;
		Real vgNorm = 0;
		int itCount = 1;

		bool valid = this->m_pOptions->m_pProblem->GetGradient(vg);

		logSimu("\n\n");
		logSimu("\n-------------------------------------------------------------------------------------------");
		logSimu("\n-------------------------------------------------------------------------------------------");
		logSimu("\nStarting problem %s SOLVE. Gradient: %f", this->m_pOptions->m_pProblem->GetName().c_str(), vg.norm());
		logSimu("\n-------------------------------------------------------------------------------------------");
		logSimu("\n-------------------------------------------------------------------------------------------");

		SolveResult fullResult = SolveResult::Success;

		// While the gradient is not zero, iterate

		while (true)
		{
			if (this->m_forceSolveStop)
			{
				fullResult = SolveResult::Stopped;
				break; // User requested solve stop
			}

			// Solve step

			SolveResult stepResult = this->SolveStep();

			this->OnEvent_AfterStep();

			if (stepResult == SolveResult::NonDesc)
			{
				fullResult = stepResult;
				break; // Cannnot improve
			}

			if (!this->m_pOptions->m_pProblem->IsFullyConstrained())
				continue; // If not fully loaded, continue

			// Compute forces

			bool validGrad = this->m_pOptions->m_pProblem->GetGradient(vg);
			if (validGrad && vg.norm() < this->m_pOptions->m_maxError)
			{
				fullResult = SolveResult::Success;
				break; // Already found the solution
			}

			// Increment step

			itCount++;

			if (itCount > this->m_pOptions->m_maxIters)
			{
				fullResult = SolveResult::MaxIter;
				break; // Reached max iterations
			}
		}

		this->m_pOptions->m_pProblem->EndSolveCallback();

		this->OnEvent_AfterFull();

		this->m_pOptions->m_pProblem->GetGradient(vg);

		logSimu("\n-------------------------------------------------------------------------------------------");
		logSimu("\n-------------------------------------------------------------------------------------------");
		logSimu("\nFinishing problem. Steps: %d, Gradient: %f", itCount, vg.norm());
		logSimu("\n-------------------------------------------------------------------------------------------");
		logSimu("\n-------------------------------------------------------------------------------------------");
		logSimu("\n\n");

		return fullResult;
	}

	SolveResult OptimSolver::SolveStep()
	{
		// To be inherit and implemented

		return SolveResult::Failure;
	}

}
