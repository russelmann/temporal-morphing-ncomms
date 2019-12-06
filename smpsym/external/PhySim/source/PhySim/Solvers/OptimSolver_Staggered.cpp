//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Solvers/OptimSolver_Staggered.h>

#include <PhySim/Solvers/LinearSolver_EigenLDLT.h>
#include <PhySim/Solvers/LinearSolver_EigenCG.h>
#include <PhySim/Solvers/LinearSolver_BiCGSTAB.h>
#include <PhySim/Solvers/LinearSolver_CholmodLDLT.h>
#include <PhySim/Solvers/LinearSolver_SSparseSPQR.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	OptimSolver_Staggered::OptimSolver_Staggered()
	{
		this->m_pOptions.reset(new Options());
	}

	OptimSolver_Staggered::~OptimSolver_Staggered()
	{
		// Nothing to do here...
	}

	void OptimSolver_Staggered::Init()
	{
		OptimSolver::Init();

		const Options& options = GetOptions();

		this->m_numProblem = (int)options.m_vpSolvers.size();
		this->m_curProblem = 0;
	}


	SolveResult OptimSolver_Staggered::SolveStep()
	{
		const Options& options = GetOptions();

		if (options.m_vpSolvers[m_curProblem]->IsOptimal())
			m_curProblem = (m_curProblem + 1) % m_numProblem;

		return options.m_vpSolvers[m_curProblem]->SolveStep();
	}

}
