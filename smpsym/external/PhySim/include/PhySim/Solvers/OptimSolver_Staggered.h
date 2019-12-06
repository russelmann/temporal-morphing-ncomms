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

#include <PhySim/Solvers/OptimSolver.h>

#include <PhySim/Solvers/OptimProblem_Subproblem.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;


	class OptimSolver_Staggered : public OptimSolver
	{
	public:
		struct Options : public OptimSolver::Options
		{
			vector<shared_ptr<IOptimProblem>>		m_vpProblem;
			vector<shared_ptr<IOptimSolver>>		m_vpSolvers;

			Options()
			{
				
			}
		};

	protected:
		int										m_numProblem;
		int										m_curProblem;

	public:

		OptimSolver_Staggered();
		virtual ~OptimSolver_Staggered();
		virtual void Init();

		virtual SolveResult SolveStep();

		virtual Options& GetOptions() override { return *static_cast<Options*>(this->m_pOptions.get()); }

	};
}
