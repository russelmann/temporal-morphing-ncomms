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

#include <PhySim/Solvers/OptimSolver_SQP.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class OptimSolver_USQP_LS : public OptimSolver_SQP
	{

	public:

		OptimSolver_USQP_LS()
		{
			// Nothing to do here...
		}

		SolveResult SolveStep();

	};
}
