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

#include <PhySim/Solvers/LinearSolver.h>

#ifdef USE_SSPARSE

#include <Eigen/CholmodSupport>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class LinearSolver_CholmodLDLT : public LinearSolver
	{
		CholmodSimplicialLDLT<MatrixSd> m_solver;

	public:

		LinearSolver_CholmodLDLT();
		LinearSolver_CholmodLDLT(const MatrixSd& mA, const LinearSolverOptions& options);
		virtual void Init(const MatrixSd& mA, const LinearSolverOptions& options);
		virtual ~LinearSolver_CholmodLDLT();

	protected:

		virtual SolveResult SolveInternal(MatrixSd& mA, const VectorXd& vb, VectorXd& vx);
		virtual SolveResult SolveInternal(MatrixSd& mA, const MatrixXd& mB, MatrixXd& mX);
		virtual SolveResult SolveInternal(MatrixSd& mA, const MatrixSd& mB, MatrixSd& mX);

		virtual void Free();

	};
}

#endif