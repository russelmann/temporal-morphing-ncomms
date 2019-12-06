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

#include <PhySim/Utils/CustomTimer.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class LinearSolver: public ILinearSolver
	{
	protected:

		LinearSolverOptions m_options;

		CustomTimer m_timerSolve;

		bool m_isInit;

		MatrixSd m_mR;

	public:

		LinearSolver();
		LinearSolver(const MatrixSd& mA, const LinearSolverOptions& options);
		virtual void Init(const MatrixSd& mA, const LinearSolverOptions& options);
		virtual ~LinearSolver();

		virtual LinearSolverOptions& GetOptions() { return this->m_options; }
		virtual void SetOptions(LinearSolverOptions& op) { m_options = op; }

		virtual void GetMatrixData(const MatrixSd& mA, Real& normA, Real& sinT, Real& regT);

		virtual SolveResult Solve(MatrixSd& mA, const VectorXd& vb, VectorXd& vx);
		virtual SolveResult Solve(MatrixSd& mA, const MatrixXd& mB, MatrixXd& mX);
		virtual SolveResult Solve(MatrixSd& mA, const MatrixSd& mB, MatrixSd& mX);

	protected:

		virtual SolveResult SolveInternal(MatrixSd& mA, const VectorXd& vb, VectorXd& vx) = 0;
		virtual SolveResult SolveInternal(MatrixSd& mA, const MatrixXd& mB, MatrixXd& mX) = 0;
		virtual SolveResult SolveInternal(MatrixSd& mA, const MatrixSd& mB, MatrixSd& mX) = 0;

		virtual void Free();

	};
}
