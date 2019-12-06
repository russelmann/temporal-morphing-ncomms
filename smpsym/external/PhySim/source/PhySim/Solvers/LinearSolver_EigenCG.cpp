//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Solvers/LinearSolver_EigenCG.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	LinearSolver_EigenCG::LinearSolver_EigenCG() : LinearSolver()
	{
		// Nothing to do here...
	}

	LinearSolver_EigenCG::LinearSolver_EigenCG(const MatrixSd& mA, const LinearSolverOptions& options) : LinearSolver(mA, options)
	{
		this->Init(mA, options);
	}

	LinearSolver_EigenCG::~LinearSolver_EigenCG()
	{
		// Nothing to do here...
	}

	void LinearSolver_EigenCG::Init(const MatrixSd& mA, const LinearSolverOptions& options)
	{
		LinearSolver::Init(mA, options);

		MatrixSd mAFull = mA.selfadjointView<Lower>();

		this->m_solver.setTolerance(options.maxError);
		this->m_solver.setMaxIterations(options.maxIters);
		this->m_solver.analyzePattern(mAFull);

		if (m_solver.info() != Success)
		{
			logSimu("\n[WARNING] Linear solve: error during preconditioner analysis");
		}
	}

	void LinearSolver_EigenCG::Free()
	{
		// Nothing to do here...
	}

	SolveResult LinearSolver_EigenCG::SolveInternal(MatrixSd& mA, const VectorXd& vb, VectorXd& vx)
	{
		MatrixSd mAFull = mA.selfadjointView<Lower>();

		int N = (int) mAFull.rows(); 
		int M = (int) mAFull.cols();

		Real normA;
		Real sinThres;
		Real regThres;
		this->GetMatrixData(mAFull, normA, sinThres, regThres);

		VectorXd vxPrev = vx;

		// Solve the specified definite positive linear system. 
		// Regularize the matrix if needed to create a DP system.

		int i = 0;
		for (i = 0; i < this->m_options.regIters; ++i)
		{
			if (i != 0)
			{
				logSimu("\n[INFO] Regularizing system, iteration %u", i);

				mAFull += m_mR*regThres*pow(10, i - 1);
			}

			// Pattern analyzed

			m_solver.factorize(mAFull);

			// Check computation

			if (m_solver.info() != Success)
			{
				logSimu("\n[FAILURE] Linear solve: error during preconditioner factorization");
				continue; // Iterate again
			}

			vx = m_solver.solveWithGuess(vb, vxPrev);

			// Check calculation

			if (m_solver.info() != Success)
			{
				logSimu("\n[WARNING] Linear solve: it was impossible to solve accurately");
			}

			//// Check indefiniteness

			//double dot = vb.dot(vx);

			//if (this->m_options.regSign > 0 && dot < 0.0)
			//{
			//	logSimu("\n[FAILURE] Linear solve: indefinite matrix, dot: %.9f", dot);
			//	continue;
			//}
			//if (this->m_options.regSign < 0 && dot > 0.0)
			//{
			//	logSimu("[FAILURE] Linear solve: indefinite matrix, dot: %.9f", dot);
			//	continue;
			//}

			//// Check exact solution

			//VectorXd vbTest = mA.selfadjointView<Lower>()*vx;
			//double absError = (vb - vbTest).norm();
			//double relError = absError / vb.norm();
			//if (absError > this->m_options.maxError && relError > this->m_options.maxError)
			//{
			//	logSimu("\n[FAILURE] Linear solve: inexact solution, error: %.9f", relError);
			//	continue; // Iterate again
			//}

			logSimu("\n[SUCCESS] Linear solve: solved using Eigen CG. Reg: %d, Iter: %d, Error: %f", i, m_solver.iterations(), m_solver.error());

			return SolveResult::Success;
		}

		return SolveResult::Failure;
	}

	SolveResult LinearSolver_EigenCG::SolveInternal(MatrixSd& mA, const MatrixXd& mB, MatrixXd& mX)
	{
		assert(false);

		return SolveResult::Failure;
	}

	SolveResult LinearSolver_EigenCG::SolveInternal(MatrixSd& mA, const MatrixSd& mB, MatrixSd& mX)
	{
      assert(false);

		return SolveResult::Failure;
	}


	//	// Solve the specified definite positive linear system. 
	//	// Regularize the matrix if needed to create a DP system.

	//	int i = 0;
	//	for (i = 0; i < this->options.regIters; ++i)
	//	{
	//		if (i != 0)
	//		{
	//			logSimu("\n[INFO] Regularizing system, iteration %u", i);

	//			mA += mR*regThres*pow(10, i - 1);
	//		}

	//		ConjugateGradient<MatrixSd> solver;
	//		solver.setTolerance(this->options.maxError);
	//		solver.setMaxIterations(this->options.maxIters);

	//		solver.compute(mA);

	//		// Check computation

	//		if (solver.info() != Success)
	//		{
	//			logSimu("\n[FAILURE] Linear solve: error during computation");
	//			continue; // Iterate again
	//		}

	//		vx = solver.solve(vb);

	//		// Check calculation

	//		if (solver.info() != Success)
	//		{
	//			logSimu("\n[FAILURE] Linear solve: error during calculation");
	//			continue; // Iterate again
	//		}

	//		// Check indefiniteness

	//		double dot = vb.dot(vx);

	//		if (this->options.regSign > 0 && dot < 0.0)
	//		{
	//			logSimu("\n[FAILURE] Linear solve: indefinite matrix, dot: %.9f", dot);
	//			continue;
	//		}
	//		if (this->options.regSign < 0 && dot > 0.0)
	//		{
	//			logSimu("\n[FAILURE] Linear solve: indefinite matrix, dot: %.9f", dot);
	//			continue;
	//		}

	//		// Check exact solution

	//		VectorXd vbTest = mA.selfadjointView<Lower>()*vx;
	//		double absError = (vb - vbTest).norm();
	//		double relError = absError / vb.norm();
	//		if (absError > this->options.maxError && relError > this->options.maxError)
	//		{
	//			logSimu("\n[FAILURE] Linear solve: inexact solution, error: %.9f", relError);
	//			continue; // Iterate again
	//		}

	//		logSimu("\n[SUCCESS] Linear solve: solved using Conjugate GradientFull after %i regularization steps", i);

	//		return SolveResult::Success;
	//	}

	//	return SolveResult::Failure;
	//}

	//SolveResult LinearSolver_EigenCG::Solve_LU(MatrixSd& mA, const VectorXd& vb, VectorXd& vx)
	//{
	//	int N;
	//	Real normA;
	//	MatrixSd mR;
	//	Real sinThres;
	//	Real regThres;
	//	this->PrepareLinearSolve(mA, N, normA, sinThres, regThres, mR);

	//	// Solve the specified definite positive linear system. 
	//	// Regularize the matrix if needed to create a DP system.

	//	int i = 0;
	//	for (i = 0; i < this->options.regIters; ++i)
	//	{
	//		if (i != 0)
	//		{
	//			logSimu("\n[INFO] Regularizing system, iteration %u", i);

	//			mA += mR*regThres*pow(10, i - 1);
	//		}

	//		SparseLU<MatrixSd> solver;

	//		solver.compute(mA);

	//		// Check factorization

	//		if (solver.info() != Success)
	//		{
	//			logSimu("\n[FAILURE] Linear solve: error during factorization");
	//			continue; // Iterate again
	//		}

	//		vx = solver.solve(vb);

	//		// Check calculation

	//		if (solver.info() != Success)
	//		{
	//			logSimu("\n[FAILURE] Linear solve: error during calculation");
	//			continue; // Iterate again
	//		}

	//		// Check indefiniteness

	//		double dot = vb.dot(vx);

	//		if (this->options.regSign > 0 && dot < 0.0)
	//		{
	//			logSimu("\n[FAILURE] Linear solve: indefinite matrix, dot: %.9f", dot);
	//			continue;
	//		}
	//		if (this->options.regSign < 0 && dot > 0.0)
	//		{
	//			logSimu("\n[FAILURE] Linear solve: indefinite matrix, dot: %.9f", dot);
	//			continue;
	//		}

	//		// Check exact solution

	//		VectorXd vbTest = mA.selfadjointView<Lower>()*vx;
	//		double absError = (vb - vbTest).norm();
	//		double relError = absError / vb.norm();
	//		if (absError > this->options.maxError && relError > this->options.maxError)
	//		{
	//			logSimu("\n[FAILURE] Linear solve: inexact solution, error: %.9f", relError);
	//			continue; // Iterate again
	//		}

	//		logSimu("\n[SUCCESS] Linear solve: solved using LU Factorization after %i regularization steps", i);

	//		return SolveResult::Success;
	//	}

	//	return SolveResult::Failure;
	//}

}
