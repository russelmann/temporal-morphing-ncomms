//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Solvers/LinearSolver.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	LinearSolver::LinearSolver()
	{
		this->m_isInit = false;
	}

	LinearSolver::LinearSolver(const MatrixSd& mA, const LinearSolverOptions& options)
	{
		this->Init(mA, options);
	}

	LinearSolver::~LinearSolver()
	{
		// Nothing to do here...
	}

	void LinearSolver::Init(const MatrixSd& mA, const LinearSolverOptions& options)
	{
		this->Free();

		this->m_options = options;

		this->m_timerSolve = CustomTimer(10, "LIN_SOL", "");

		// Create regularization

		int N = mA.rows();
		m_mR = MatrixSd(N,N);
		VectorTd vT(N);
		for (int i = 0; i < N; ++i)
			if (this->m_options.regSign < 0.0)
				vT[i] = Triplet<Real>(i, i, -1.0);
			else vT[i] = Triplet<Real>(i, i, 1.0);
		m_mR.setFromTriplets(vT.begin(), vT.end());
		m_mR.makeCompressed();

		this->m_isInit = true;
	}

	void LinearSolver::LinearSolver::Free()
	{
		this->m_isInit = true;
	}

	void LinearSolver::GetMatrixData(const MatrixSd& mA, Real& normA, Real& sinT, Real& regT)
	{
		normA = mA.norm();
		sinT = min(1e-6, max(1e-9, 1e-9*normA));
		regT = min(1e-3, max(1e-6, 1e-6*normA));
	}

	SolveResult LinearSolver::Solve(MatrixSd& mA, const VectorXd& vb, VectorXd& vx)
	{
		if (!this->m_isInit)
			this->Init(mA, m_options);

		SolveResult solution;

		if (m_options.profileTime) this->m_timerSolve.resume();

		solution = this->SolveInternal(mA, vb, vx);

		if (m_options.profileTime) this->m_timerSolve.stopStoreLog();

		return solution;
	}

	SolveResult LinearSolver::Solve(MatrixSd& mA, const MatrixXd& mB, MatrixXd& mX)
	{
		if (!this->m_isInit)
			this->Init(mA, m_options);

		SolveResult solution;

		if (m_options.profileTime) this->m_timerSolve.resume();

		solution = this->SolveInternal(mA, mB, mX);

		if (m_options.profileTime) this->m_timerSolve.stopStoreLog();

		return solution;
	}

	SolveResult LinearSolver::Solve(MatrixSd& mA, const MatrixSd& mB, MatrixSd& mX)
	{
		if (!this->m_isInit)
			this->Init(mA, m_options);

		SolveResult solution;

		if (m_options.profileTime) this->m_timerSolve.resume();

		solution = this->SolveInternal(mA, mB, mX);

		if (m_options.profileTime) this->m_timerSolve.stopStoreLog();

		return solution;
	}

}
