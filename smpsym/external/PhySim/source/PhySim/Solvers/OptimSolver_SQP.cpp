//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Solvers/OptimSolver_SQP.h>

#include <PhySim/Solvers/LinearSolver_EigenLDLT.h>
#include <PhySim/Solvers/LinearSolver_EigenCG.h>
#include <PhySim/Solvers/LinearSolver_BiCGSTAB.h>
#include <PhySim/Solvers/LinearSolver_CholmodLDLT.h>
#include <PhySim/Solvers/LinearSolver_SSparseSPQR.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	OptimSolver_SQP::OptimSolver_SQP()
	{
		this->m_pOptions.reset(new Options());

		this->m_pLSolver = NULL;
	}

	OptimSolver_SQP::~OptimSolver_SQP()
	{
		if (this->m_pLSolver != NULL)
			delete this->m_pLSolver;
	}

	void OptimSolver_SQP::Init()
	{
		OptimSolver::Init();

		// Initialize linear solver

		const Options& options = GetOptions();

		LinearSolverOptions linearOptions;
		linearOptions.type = options.m_lsSolverType;

		MatrixSd mA;
		options.m_pProblem->GetHessian(mA);

		if (this->m_pLSolver != NULL)
			delete this->m_pLSolver;

		switch (linearOptions.type)
		{
			case LSSolverType::LS_EigenLDLT: 
				this->m_pLSolver = new LinearSolver_EigenLDLT(mA, linearOptions);
				break;
			case LSSolverType::LS_EigenCG:
				this->m_pLSolver = new LinearSolver_EigenCG(mA, linearOptions);
				break;
			case LSSolverType::LS_BiCGSTAB:
				this->m_pLSolver = new LinearSolver_BiCGSTAB(mA, linearOptions);
				break;
#ifdef USE_SSPARSE
			case LSSolverType::LS_CholmodLDLT:
				this->m_pLSolver = new LinearSolver_CholmodLDLT(mA, linearOptions);
				break;
			case LSSolverType::LS_SSparseSPQR:
				this->m_pLSolver = new LinearSolver_SSparseSPQR(mA, linearOptions);
				break;
#endif
#ifdef USE_CUDA
			case LSSolverType::LS_CUDASC:
				this->m_pLSolver = new LinearSolver_CudaSC(mA, linearOptions);
				break;
#endif
		}

		// Initialize previous step

		this->m_dxPrev = VectorXd::Zero(options.m_pProblem->GetNumVariables());

		// Initialize BFGS stuff

		this->m_isBFGSInit = false;
	}

	bool OptimSolver_SQP::ComputeStep(const VectorXd& vx, const VectorXd& vg, VectorXd& dx)
	{
		//return this->ComputeStep_NewtonAligned(vx, vg, dx);

		switch (this->GetOptions().m_qpSolverType)
		{
		case QPSolverType::QP_Newton: return this->ComputeStep_NewtonDefault(vx, vg, dx); break;
		case QPSolverType::QP_BFGS_D: return this->ComputeStep_BFGSDirect(vx, vg, dx); break;
		case QPSolverType::QP_BFGS_I: return this->ComputeStep_BFGSInverse(vx, vg, dx); break;
		case QPSolverType::QP_LBFGS: return this->ComputeStep_LBFGS(vx, vg, dx); break;
		case QPSolverType::QP_Steepest: dx = -vg; return true;
		default:
         assert(false); ;
		}

		return false;
	}

	bool OptimSolver_SQP::ComputeStep_NewtonScaled(const VectorXd& vx, const VectorXd& vg, VectorXd& dx)
	{
		return false;
	}

	bool OptimSolver_SQP::ComputeStep_NewtonAligned(const VectorXd& vx, const VectorXd& vg, VectorXd& dx)
	{
		MatrixSd mH;

		const Options& options = this->GetOptions();

		options.m_pProblem->GetHessian(mH);

		MatrixSd mHs = mH.selfadjointView<Lower>();
		MatrixXd mHd = mHs.toDense();
		JacobiSVD<MatrixXd> svdSolver(mHd, ComputeFullV);
		const VectorXd& vs = svdSolver.singularValues();
		const MatrixXd& mV = svdSolver.matrixV();

		options.m_pProblem->PreComputeStepCallback();
		
		dx = -mV.transpose()*vg;
	
		for (int i = 0; i < dx.size(); ++i)
			dx(i) /= vs(i);

		dx = mV*dx;

		options.m_pProblem->PosComputeStepCallback();

		return true;
	}

	bool OptimSolver_SQP::ComputeStep_NewtonDefault(const VectorXd& vx, const VectorXd& vg, VectorXd& dx)
	{
		const Options& options = this->GetOptions();

		MatrixSd mH;

		options.m_pProblem->GetHessian(mH);

		// Solve the system H*dx = -g

		dx = this->m_dxPrev;
		options.m_pProblem->PreComputeStepCallback();
		this->m_pLSolver->Solve(mH, -(1.0)*vg, dx);
		options.m_pProblem->PosComputeStepCallback();
		this->m_dxPrev = dx;

		return true;
	}

	bool OptimSolver_SQP::ComputeStep_BFGSDirect(const VectorXd& vx, const VectorXd& vg, VectorXd& dx)
	{
		const Options& options = this->GetOptions();

		this->UpdateBFGS_Direct(vx, vg);

		// Solve the system H*dx = -g

		dx = this->m_dxPrev;
		MatrixSd& mH = this->m_mBFGS.m_msparseMatrix;
		options.m_pProblem->PreComputeStepCallback();
		this->m_pLSolver->Solve(mH, -(1.0)*vg, dx);
		options.m_pProblem->PosComputeStepCallback();
		this->m_dxPrev = dx;

		return true;
	}

	bool OptimSolver_SQP::ComputeStep_BFGSInverse(const VectorXd& vx, const VectorXd& vg, VectorXd& dx)
	{
		const Options& options = this->GetOptions();

		this->UpdateBFGS_Direct(vx, vg);

		// Solve multiplication -H^-1*g

		dx = this->m_dxPrev;
		options.m_pProblem->PreComputeStepCallback();
		dx = -this->m_mBFGS.m_msparseMatrix*vg;
		options.m_pProblem->PosComputeStepCallback();
		this->m_dxPrev = dx;
		
		return true;
	}

	bool OptimSolver_SQP::ComputeStep_LBFGS(const VectorXd& vx, const VectorXd& vg, VectorXd& vdx)
	{
		if (!this->m_isBFGSInit)
		{
			this->m_vx_P = vx;
			this->m_vg_P = vg;

			m_M = 10;
			m_N = vx.size();
			m_mS.resize(m_N, m_M);
			m_mY.resize(m_N, m_M);
			m_vys.resize(m_M);
			m_va.resize(m_M);
			m_K = 0;

			// Steep!
			vdx = -vg;

			this->m_isBFGSInit = true;

			logSimu("\n[INFO] Initializing BFGS step to gradient");

			return true;
		}

		// Update s and y
		// s_{k+1} = x_{k+1} - x_k
		// y_{k+1} = g_{k+1} - g_k

		MapVec svec(&m_mS(0, m_K), m_N);
		MapVec yvec(&m_mY(0, m_K), m_N);
		svec.noalias() = vx - m_vx_P;
		yvec.noalias() = vg - m_vg_P;

		// ys = y's = 1/rho
		// yy = y'y

		Real ys = yvec.dot(svec);
		Real yy = yvec.squaredNorm();
		m_vys[m_K] = ys;

		logSimu("\n[INFO] Updated BFGS step with discriminant %f", ys);

		vdx.noalias() = -vg;

		// Recursive formula to compute d = -H * g

		int bound = std::min(m_M, m_K+1);
		m_K = (m_K + 1) % m_M;
		int j = m_K;
		for (int i = 0; i < bound; i++)
		{
			j = (j + m_M - 1) % m_M;
			MapVec sj(&m_mS(0, j), m_N);
			MapVec yj(&m_mY(0, j), m_N);
			m_va[j] = sj.dot(vdx) / m_vys[j];
			vdx.noalias() -= m_va[j] * yj;
		}

		vdx *= (ys / yy);

		for (int i = 0; i < bound; i++)
		{
			MapVec sj(&m_mS(0, j), m_N);
			MapVec yj(&m_mY(0, j), m_N);
			Real beta = yj.dot(vdx) / m_vys[j];
			vdx.noalias() += (m_va[j] - beta) * sj;
			j = (j + 1) % m_N;
		}

		return true;
	}

	bool OptimSolver_SQP::UpdateBFGS_Inverse(const VectorXd& vxNew, const VectorXd& vgNew)
	{
		if (!this->m_isBFGSInit) return this->InitialBFGS(vxNew, vgNew);

		// TODO

		return false;
	}

	bool OptimSolver_SQP::UpdateBFGS_Direct(const VectorXd& vxNew, const VectorXd& vgNew)
	{
		if (!this->m_isBFGSInit) return this->InitialBFGS(vxNew, vgNew);

		m_vs = vxNew - this->m_vx_P;
		m_vy = vgNew - this->m_vg_P;
		
		this->m_vx_P = vxNew;
		this->m_vg_P = vgNew;
		
		double D = m_vy.dot(m_vs);
		
		if (D > 1e-6) // ?
		{
			logSimu("\n[SUCCESS] Updating BFGS matrix with discriminant %f", D);

			Real rho = 1.0 / D; // For too small discriminant, the update is invalid

			// Update the HessianFull matrix (direct one) (Numerical Optimization, Nocedal, Wright, page 25)
			// Hk+1 = (I - rho * m_sk * m_yk^T) * Hk * (I - rho * m_yk * m_sk^T) + rho * m_yk * m_yk^T

			MatrixSd& W = this->m_mBFGS.m_msparseMatrix;
			int N = this->m_mBFGS.m_msparseMatrix.rows();

			m_mBFGSMod_A(N, N);
			m_mBFGSMod_B(N, N);

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j <= i; j++)
				{
					if (i == j)
					{
						m_mBFGSMod_A(i, i) = 1.0 - rho * m_vy[i] * m_vs[i];
					}
					else
					{
						m_mBFGSMod_A(i, j) = -rho * m_vs[i] * m_vy[j];
						m_mBFGSMod_A(j, i) = -rho * m_vs[j] * m_vy[i];
					}
			
					m_mBFGSMod_B(i, j) = rho * m_vy[i] * m_vy[j];
					m_mBFGSMod_B(j, i) = rho * m_vy[j] * m_vy[i];
				}
			}
			
			W = (m_mBFGSMod_A*this->m_mBFGS.m_msparseMatrix.toDense()*m_mBFGSMod_A.transpose() + m_mBFGSMod_B).sparseView();

			//size_t numCoef = m_mBFGS.m_vpointTriplets.size();

			//for (size_t c0 = 0; c0 < numCoef; ++c0)
			//{
			//	int i = this->m_mBFGS.m_vpointTriplets[c0].row();
			//	int j = this->m_mBFGS.m_vpointTriplets[c0].col();
			//	
			//	Real value = rho * m_yk[i] * m_yk[j];

			//	for (size_t c1 = 0; c1 < numCoef; ++c1)
			//	{
			//		int p = this->m_mBFGS.m_vpointTriplets[c1].row();
			//		int q = this->m_mBFGS.m_vpointTriplets[c1].col();

			//		Real Apq = this->m_mBFGS.m_vvalueTriplets[c1].value();

			//		Real Bip = - rho * m_sk[i] * m_yk[p];
			//		Real Bjq = - rho * m_sk[j] * m_yk[q];
			//		if (i == p) Bip += 1.0;
			//		if (j == q) Bjq += 1.0;

			//		value += Bip*Apq*Bjq;

			//		//if (p != q)
			//		//{
			//		//	int p = this->m_mBFGS.m_vpointTriplets[c1].col();
			//		//	int q = this->m_mBFGS.m_vpointTriplets[c1].row();

			//		//	Real Apq = this->m_mBFGS.m_vvalueTriplets[c1].value();

			//		//	Real Bip = -rho * m_sk[i] * m_yk[p];
			//		//	Real Bjq = -rho * m_sk[j] * m_yk[q];
			//		//	if (i == p) Bip += 1.0;
			//		//	if (j == q) Bjq += 1.0;

			//		//	value += Bip*Apq*Bjq;
			//		//}
			//	}

			//	*this->m_mBFGS.m_vpointTriplets[c0].value() = value;
			//}
		
			// Store current value
		
			//for (size_t c0 = 0; c0 < numCoef; ++c0)
			//{
			//	this->m_mBFGS.m_vvalueTriplets[c0] =
			//		Triplet<Real>(
			//			this->m_mBFGS.m_vvalueTriplets[c0].row(),
			//			this->m_mBFGS.m_vvalueTriplets[c0].col(),
			//			*this->m_mBFGS.m_vpointTriplets[c0].value());
			//}
		}
		else
		{
			logSimu("\n[FAILURE] BFGS discriminant to low %f", D);

			this->InitialBFGS(vxNew, vgNew);
		}

		return true;
	}

	bool OptimSolver_SQP::InitialBFGS(const VectorXd& vxNew, const VectorXd& vgNew)
	{
		const Options& options = this->GetOptions();

		logSimu("\n[INFO] Initializing BFGS matrix to Hessian");

		this->m_vx_P = vxNew;
		this->m_vg_P = vgNew;

		options.m_pProblem->GetHessian(this->m_mBFGS.m_msparseMatrix);

		MatrixSd mHfull = this->m_mBFGS.m_msparseMatrix.selfadjointView<Lower>();
		this->m_mBFGS.m_msparseMatrix = mHfull; // Copy ncessary to do it effective

		eigenSparseMatrixToPointers(this->m_mBFGS.m_msparseMatrix, this->m_mBFGS.m_vpointTriplets);
		eigenSparseMatrixToTriplets(this->m_mBFGS.m_msparseMatrix, this->m_mBFGS.m_vvalueTriplets);

		// Set identity

		this->m_mBFGS.m_msparseMatrix *= 0.0;
		int N = m_mBFGS.m_msparseMatrix.rows();
		for (int i = 0; i < N; ++i)
			this->m_mBFGS.m_msparseMatrix.coeffRef(i, i) = 1.0;

		this->m_isBFGSInit = true;

		return true;
	}


}
