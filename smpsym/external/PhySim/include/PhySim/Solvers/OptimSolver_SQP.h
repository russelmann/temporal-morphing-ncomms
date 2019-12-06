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

namespace PhySim
{
	using namespace std;
	using namespace Eigen;


	class OptimSolver_SQP : public OptimSolver
	{
	public:
		struct Options : public OptimSolver::Options
		{
			LSSolverType		m_lsSolverType;
			QPSolverType		m_qpSolverType;
			StepSelType			m_stepSelType;
			int					m_lineSearchIters;
			Real				m_lineSearchFactor;
			Real				m_maxStepSize;

			Options()
			{
				m_stepSelType = StepSelType::SS_LineSearch;
				m_lsSolverType = LSSolverType::LS_EigenLDLT;
				m_qpSolverType = QPSolverType::QP_Newton;
				m_lineSearchIters = 10;
				m_lineSearchFactor = 0.5;
				m_maxStepSize = 1.0;
			}
		};

	protected:

		ILinearSolver* m_pLSolver;
		VectorXd m_dxPrev;

		// For matrix-based BFGS

		FastMatrixSd m_mBFGS;
		VectorXd m_vx_P;
		VectorXd m_vg_P;
		VectorXd m_vs;
		VectorXd m_vy;
		MatrixXd m_mBFGSMod_A;
		MatrixXd m_mBFGSMod_B;
		bool m_isBFGSInit;

		// For L-BFGS step

		typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> Matrix;
		typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> Vector;
		typedef Eigen::Map<Vector> MapVec;

		Matrix                    m_mS;      // History of the s vectors
		Matrix                    m_mY;      // History of the y vectors
		Vector                    m_vys;     // History of the s'y values
		Vector                    m_va;		 // History of the step lengths
		int						  m_K;		 // Last position
		int						  m_N;		 // Problem size
		int						  m_M;		 // Buffer size

	public:

		OptimSolver_SQP();
		virtual void Init();
		virtual ~OptimSolver_SQP();

		virtual Options& GetOptions() override { return *static_cast<Options*>(this->m_pOptions.get()); }

	protected:

		virtual bool ComputeStep(const VectorXd& vx, const VectorXd& vg, VectorXd& dx);
		virtual bool ComputeStep_NewtonScaled(const VectorXd& vx, const VectorXd& vg, VectorXd& dx);
		virtual bool ComputeStep_NewtonAligned(const VectorXd& vx, const VectorXd& vg, VectorXd& dx);
		virtual bool ComputeStep_NewtonDefault(const VectorXd& vx, const VectorXd& vg, VectorXd& dx);
		virtual bool ComputeStep_BFGSDirect(const VectorXd& vx, const VectorXd& vg, VectorXd& dx);
		virtual bool ComputeStep_BFGSInverse(const VectorXd& vx, const VectorXd& vg, VectorXd& dx);
		virtual bool ComputeStep_LBFGS(const VectorXd& vx, const VectorXd& vg, VectorXd& dx);

		virtual bool UpdateBFGS_Inverse(const VectorXd& vxNew, const VectorXd& vgNew);
		virtual bool UpdateBFGS_Direct(const VectorXd& vxNew, const VectorXd& vgNew);
		virtual bool InitialBFGS(const VectorXd& vxNew, const VectorXd& vgNew);

	};
}
