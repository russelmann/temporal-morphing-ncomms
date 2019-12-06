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

#include <PhySim/Solvers/OptimProblem.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class OptimProblem_Subproblem : public OptimProblem
	{
	public:
		IOptimProblem*	m_pProblem;
		iVector			m_varSel;
		MatrixSd		m_selMat;

		OptimProblem_Subproblem()
		{
			// TODO: Remove
		}

		OptimProblem_Subproblem(IOptimProblem* pProblem, const iVector& varSel)
		{
			assert(pProblem != NULL);
			assert(varSel.size() != 0);

			this->m_pProblem = pProblem;
			this->m_varSel = varSel;

			int numVars = (int)varSel.size();

			VectorTd vT;
			vT.resize(numVars);
			for (int j = 0; j < numVars; ++j)
				vT[j] = Triplet<Real>(j, varSel[j], 1.0);

			m_selMat = MatrixSd(numVars, m_pProblem->GetNumVariables());
			m_selMat.setFromTriplets(vT.begin(), vT.end());
			m_selMat.makeCompressed();
		}

		virtual int GetNumVariables() const { return (int) this->m_varSel.size(); }

		virtual void GetVariables(VectorXd& vx) const
		{
			VectorXd vxFull;
			m_pProblem->GetVariables(vxFull);
			vx = this->m_selMat*vxFull;
		}

		virtual void SetVariables(const VectorXd& vx)
		{
			VectorXd v0Full;
			m_pProblem->GetVariables(v0Full);
			VectorXd dx = vx - this->m_selMat*v0Full;
			VectorXd dxFull = m_selMat.transpose()*dx;
			m_pProblem->SetVariables(v0Full + dxFull);
		}

		virtual bool GetEnergy(Real& e)
		{
			bool result = this->m_pProblem->GetEnergy(e);

			return result;
		}

		virtual bool GetGradient(VectorXd& vg)
		{
			VectorXd vgFull;
			bool result = m_pProblem->GetGradient(vgFull);
			double normFull = vgFull.norm();
			vg = this->m_selMat*vgFull;
			double norm = vg.norm();
			return result;
		}

		virtual bool GetHessian(MatrixSd& mH)
		{
			MatrixSd mHFull;
			bool result = m_pProblem->GetHessian(mHFull);
			mH = m_selMat*mHFull*m_selMat.transpose();
			return result;
		}

		virtual bool IsFullyConstrained() override
		{
			return m_pProblem->IsFullyConstrained();
		}

		virtual bool StartSolveCallback() override
		{
			return m_pProblem->StartSolveCallback();
		};

		virtual bool StartStepCallback() override
		{
			return m_pProblem->StartStepCallback();
		};

		virtual bool EndSolveCallback() override
		{
			return m_pProblem->EndSolveCallback();
		};

		virtual bool EndStepCallback() override
		{
			return m_pProblem->EndStepCallback();
		};

		virtual bool PreComputeStepCallback() override
		{
			return m_pProblem->PreComputeStepCallback();
		};

		virtual bool PosComputeStepCallback() override
		{
			return m_pProblem->PosComputeStepCallback();
		};

		virtual bool PrePerformStepCallback() override
		{
			return m_pProblem->PrePerformStepCallback();
		};

		virtual bool PosPerformStepCallback() override
		{
			return m_pProblem->PosPerformStepCallback();
		};

	};
}