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

	class OptimProblem : public IOptimProblem
	{
	protected:
		int m_N;
		int m_M;
		string m_name;

	public:

		virtual const string& GetName() const { return this->m_name; }

		virtual int GetNumVariables() const { return this->m_N; }
		virtual int GetNumConstraints() const { return this->m_M; }

		virtual void GetVariables(VectorXd& vx) const = 0;
		virtual void SetVariables(const VectorXd& vx) = 0;
		
		virtual void GetUpperBound(VectorXd& vub) const;
		virtual void GetLowerBound(VectorXd& vlb) const;

		virtual bool GetEnergy(Real& e);
		virtual bool GetGradient(VectorXd& vg);
		virtual bool GetHessian(MatrixSd& mH);

		virtual bool GetConstraints(VectorXd& vc);
		virtual bool GetJacobian(MatrixSd& mJ);

		virtual bool IsFullyConstrained() override
		{
			return true;
		}

		virtual bool StartSolveCallback() override
		{
			// Default implementation, do nothing
			return true;
		};

		virtual bool StartStepCallback() override
		{
			// Default implementation, do nothing
			return true;
		};

		virtual bool EndSolveCallback() override
		{
			// Default implementation, do nothing
			return true;
		};

		virtual bool EndStepCallback() override
		{
			// Default implementation, do nothing
			return true;
		};

		virtual bool PreComputeStepCallback() override
		{
			// Default implementation, do nothing
			return true;
		};

		virtual bool PosComputeStepCallback() override
		{
			// Default implementation, do nothing
			return true;
		};

		virtual bool PrePerformStepCallback() override
		{
			// Default implementation, do nothing
			return true;
		};

		virtual bool PosPerformStepCallback() override
		{
			// Default implementation, do nothing
			return true;
		};

	};
}