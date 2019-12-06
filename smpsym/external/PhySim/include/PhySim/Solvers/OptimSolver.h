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

	class OptimSolver : public IOptimSolver
	{
	public:
		struct Options
		{
			IOptimProblem*		m_pProblem;

			Real				m_maxError;
			int					m_maxIters;

			bool profileTime;

			Options()
			{
				m_pProblem = NULL;

				m_maxError = 1e-9;
				m_maxIters = 1000;

				profileTime = true;
			}
		};

	protected:
		shared_ptr<Options> m_pOptions;

		bool m_isInit;

		CustomTimer m_timerSolve;

		void(*m_stepCallback)(IOptimSolver* pSender, void* pReceiver);
		void(*m_fullCallback)(IOptimSolver* pSender, void* pReceiver);
		void* m_stepCallbackReceiver;
		void* m_fullCallbackReceiver;

		bool m_forceSolveStop;

	public:

		OptimSolver();
		virtual void Init();
		virtual ~OptimSolver();

		virtual IOptimProblem& Problem() override { return *this->m_pOptions->m_pProblem; }

		virtual Options& GetOptions();

		virtual SolveResult SolveFull();
		virtual SolveResult SolveStep();

		virtual Real ComputeOptimality();
		virtual Real ComputeFeasibility();
		virtual bool IsOptimal();
		virtual bool IsFeasible();

		virtual void RegisterCallback_Step(void(*amazingCallback)(IOptimSolver*, void*), void* pReceiver) override
		{
			this->m_stepCallback = amazingCallback;
			this->m_stepCallbackReceiver = pReceiver;
		}

		virtual void RegisterCallback_Full(void(*amazingCallback)(IOptimSolver*, void*), void* pReceiver) override
		{
			this->m_fullCallback = amazingCallback;
			this->m_fullCallbackReceiver = pReceiver;
		}

		virtual void OnEvent_AfterStep()
		{
			if (m_stepCallback != NULL)
				m_stepCallback(this, m_stepCallbackReceiver);
		}

		virtual void OnEvent_AfterFull()
		{
			if (m_fullCallback != NULL)
				m_fullCallback(this, m_stepCallbackReceiver);
		}

		virtual bool& ForceSolveStop() { return this->m_forceSolveStop; }

	};
}
