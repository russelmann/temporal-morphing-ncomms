//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Solvers/OptimProblem_BasicStatic.h>

#include <PhySim/Models/Model.h>

#include <PhySim/Models/Model_Particles.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	OptimProblem_BasicStatic::OptimProblem_BasicStatic(IModel* pM)
	{
		this->Init(pM);
	}

	OptimProblem_BasicStatic::~OptimProblem_BasicStatic()
	{
		this->m_pS.reset();
	}

	void OptimProblem_BasicStatic::Init(IModel* pM)
	{
		this->m_pM = pM;
		this->m_pS = m_pM->CreateState();
		this->m_N = m_pM->GetNumFreeDOF();
		this->m_M = 0;

		this->m_name = "[BasicStatic]";
	}

	void OptimProblem_BasicStatic::GetVariables(VectorXd& vx) const
	{
		m_pM->GetFreeDOFPosition(vx);
	}

	void OptimProblem_BasicStatic::SetVariables(const VectorXd& vx)
	{
		m_pM->SetFreeDOFPosition(vx);
	}

	bool OptimProblem_BasicStatic::GetEnergy(Real& e)
	{
		e = m_pM->GetEnergy();

		return true;
	}

	bool OptimProblem_BasicStatic::GetGradient(VectorXd& vg)
	{
		vg = m_pM->GetGradient();

		return true;
	}

	bool OptimProblem_BasicStatic::GetHessian(MatrixSd& mH)
	{
		mH = m_pM->GetHessian();

		return true;
	}

	bool OptimProblem_BasicStatic::StartStepCallback()
	{
		m_pM->StepBoundaryConditions();

		m_pM->GetState(m_pS);

		return true;
	}

	bool OptimProblem_BasicStatic::PrePerformStepCallback()
	{
		if (!m_pM->HasState(this->m_pS))
			m_pM->SetState(this->m_pS);

		return true;
	}

}
