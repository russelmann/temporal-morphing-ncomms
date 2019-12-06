//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Solvers/OptimProblem_BasicDynamic.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	OptimProblem_BasicDynamic::OptimProblem_BasicDynamic(IModel* pM, Real dt)
	{
		this->Init(pM, dt);
	}

	OptimProblem_BasicDynamic::~OptimProblem_BasicDynamic()
	{
		// Nothing to do here..
	}

	void OptimProblem_BasicDynamic::Init(IModel* pM, Real dt)
	{
		this->m_timeStep = dt;
		this->m_pModel = pM;

		this->m_name = "[BasicSynamic]";
	}

	void OptimProblem_BasicDynamic::GetVariables(VectorXd& vx) const
	{
		m_pModel->GetFreeDOFPosition(vx);
	}

	void OptimProblem_BasicDynamic::SetVariables(const VectorXd& vx)
	{
		m_pModel->SetFreeDOFPosition(vx);
	}

	bool OptimProblem_BasicDynamic::GetEnergy(Real& e)
	{
		// TODO

		return false;
	}

	bool OptimProblem_BasicDynamic::GetGradient(VectorXd& vg)
	{
		// TODO

		return false;
	}

	bool OptimProblem_BasicDynamic::GetHessian(MatrixSd& mH)
	{
		// TODO

		return false;
	}

}
