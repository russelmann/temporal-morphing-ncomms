//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Models/BC_Force.h>

#include <PhySim/Models/Model.h>
#include <PhySim/Models/BCondition.h>

#include <PhySim/Energies/EnergyElement_Force.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	BC_Force::BC_Force(Model* pModel, const BCSetup& bcSetup) : BCondition(pModel, bcSetup)
	{
		this->m_pEle = this->m_pEleForce = new EnergyElement_Force((Model_Particles*)pModel, bcSetup.m_vDoF);
	}

	BC_Force::~BC_Force(void)
	{
		// Nothing to do here...
	}

	void BC_Force::Init()
	{
		this->Apply();
	}

	void BC_Force::Apply()
	{
		vector<VectorXd> vval;
		this->GetCurValues(vval);

		m_pEleForce->SetForces(vval);
		this->m_pModel->DirtyDeformed();
	}

}

