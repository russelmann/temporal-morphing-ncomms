//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Models/BC_Gravity.h>

#include <PhySim/Models/Model.h>
#include <PhySim/Models/BCondition.h>

#include <PhySim/Energies/EnergyElement_Gravity.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	BC_Gravity::BC_Gravity(Model* pModel, const BCSetup& bcSetup) : BCondition(pModel, bcSetup)
	{
		this->m_pEle = this->m_pEleGravity = new EnergyElement_Gravity((Model_Particles*)pModel);
	}

	BC_Gravity::~BC_Gravity(void)
	{
		// Nothing to do here...
	}

	void BC_Gravity::Init()
	{
		this->Apply();
	}

	void BC_Gravity::Apply()
	{
		vector<VectorXd> vval;
		this->GetCurValues(vval);

		this->m_pEleGravity->SetGravityAcceleration(vval[0]);

		this->m_pModel->DirtyDeformed();
	}

}

