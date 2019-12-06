//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Models/BC_Fixed.h>

#include <PhySim/Models/Model.h>
#include <PhySim/Models/BCondition.h>
#include <PhySim/Geometry/DoFSet.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	BC_Fixed::BC_Fixed(Model* pModel, const BCSetup& bcSetup) : BCondition(pModel, bcSetup)
	{
		this->m_pEle = NULL;
	}

	BC_Fixed::~BC_Fixed(void)
	{
		// Nothing to do here...
	}

	void BC_Fixed::Init()
	{
		for (size_t i = 0; i < this->m_setup.m_vDoF.size(); ++i)
		{
			Vector3d pos = this->m_setup.m_vDoF[i]->GetPosition_x();
			this->m_setup.m_vDoF[i]->Fix();
			//this->m_setup.m_vini[i] = pos;
		}

		this->Apply();
	}

	void BC_Fixed::Apply()
	{
		vector<VectorXd> vval;
		this->GetCurValues(vval);

		for (size_t i = 0; i < m_setup.m_vDoF.size(); ++i)
			this->m_setup.m_vDoF[i]->SetPosition_x(vval[i]);
			
		this->m_pModel->DirtyDeformed();
	}

}

