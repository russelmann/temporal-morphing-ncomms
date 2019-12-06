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

	bool BCondition::StepLoading()
	{
		if (m_curStep == this->m_setup.m_maxStep)
			return false; // Already at max step

		if (this->m_pModel->GetGradient().norm() < this->m_setup.m_maxError)
		{
			// Next BC counter
			this->m_curStep++;

			// Is this a critical step? Should we apply BC?
			return m_curStep % this->m_setup.m_incStep == 0;
		}
		else
		{
			logSimu("\n[TRACE] Not ready to load...");

			return false;
		}
	}

	bool BCondition::ResetLoading()
	{
		if (this->m_curStep == 0)
			return false;
			
		this->m_curStep = 0;
		return true;
	}

	bool BCondition::FullLoading()
	{
		if (this->m_curStep == this->m_setup.m_maxStep)
			return false;

		this->m_curStep = this->m_setup.m_maxStep;
		return true;
	}

}

