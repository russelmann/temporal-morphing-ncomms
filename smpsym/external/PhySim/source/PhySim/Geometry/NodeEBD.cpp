//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Geometry/NodeEBD.h>

#include <PhySim/Geometry/Node.h>
#include <PhySim/Geometry/DoFSet.h>
#include <PhySim/Geometry/Polytope.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	NodeEBD::NodeEBD()
	{
		m_valid = false;

		this->m_pMaster = NULL;
		this->m_visopar.setZero();
	}

	NodeEBD::NodeEBD(const NodeEBD& toCopy)
	{
		this->m_pMaster = toCopy.m_pMaster;
		this->m_visopar = toCopy.m_visopar;
		this->m_valid = toCopy.m_valid;
	}

	NodeEBD::NodeEBD(Polytope* pMaster, const VectorXd& visopar)
	{
		this->m_pMaster = pMaster;
		this->m_visopar = visopar;
		this->m_valid = this->m_pMaster->IsValidParametric(m_visopar);
	}


	NodeEBD::~NodeEBD(void)
	{
		// Nothing to do here...
	}

	Vector3d NodeEBD::GetPosition(Space s) const
	{
		return this->m_pMaster->InterpolatePosition(this->m_visopar, s);
	}

}