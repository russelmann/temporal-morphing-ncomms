//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Geometry/Node.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	Node::Node(const Vector3d& pos)
	{
		Vector3d zero;
		zero.setZero();
		this->m_dofs = new DoFSet(3);
		this->m_dofs->SetPositions(pos);
		this->m_dofs->SetVelocities(zero);
		this->m_vnodes.push_back(this);
	}

	Node::~Node(void)
	{
		delete this->m_dofs;
	}

}