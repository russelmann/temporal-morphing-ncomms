//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================


#include <PhySim/Geometry/Edge.h>
#include <PhySim/Geometry/Node.h>
#include <PhySim/Geometry/Face.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	Edge::Edge(const vector<Node*>& m_vnodes)
	{
		this->m_vnodes = m_vnodes;
	}

	Edge::Edge(Node* node1, Node* node2)
	{
		this->m_vnodes.push_back(node1);
		this->m_vnodes.push_back(node2);
		this->m_vedges.push_back(this);
	}

	Edge::~Edge(void) 
	{
		// Nothing to do here...
	}

	inline Real Edge::ComputeVolume(Space s) const
	{
		return (m_vnodes[1]->DoF()->GetPosition(s) - m_vnodes[0]->DoF()->GetPosition(s)).norm();
	}

}