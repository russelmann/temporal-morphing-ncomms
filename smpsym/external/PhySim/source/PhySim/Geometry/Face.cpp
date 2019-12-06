//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Geometry/Face.h>

#include <PhySim/Geometry/Node.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	Face::Face(const vector<Node*>& m_vnodes) : Polytope()
	{
		this->m_vnodes = m_vnodes;
	}

	Face::~Face(void)
	{
		// Nothing to do here...
	}

}