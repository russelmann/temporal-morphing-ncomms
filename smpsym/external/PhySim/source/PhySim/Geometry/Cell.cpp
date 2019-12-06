//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Geometry/Cell.h>

#include <PhySim/Geometry/Node.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	Cell::Cell(const vector<Node*>& m_vnodes) : Polytope()
	{
		this->m_vnodes = m_vnodes;
	}

	Cell::~Cell(void)
	{
		// Nothing to do here...
	}

}