//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#pragma once

#include <PhySim/CommonIncludes.h>
#include <PhySim/PhySimInterface.h>

#include <PhySim/Geometry/Polytope.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class Cell : public Polytope
	{

	public:
		Cell(const vector<Node*>& vnodes);

		virtual ~Cell(void);

		inline virtual int Order() { return 3; }

	};
}
