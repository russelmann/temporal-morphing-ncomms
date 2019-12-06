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
#include <PhySim/Geometry/DoFSet.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class Node : public Polytope
	{

	public:
		Node(const Vector3d& pos = Vector3d(0.0, 0.0, 0.0));
		virtual ~Node(void);

		inline int Order() { return 1; }

		virtual Real ComputeVolume(Space s = Space::DEF) const { return 0; }

	};
}
