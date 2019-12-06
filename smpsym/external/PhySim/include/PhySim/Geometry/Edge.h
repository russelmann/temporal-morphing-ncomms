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

	class Edge : public Polytope
	{

	public:
		Edge(const vector<Node*>& vnodes);
		Edge(Node* node1, Node* node2);

		~Edge(void);

		inline Node* GetOrigin(void) const { return m_vnodes[0]; }
		inline void SetOrigin(Node* node) { m_vnodes[0] = node; }

		inline Node* GetHead(void) const { return m_vnodes[1]; }
		inline void SetHead(Node* node) { m_vnodes[1] = node; }

		inline virtual int Order() { return 1; }

		virtual Real ComputeVolume(Space s = Space::DEF) const;

	};
}
