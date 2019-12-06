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

#include <PhySim/Geometry/Face.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class Face_Tri : public Face
	{

	public:
		Face_Tri(const vector<Node*>& vnodes);

		~Face_Tri(void);

		inline Node* GetNode0(void) const { return m_vnodes[0]; }
		inline void SetNode0(Node* node) { m_vnodes[0] = node; }

		inline Node* GetNode1(void) const { return m_vnodes[1]; }
		inline void SetNode1(Node* node) { m_vnodes[1] = node; }

		inline Node* GetNode2(void) const { return m_vnodes[2]; }
		inline void SetNode2(Node* node) { m_vnodes[2] = node; }

		virtual Real ComputeVolume(Space s = Space::DEF) const;
		virtual void ComputeShapeFunction(const VectorXd& vp, VectorXd& vN, Space s = Space::MAT) const;
		virtual void ComputeShapeDerivative(const VectorXd& vp, MatrixXd& vN, Space s = Space::MAT) const;

		virtual void ComputeNat2IsoTransform(VectorXd& b, MatrixXd& A, Space s = Space::MAT) const;
		virtual void ComputeIso2NatTransform(VectorXd& b, MatrixXd& A, Space s = Space::MAT) const;

		virtual void GetQuadrature(int num, vector<VectorXd>& vQPos, vector<Real>& vQWei) const;

	};
}
