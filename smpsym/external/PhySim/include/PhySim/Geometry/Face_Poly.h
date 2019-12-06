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

	class Face_Poly : public Face
	{

	public:
		Face_Poly(const vector<Node*>& vnodes);

		virtual ~Face_Poly(void);

		virtual Real ComputeVolume(Space s = Space::DEF) const;
		virtual void ComputeShapeFunction(const VectorXd& vp, VectorXd& vN, Space s = Space::MAT) const;
		virtual void ComputeShapeDerivative(const VectorXd& vp, MatrixXd& vN, Space s = Space::MAT) const;

		virtual void GetQuadrature(int num, vector<VectorXd>& vQPos, vector<Real>& vQWei) const;

	};
}