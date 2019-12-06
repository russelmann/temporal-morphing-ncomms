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

#include <PhySim/Geometry/Cell.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class Cell_Poly : public Cell
	{

	public:
		Cell_Poly(const vector<Node*>& vnodes);

		virtual ~Cell_Poly(void);

		virtual Real ComputeVolume(Space s = Space::DEF) const;
		virtual void ComputeShapeFunction(const VectorXd& vp, VectorXd& vN, Space s = Space::MAT) const;
		virtual void ComputeShapeDerivative(const VectorXd& vp, MatrixXd& vN, Space s = Space::MAT) const;

	};
}