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

	class Cell_Hexa : public Cell
	{

	public:
		Cell_Hexa(const vector<Node*>& mvnodes);

		virtual ~Cell_Hexa(void);

		virtual Real ComputeVolume(Space s = Space::DEF) const;
		virtual void ComputeShapeFunction(const VectorXd& vp, VectorXd& vN, Space s = Space::MAT) const;
		virtual void ComputeShapeDerivative(const VectorXd& vp, MatrixXd& vN, Space s = Space::MAT) const;

		virtual void GetQuadrature(int num, vector<VectorXd>& vQPos, vector<Real>& vQWei) const;

	};
}