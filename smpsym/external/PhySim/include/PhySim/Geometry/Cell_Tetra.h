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

	class Cell_Tetra : public Cell
	{

	public:
		Cell_Tetra(const vector<Node*>& vnodes);

		virtual ~Cell_Tetra(void);

		virtual Real ComputeVolume(Space s = Space::DEF) const;
		virtual void ComputeShapeFunction(const VectorXd& vs, VectorXd& vN, Space s = Space::MAT) const;
		virtual void ComputeShapeDerivative(const VectorXd& vs, MatrixXd& vN, Space s = Space::MAT) const;

		virtual void ComputeNat2IsoTransform(VectorXd& vb, MatrixXd& mA, Space s = Space::MAT) const override;
		virtual void ComputeIso2NatTransform(VectorXd& vb, MatrixXd& mA, Space s = Space::MAT) const override;

		virtual void GetQuadrature(int num, vector<VectorXd>& vQPos, vector<Real>& vQWei) const;

		virtual bool IsValidParametric(const VectorXd& vp) const override;

		virtual void SetSubelementPositions(Space s) override;

	};
}