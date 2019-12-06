//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Geometry/Cell_Hexa.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	Cell_Hexa::Cell_Hexa(const vector<Node*>& m_vnodes) : Cell(m_vnodes)
	{
		// Nothing to do here...
	}

	Cell_Hexa::~Cell_Hexa(void)
	{
		// Nothing to do here...
	}

	Real Cell_Hexa::ComputeVolume(Space s) const
	{
      assert(false);

		return -1;
	}

	void Cell_Hexa::ComputeShapeFunction(const VectorXd& vp, VectorXd& vN, Space s) const
	{
      assert(false);
	}

	void Cell_Hexa::ComputeShapeDerivative(const VectorXd& vp, MatrixXd& mB, Space s) const
	{
      assert(false);
	}

	void Cell_Hexa::GetQuadrature(int num, vector<VectorXd>& vQPos, vector<Real>& vQWei) const
	{
		if (num == 1)
		{
			vQPos.resize(1);
			vQWei.resize(1);

			vQWei[0] = 1;

			// Centroid in iso-parametric
			vQPos[0] = Vector3d::Zero(3);
		}
		else if (num == 8)
		{
         assert(false);
		}
		else
		{
         assert(false);
		}
	}

}
