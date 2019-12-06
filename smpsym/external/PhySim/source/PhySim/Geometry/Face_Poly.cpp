//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Geometry/Face_Poly.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	Face_Poly::Face_Poly(const vector<Node*>& m_vnodes) : Face(m_vnodes)
	{
		// Nothing to do here...
	}

	Face_Poly::~Face_Poly(void)
	{
		// Nothing to do here...
	}

	Real Face_Poly::ComputeVolume(Space s) const
	{
      assert(false);

		return -1;
	}

	void Face_Poly::ComputeShapeFunction(const VectorXd& vp, VectorXd& vN, Space s) const
	{
		assert(false);
	}

	void Face_Poly::ComputeShapeDerivative(const VectorXd& vp, MatrixXd& mB, Space s) const
	{
		assert(false);
	}

	void Face_Poly::GetQuadrature(int num, vector<VectorXd>& vQPos, vector<Real>& vQWei) const
	{
		if (num == 1)
		{
			vQPos.resize(1);
			vQWei.resize(1);

			vQWei[0] = 1;

			// Centroid in iso-parametric
			vQPos[0] = Vector3d::Zero();
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
