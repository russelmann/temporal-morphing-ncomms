//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Geometry/Cell_Poly.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	Cell_Poly::Cell_Poly(const vector<Node*>& m_vnodes) : Cell(m_vnodes)
	{
		// Nothing to do here...
	}

	Cell_Poly::~Cell_Poly(void)
	{
		// Nothing to do here...
	}

	Real Cell_Poly::ComputeVolume(Space s) const
	{
      assert(false);

		return -1;
	}

	void Cell_Poly::ComputeShapeFunction(const VectorXd& vp, VectorXd& vN, Space s) const
	{
      assert(false);
	}

	void Cell_Poly::ComputeShapeDerivative(const VectorXd& vp, MatrixXd& mB, Space s) const
	{
      assert(false);
	}

}
