//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================


#include <PhySim/Geometry/Polytope.h>
#include <PhySim/Geometry/Node.h>

#include <PhySim/Geometry/NodeEBD.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	Real Polytope::ComputeVolume(Space s) const
	{
		assert(false);

		return -1;
	}

	void Polytope::GetNodeMatrix(MatrixXd& mN, Space s) const
	{
		int N = (int) this->m_vnodes.size();

		mN.resize(3, N);
		for (int i = 0; i < N; ++i)
			mN.col(i) = m_vnodes[i]->DoF()->GetPosition(s);
	}

	void Polytope::SetNodeMatrix(const MatrixXd& mN, Space s)
	{
		int N = (int) this->m_vnodes.size();

		assert(mN.rows() == 3);
		assert(mN.cols() == N);

		for (int i = 0; i < N; ++i)
			m_vnodes[i]->DoF()->SetPosition(mN.col(i), s);
	}

	Vector3d Polytope::ComputeCentroid(Space s) const
	{
		Vector3d vsum = Vector3d::Zero();

		size_t numNode = this->m_vnodes.size();
		for (size_t i = 0; i < numNode; ++i)
			vsum += this->m_vnodes[i]->DoF()->GetPosition(s);

		return (1.0 / numNode)*vsum;
	}

	Matrix3d Polytope::ComputeRotation(Space f, Space t) const
	{
		assert(false);

		return Matrix3d();
	}

	Vector3d Polytope::InterpolatePosition(const VectorXd& vp, Space s) const
	{
		VectorXd vN;		
		Vector3d vx;
		this->ComputeShapeFunction(vp, vN);

		vx.setZero();

		size_t numNode = this->m_vnodes.size();
		for (size_t i = 0; i < numNode; ++i)
			vx += vN(i)*this->m_vnodes[i]->DoF()->GetPosition(s);

		return vx;
	}

	void Polytope::GetQuadrature(int num, vector<VectorXd>& vQPos, vector<Real>& vQWei) const
	{
		assert(false);
	}

	void Polytope::ComputeShapeFunction(const VectorXd& vp, VectorXd& vN, Space s) const
	{
		assert(false);
	}

	void Polytope::ComputeShapeDerivative(const VectorXd& vp, MatrixXd& mB, Space s) const
	{
		assert(false);
	}

	void Polytope::ComputeNat2IsoTransform(VectorXd& vb, MatrixXd& mA, Space s) const
	{
		assert(false);
	}

	void Polytope::ComputeIso2NatTransform(VectorXd& vb, MatrixXd& mA, Space s) const
	{
		assert(false);
	}

	NodeEBD Polytope::ComputeEmbedding(const Vector3d& vp, Space s)
	{
		VectorXd vb;
		MatrixXd mA;
		ComputeNat2IsoTransform(vb, mA, s);
		VectorXd viso = mA*vp + vb;
		return NodeEBD(this, viso);
	}

	bool Polytope::IsValidParametric(const VectorXd& vp) const
	{
		assert(false);

		return false;
	}

	void Polytope::SetSubelementPositions(Space s)
	{
		assert(false);
	}

}
