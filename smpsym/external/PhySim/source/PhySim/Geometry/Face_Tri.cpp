//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================


#include <PhySim/Geometry/Edge.h>
#include <PhySim/Geometry/Node.h>
#include <PhySim/Geometry/Face_Tri.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	Face_Tri::Face_Tri(const vector<Node*>& m_vnodes) : Face(m_vnodes)
	{
		m_vfaces.push_back(this);
	}

	Face_Tri::~Face_Tri(void) { }

	Real Face_Tri::ComputeVolume(Space s) const
	{
		Vector3d e0 = (m_vnodes[1]->DoF()->GetPosition(s) - m_vnodes[0]->DoF()->GetPosition(s));
		Vector3d e1 = (m_vnodes[2]->DoF()->GetPosition(s) - m_vnodes[0]->DoF()->GetPosition(s));
		return 0.5*e0.cross(e1).norm();
	}

	void Face_Tri::ComputeShapeFunction(const VectorXd& vp, VectorXd& vN, Space s) const
	{
		// Basic implementation: linear interpolation

		vN.resize(3);
		const Vector3d& va = this->GetNode0()->DoF()->GetPosition(s);
		const Vector3d& vb = this->GetNode1()->DoF()->GetPosition(s);
		const Vector3d& vc = this->GetNode2()->DoF()->GetPosition(s);
		Vector3d v0 = vb - va, v1 = vc - va, v2 = vp - va;
		Real d00 = v0.dot(v0);
		Real d01 = v0.dot(v1);
		Real d11 = v1.dot(v1);
		Real d20 = v2.dot(v0);
		Real d21 = v2.dot(v1);
		Real D = d00 * d11 - d01 * d01;
		vN(1) = (d11 * d20 - d01 * d21) / D;
		vN(2) = (d00 * d21 - d01 * d20) / D;
		vN(0) = 1.0f - vN(1) - vN(2);
	}

	void Face_Tri::ComputeShapeDerivative(const VectorXd& vp, MatrixXd& mB, Space s) const
	{
		mB.resize(3, 2);
		mB.row(0) = Vector2d(-1, -1);
		mB.row(1) = Vector2d(1, 0);
		mB.row(2) = Vector2d(0, 1);
	}

	void Face_Tri::ComputeNat2IsoTransform(VectorXd& b, MatrixXd& A, Space s) const
	{
		MatrixXd mNB;
		mNB.resize(3, 3);
		Vector3d e0 = this->m_vnodes[1]->DoF()->GetPosition(s) - this->m_vnodes[0]->DoF()->GetPosition(s);
		Vector3d e1 = this->m_vnodes[2]->DoF()->GetPosition(s) - this->m_vnodes[0]->DoF()->GetPosition(s);
		mNB.col(0) = e0;
		mNB.col(1) = e1;
		mNB.col(2) = e0.cross(e1);
		A = mNB.inverse().block(0, 0, 2, 3);
		b = -A*this->m_vnodes[0]->DoF()->GetPosition(s);
	}

	void Face_Tri::ComputeIso2NatTransform(VectorXd& b, MatrixXd& A, Space s) const
	{
		A.resize(3, 2);
		A.col(0) = this->m_vnodes[1]->DoF()->GetPosition(s) - this->m_vnodes[0]->DoF()->GetPosition(s);
		A.col(1) = this->m_vnodes[2]->DoF()->GetPosition(s) - this->m_vnodes[0]->DoF()->GetPosition(s);
		b = this->m_vnodes[0]->DoF()->GetPosition(s);
	}

	void Face_Tri::GetQuadrature(int num, vector<VectorXd>& vQPos, vector<Real>& vQWei) const
	{
		if (num == 1)
		{
			vQPos.resize(1);
			vQWei.resize(1);

			vQWei[0] = 1.0/2.0;

			// Centroid in iso-parametric coordinates
			vQPos[0] = (1.0/3.0)*Vector3d::Ones();
		}
		else
		{
			assert(false);
		}
	}

}
