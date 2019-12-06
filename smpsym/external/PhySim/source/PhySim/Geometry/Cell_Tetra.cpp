//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Geometry/Cell_Tetra.h>

#include <PhySim/Geometry/Node.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	Cell_Tetra::Cell_Tetra(const vector<Node*>& m_vnodes) : Cell(m_vnodes)
	{
		// Nothing to do here...
	}

	Cell_Tetra::~Cell_Tetra(void)
	{
		// Nothing to do here...
	}

	Real Cell_Tetra::ComputeVolume(Space s) const
	{
		Matrix3d mV;
		mV.col(0) = this->m_vnodes[1]->DoF()->GetPosition(s) - this->m_vnodes[0]->DoF()->GetPosition(s);
		mV.col(1) = this->m_vnodes[2]->DoF()->GetPosition(s) - this->m_vnodes[0]->DoF()->GetPosition(s);
		mV.col(2) = this->m_vnodes[3]->DoF()->GetPosition(s) - this->m_vnodes[0]->DoF()->GetPosition(s);
		return (1.0/6.0)*mV.determinant();
	}

	void Cell_Tetra::ComputeShapeFunction(const VectorXd& vp, VectorXd& vN, Space s) const
	{
		// 4-Node Linear tetrahedron

		if (this->m_vnodes.size() == 4)
		{
			vN.resize(4);
			vN(0) = 1.0 - vp(0) - vp(1) - vp(2);
			vN(1) = vp(0);
			vN(2) = vp(1);
			vN(3) = vp(2);
			return;
		}

		// 10-Node Quadratic tetrahedron

		if (this->m_vnodes.size() == 10)
		{
			vN.resize(10);
			Real t0 = 1.0 - vp(0) - vp(1) - vp(2);
			Real t1 = vp(0);
			Real t2 = vp(1);
			Real t3 = vp(2);
			vN(0) = t0*(2*t0 - 1);
			vN(1) = t1*(2*t1 - 1);
			vN(2) = t2*(2*t2 - 1);
			vN(3) = t3*(2*t3 - 1);
			vN(4) = 4*t0*t1;
			vN(5) = 4*t1*t2;
			vN(6) = 4*t2*t0;
			vN(7) = 4*t0*t3;
			vN(8) = 4*t1*t3;
			vN(9) = 4*t2*t3;
			return;
		}

	}

	void Cell_Tetra::ComputeShapeDerivative(const VectorXd& vp, MatrixXd& mB, Space s) const
	{
		// 4-Node Linear tetrahedron

		if (this->m_vnodes.size() == 4)
		{
			mB.resize(4, 3);
			mB.row(0) = Vector3d(-1, -1, -1);
			mB.row(1) = Vector3d(1, 0, 0);
			mB.row(2) = Vector3d(0, 1, 0);
			mB.row(3) = Vector3d(0, 0, 1);

			return;
		}

		// 10-Node Quadratic tetrahedron

		if (this->m_vnodes.size() == 10)
		{
			Real t0 = 1.0 - vp(0) - vp(1) - vp(2);
			Real t1 = vp(0);
			Real t2 = vp(1);
			Real t3 = vp(2);

			MatrixXd DNDt(10, 4);
			DNDt.row(0) = Vector4d(4*t0 - 1, 0, 0, 0);
			DNDt.row(1) = Vector4d(0, 4*t1 - 1, 0, 0);
			DNDt.row(2) = Vector4d(0, 0, 4*t2 - 1, 0);
			DNDt.row(3) = Vector4d(0, 0, 0, 4*t3 - 1);
			DNDt.row(4) = Vector4d(4*t1, 4*t0, 0, 0);
			DNDt.row(5) = Vector4d(0, 4*t2, 4*t1, 0);
			DNDt.row(6) = Vector4d(4*t2, 0, 4*t0, 0);
			DNDt.row(7) = Vector4d(4*t3, 0, 0, 4*t0);
			DNDt.row(8) = Vector4d(0, 4*t3, 0, 4*t1);
			DNDt.row(9) = Vector4d(0, 0, 4*t3, 4*t2);

			MatrixXd DtDp(4, 3);
			DtDp.row(0) = Vector3d(-1, -1, -1);
			DtDp.row(1) = Vector3d(1, 0, 0);
			DtDp.row(2) = Vector3d(0, 1, 0);
			DtDp.row(3) = Vector3d(0, 0, 1);

			mB = DNDt*DtDp;

			return;
		}
	}

	void Cell_Tetra::ComputeNat2IsoTransform(VectorXd& b, MatrixXd& A, Space s) const
	{
		MatrixXd mNB;
		mNB.resize(3, 3);
		mNB.col(0) = this->m_vnodes[1]->DoF()->GetPosition(s) - this->m_vnodes[0]->DoF()->GetPosition(s);
		mNB.col(1) = this->m_vnodes[2]->DoF()->GetPosition(s) - this->m_vnodes[0]->DoF()->GetPosition(s);
		mNB.col(2) = this->m_vnodes[3]->DoF()->GetPosition(s) - this->m_vnodes[0]->DoF()->GetPosition(s);
		A = mNB.inverse();
		b = -A*this->m_vnodes[0]->DoF()->GetPosition(s);
	}

	void Cell_Tetra::ComputeIso2NatTransform(VectorXd& b, MatrixXd& A, Space s) const
	{
		A.resize(3, 3);
		A.col(0) = this->m_vnodes[1]->DoF()->GetPosition(s) - this->m_vnodes[0]->DoF()->GetPosition(s);
		A.col(1) = this->m_vnodes[2]->DoF()->GetPosition(s) - this->m_vnodes[0]->DoF()->GetPosition(s);
		A.col(2) = this->m_vnodes[3]->DoF()->GetPosition(s) - this->m_vnodes[0]->DoF()->GetPosition(s);
		b = this->m_vnodes[0]->DoF()->GetPosition(s);
	}

	void Cell_Tetra::GetQuadrature(int num, vector<VectorXd>& vQPos, vector<Real>& vQWei) const
	{
		// 4-Node Linear tetrahedron

		if (this->m_vnodes.size() == 4)
		{
			vQPos.resize(1);
			vQWei.resize(1);

			vQWei[0] = 1.0/6.0;

			// Iso-parametric centroid, for instance
			vQPos[0] = (1.0 / 4.0)*Vector3d::Ones(3);

			return;
		}

		// 10-Node Linear tetrahedron

		if (this->m_vnodes.size() == 10)
		{
			if (num == 1)
			{
				vQPos.resize(1);
				vQWei.resize(1);

				vQWei[0] = 1.0 / 6.0;

				// Iso-parametric centroid, for instance
				vQPos[0] = (1.0 / 4.0)*Vector3d::Ones(3);

				return;
			}

			// Refs 4,5,11
			// P Keast, Moderate degree tetrahedral quadrature formulas, CMAME 55: 
			// 339 - 348 (1986) O.C.Zienkiewicz, The Finite Element Method, Sixth Edition,
			
			// Refs 8
			// Gauss Legendre Quadrature Formulae for Tetrahedra
			// H. T. Rathod , B. Venkatesudu & K. V. Nagaraja

			if (num == 4)
			{
				vQPos.resize(4);
				vQWei.resize(4);

				vQPos[0] = Vector3d(0.5854101966249685, 0.1381966011250105, 0.1381966011250105);
				vQPos[1] = Vector3d(0.1381966011250105, 0.1381966011250105, 0.1381966011250105);
				vQPos[2] = Vector3d(0.1381966011250105, 0.1381966011250105, 0.5854101966249685);
				vQPos[3] = Vector3d(0.1381966011250105, 0.5854101966249685, 0.1381966011250105);
				vQWei[0] = 0.25 / 6.0;
				vQWei[1] = 0.25 / 6.0;
				vQWei[2] = 0.25 / 6.0;
				vQWei[3] = 0.25 / 6.0;

				return;
			}

			if (num == 5)
			{
				vQPos.resize(5);
				vQWei.resize(5);
				vQPos[0] = Vector3d(0.250000000000, 0.250000000000, 0.250000000000);
				vQPos[1] = Vector3d(0.500000000000, 0.166666666667, 0.166666666667);
				vQPos[2] = Vector3d(0.166666666667, 0.166666666667, 0.166666666667);
				vQPos[3] = Vector3d(0.166666666667, 0.166666666667, 0.500000000000);
				vQPos[4] = Vector3d(0.166666666667, 0.500000000000, 0.166666666667);
				vQWei[0] = -0.8 / 6.0;
				vQWei[1] = 0.45 / 6.0;
				vQWei[2] = 0.45 / 6.0;
				vQWei[3] = 0.45 / 6.0;
				vQWei[4] = 0.45 / 6.0;

				return; 
			}

			if (num == 8)
			{
				vQPos.resize(8);
				vQWei.resize(8);
				vQPos[0] = Vector3d(0.009437387888358, 0.035220811090087, 0.166666666666667);
				vQPos[1] = Vector3d(0.035220811090087, 0.009437387888358, 0.166666666666667);
				vQPos[2] = Vector3d(0.035220811090087, 0.131445856471988, 0.044658198978444);
				vQPos[3] = Vector3d(0.131445856471988, 0.035220811090087, 0.044658198978444);
				vQPos[4] = Vector3d(0.035220810850163, 0.131445855576580, 0.622008467032738);
				vQPos[5] = Vector3d(0.131445855576580, 0.035220810850163, 0.622008467032738);
				vQPos[6] = Vector3d(0.131445855576580, 0.490562611456158, 0.166666666666667);
				vQPos[7] = Vector3d(0.490562611456158, 0.131445855576580, 0.166666666666667);

				vQWei[0] = 0.001179673492382;
				vQWei[1] = 0.001179673492382;
				vQWei[2] = 0.004402601409914;
				vQWei[3] = 0.004402601409914;
				vQWei[4] = 0.016430731923420;
				vQWei[5] = 0.016430731923420;
				vQWei[6] = 0.061320326343747;
				vQWei[7] = 0.061320326343747;

				return;
			}

			if (num == 11)
			{
				vQPos.resize(11);
				vQWei.resize(11);
				vQPos[0] = Vector3d(0.250000000000000, 0.250000000000000, 0.250000000000000);
				vQPos[1] = Vector3d(0.785714285714286, 0.071428571428571, 0.071428571428571);
				vQPos[2] = Vector3d(0.071428571428571, 0.071428571428571, 0.071428571428571);
				vQPos[3] = Vector3d(0.071428571428571, 0.071428571428571, 0.785714285714286);
				vQPos[4] = Vector3d(0.071428571428571, 0.785714285714286, 0.071428571428571);
				vQPos[5] = Vector3d(0.100596423833201, 0.399403576166799, 0.399403576166799);
				vQPos[6] = Vector3d(0.399403576166799, 0.100596423833201, 0.399403576166799);
				vQPos[7] = Vector3d(0.399403576166799, 0.399403576166799, 0.100596423833201);
				vQPos[8] = Vector3d(0.399403576166799, 0.100596423833201, 0.100596423833201);
				vQPos[9] = Vector3d(0.100596423833201, 0.399403576166799, 0.100596423833201);
				vQPos[10] = Vector3d(0.100596423833201, 0.100596423833201, 0.399403576166799);
				vQWei[0] = -0.0789333333333333 / 6.0;
				vQWei[1] = 0.0457333333333333 / 6.0;
				vQWei[2] = 0.0457333333333333 / 6.0;
				vQWei[3] = 0.0457333333333333 / 6.0;
				vQWei[4] = 0.0457333333333333 / 6.0;
				vQWei[5] = 0.1493333333333333 / 6.0;
				vQWei[6] = 0.1493333333333333 / 6.0;
				vQWei[7] = 0.1493333333333333 / 6.0;
				vQWei[8] = 0.1493333333333333 / 6.0;
				vQWei[9] = 0.1493333333333333 / 6.0;
				vQWei[10] = 0.1493333333333333 / 6.0;

				return;
			}
		}
	}

	bool Cell_Tetra::IsValidParametric(const VectorXd& vp) const
	{
		if (vp.size() != 3)
			return false;

		Real weightSum = 0;
		bool valid = true;

		for (int i = 0; i < vp.size(); ++i)
		{
			if (vp[i] < 0 - 1e-6 || vp[i] > 1 + 1e-6)
			{
				valid = false;
				break; // Found
			}

			weightSum += vp[i];
		}

		if (weightSum > 1 + 1e-6)
			valid = false;

		return valid;
	}

	void Cell_Tetra::SetSubelementPositions(Space s)
	{
		if (this->m_vnodes.size() == 4)
		{
			// Nothing to do here...
		}

		if (this->m_vnodes.size() == 10)
		{
			MatrixXd mN; this->GetNodeMatrix(mN, s);
			this->m_vnodes[4]->DoF()->SetPosition(0.5*(mN.col(0) + mN.col(1)), s);
			this->m_vnodes[5]->DoF()->SetPosition(0.5*(mN.col(1) + mN.col(2)), s);
			this->m_vnodes[6]->DoF()->SetPosition(0.5*(mN.col(0) + mN.col(2)), s);
			this->m_vnodes[7]->DoF()->SetPosition(0.5*(mN.col(0) + mN.col(3)), s);
			this->m_vnodes[8]->DoF()->SetPosition(0.5*(mN.col(1) + mN.col(3)), s);
			this->m_vnodes[9]->DoF()->SetPosition(0.5*(mN.col(2) + mN.col(3)), s);
		}
	}

}