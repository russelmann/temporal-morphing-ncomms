//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Energies/EnergyElement_FEM.h>

#include <PhySim/Geometry/Polytope.h>
#include <PhySim/Geometry/Node.h>
#include <PhySim/Models/Model_Particles.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	EnergyElement_FEM::EnergyElement_FEM(Model_Particles* pModel, Polytope* pPoly, Material* pMaterial, int numQ) : EnergyElement(pModel, pMaterial)
	{
		this->m_pPoly = pPoly;

		// Check if this is 3D or 2D
		// in the rest configuration

		this->m_order = pPoly->Order();

		int numNode = (int) pPoly->Nodes().size();

		this->m_vDoFs.resize(numNode);
		for (int i = 0; i < numNode; ++i)
			this->m_vDoFs[i] = pPoly->Nodes()[i]->DoF();

		this->m_vgradient.resize(3*numNode);
		this->m_mHessian.resize(3*numNode, 3*numNode);
		this->m_mHessianPFree.resize(3*numNode, 3*numNode);
		this->m_mHessianPFull.resize(3*numNode, 3*numNode);

		this->m_mrestStrain.setIdentity(m_order, m_order);

		m_numQ = numQ;

		this->Init();
	}

	EnergyElement_FEM::~EnergyElement_FEM()
	{
		// Nothing to do...
	}

	void EnergyElement_FEM::Init()
	{
		m_pPoly->GetQuadrature(m_numQ, m_vQpos, m_vQwei);
		int numN = (int) m_pPoly->Nodes().size();

		MatrixXd mT = this->m_mrestStrain.inverse();

		MatrixXd mN0;
		this->m_pPoly->GetNodeMatrix(mN0, Space::MAT);
		if (m_order == 2) // Convert node matrix to 2D
		{
			Vector3d v0 = mN0.col(0);
			Vector3d v1 = mN0.col(1);
			Vector3d v2 = mN0.col(2);
			Vector3d n0 = (v1 - v0);
			Vector3d n1 = (v2 - v0);
			Vector3d n2 = n0.cross(n1);
			n1 = n2.cross(n0);
			n0.normalize();
			n1.normalize();
			n2.normalize();
			Matrix3d T;
			T.col(0) = n0;
			T.col(1) = n1;
			T.col(2) = n2;
			mN0 = (T.inverse()*mN0).block(0, 0, 2, numN);
		}

		// Transform to real rest configuration

		mN0 = mT*mN0;

		// Deformation gradient F 
		//
		// At given material point pi in parametric coordinates, Fi = Nx*Bi*(N0*Bi)^-1, where everything except Nx is cached at initialization
		// depending on the, rest state. For Ai = Bi*(N0*Bi)^-1, and Fi = [F00, F01, ..., F0m, F10, F11, ..., F1m, ..., Fn0, Fn1, ..., Fnm],
		// expressed as a vector, we will compute and store DxDF and implement the different material models as a function of F as a vector. 
		// This way, the gradient dUdx = dUdF*dFdx and the Hessian d2Udx2 = dFdx^T*d2UdF2*dFdx. For Fi = Nx*Ai, the matrix dFdx is constant.
		// For instance, for order = 2, A is a 3x2 matrix and dFdx is a 6x9 matrix and takes the form:
		//
		// A11		0		0		A21		0		0		A31		0		0
		// A12		0		0		A22		0		0		A32		0		0
		// 0		A11		0		0		A21		0		0		A31		0
		// 0		A12		0		0		A22		0		0		A32		0
		// 0		0		A11		0		0		A21		0		0		A31	
		// 0		0		A12		0		0		A22		0		0		A32

		// Energy computation
		//
		// Given some energy density definition depending on the deformation gradient U(F), the total potential is computed as the 
		// integral through the volume in material space: E = Int_X U(F(X)) dX. As shape functions are defined in parametric space,
		// the integral must be defined also in parametric space: E = Int_p U(F(x,B(X))) dVXdVp dp, where the differential change in
		// volume dVXdVp = det(dXdp) = det(N0*B), resulting: E = Int_p U(F(x,B(p))) det(N0*B(p)) dp.
		//
		// For linear basis functions B is constant and the integral is evaluated as just one quadrature point with weight w = Vp, the
		// volume of the iso-parametric element (e.g. 1/6 for the tetrahedron, 8 for the hexahedron, etc.). Note that this results in
		// E = U(F(x,B)) VX, where VX = Vp det(N0*B), the volume of the element in material space.
		//
		// For non-linear basis functions, B(p) is not constant and the integral must be numerically approximated using quadrature as
		// in E = Sum_i U(Fi(x,Bi)) det(N0*Bi) wQi. Independently on the quadrature scheme used, Sum_i wQi = Vp, the volume of the 
		// iso-parametetric element. 
		//
		// TODO: Not clear if similarly to the previous case Sum_i wQi det(N0*Bi) = VX. It should... Check!

		this->m_vmBH0i.resize(m_numQ);
		this->m_vmDFDx.resize(m_numQ);
		for (int i = 0; i < m_numQ; ++i)
		{
			MatrixXd mB;
			this->m_pPoly->ComputeShapeDerivative(m_vQpos[i], mB);
			MatrixXd mH0 = mN0*mB;

			// Precompute and store BH0i = Bi*(N0*Bi)^-1 to compute Fi = Nx*BH0i

			this->m_vmBH0i[i] = mB*mH0.inverse();

			// Precompute integration constant Wi = wQi*|DXDp|_i = wQi*det(N0*Bi)

			this->m_vQwei[i] *= mH0.determinant();

			if (m_order == 2) // It's an element of order 2, consider material thickness
				this->m_vQwei[i] *= (*this->m_pMaterial)[Material::Property::Thickness];

			this->m_vmDFDx[i].resize(3*m_order, 3*numN);
			this->m_vmDFDx[i].setZero();

			for (int ii = 0; ii < 3; ++ii)
				for (int jj = 0; jj < numN; ++jj)
					for (int kk = 0; kk < m_order; ++kk)
						this->m_vmDFDx[i](m_order * ii + kk, 3 * jj + ii) = this->m_vmBH0i[i](jj, kk);
		}

		this->m_intVolume = mT.determinant()*this->m_pPoly->ComputeVolume(Space::MAT);
		if (m_order == 2) // If it's an element of order 2, consider material thickness
			this->m_intVolume *= (*this->m_pMaterial)[Material::Property::Thickness];
	}

	MatrixXd EnergyElement_FEM::ComputeRestStrain()
	{
		// But default the rest strain is the identity, what 
		// means the rest configuration would be equal to the
		// rest configuration.

		return MatrixXd::Identity(this->m_order, this->m_order);
	}

	MatrixXd EnergyElement_FEM::ComputeDeformationGradient() const
	{
		int numQ = (int) this->m_vQpos.size();

		MatrixXd mNx;
		this->m_pPoly->GetNodeMatrix(mNx, Space::DEF);
		MatrixXd mFi = (mNx*m_vmBH0i[0]).transpose();

		for (int i = 1; i < numQ; ++i)
		{
			mFi += (mNx*m_vmBH0i[i]).transpose();
		}

		mFi /= numQ;

		return mFi;
	}

	void EnergyElement_FEM::ComputeAndStore_Energy()
	{
		if (this->m_pPoly->ComputeVolume() < -1)
		{
			this->m_energy = HUGE_VAL;
			return; // Inverted element
		}

		int numQ = (int) this->m_vQpos.size();
		int numN = (int)m_pPoly->Nodes().size();
		this->m_energy = 0;

		MatrixXd mNx;
		this->m_pPoly->GetNodeMatrix(mNx, Space::DEF);

		for (int i = 0; i < numQ; ++i)
		{
			MatrixXd mFi = (mNx*m_vmBH0i[i]).transpose();
			mFi.resize(mFi.size(), 1);

			Real energyi = 0;
			this->ComputeEnergyForF(mFi, energyi);

			//if (energyi > 1e9)
			//	logSimu("\n[WARNING] Infinite FEM energy");

			this->m_energy += energyi*this->m_vQwei[i];
		}	

//		// TEEEEST!
//
//		MatrixXd mNx;
//		this->m_pPoly->GetNodeMatrix(mNx, Space::DEF);
//
//		Real Nx[3][4];
//		for (int i = 0; i < 3; i++)
//			for (int j = 0; j < 4; ++j)
//				Nx[i][j] = mNx(i, j);
//
//		Real BHi[4][3];
//		for (int i = 0; i < 4; i++)
//			for (int j = 0; j < 3; ++j)
//				BHi[i][j] = m_vmBH0i[0](i, j);
//
//		Real lame1 = (*this->m_pMaterial)[Material::Property::Lame1];
//		Real lame2 = (*this->m_pMaterial)[Material::Property::Lame2];
//
//		Real V0 = this->m_intVolume;
//
//#include "../codegen/FEMVol_CoNH_Energy_Test.mcg"
//
//		this->m_energy = t464;
	}

	void EnergyElement_FEM::ComputeAndStore_Gradient()
	{
		int numQ = (int) this->m_vQpos.size();
		int numN = (int)m_pPoly->Nodes().size();
		this->m_vgradient.setZero(3 * numN);

		MatrixXd mNx;
		this->m_pPoly->GetNodeMatrix(mNx, Space::DEF);

		for (int i = 0; i < numQ; ++i)
		{
			MatrixXd mFi = (mNx*m_vmBH0i[i]).transpose();
			mFi.resize(mFi.size(), 1);

			VectorXd vgradienti;
			this->ComputeGradientForF(mFi, vgradienti);
			this->m_vgradient += (this->m_vmDFDx[i].transpose()*vgradienti)*this->m_vQwei[i];
		}

//		// TEEEEST!
//
//		MatrixXd mNx;
//		this->m_pPoly->GetNodeMatrix(mNx, Space::DEF);
//
//		Real Nx[3][4];
//		for (int i = 0; i < 3; i++)
//			for (int j = 0; j < 4; ++j)
//				Nx[i][j] = mNx(i, j);
//
//		Real BHi[4][3];
//		for (int i = 0; i < 4; i++)
//			for (int j = 0; j < 3; ++j)
//				BHi[i][j] = m_vmBH0i[0](i, j);
//
//		Real lame1 = (*this->m_pMaterial)[Material::Property::Lame1];
//		Real lame2 = (*this->m_pMaterial)[Material::Property::Lame2];
//
//		Real V0 = this->m_intVolume;
//
//		Real vgx[12];
//
//#include "../codegen/FEMVol_CoNH_Gradient_Test.mcg"
//
//		for (int i = 0; i < 12; ++i)
//			m_vgradient(i) = vgx[i]; // Copy
	}

	void EnergyElement_FEM::ComputeAndStore_Hessian()
	{
		int numQ = (int) this->m_vQpos.size();
		int numN = (int)m_pPoly->Nodes().size();
		m_mHessian.setZero(3 * numN, 3 * numN);

		MatrixXd mNx;
		this->m_pPoly->GetNodeMatrix(mNx, Space::DEF);

		for (int i = 0; i < numQ; ++i)
		{
			MatrixXd mFi = (mNx*m_vmBH0i[i]).transpose();
			mFi.resize(mFi.size(), 1);

			MatrixXd mHessiani;
			this->ComputeHessianForF(mFi, mHessiani);
			this->m_mHessian += (this->m_vmDFDx[i].transpose()*mHessiani*this->m_vmDFDx[i])*this->m_vQwei[i];
		}

//		// TEEEEST!
//
//		MatrixXd mNx;
//		this->m_pPoly->GetNodeMatrix(mNx, Space::DEF);
//
//		Real Nx[3][4];
//		for (int i = 0; i < 3; i++)
//			for (int j = 0; j < 4; ++j)
//				Nx[i][j] = mNx(i, j);
//
//		Real BHi[4][3];
//		for (int i = 0; i < 4; i++)
//			for (int j = 0; j < 3; ++j)
//				BHi[i][j] = m_vmBH0i[0](i, j);
//
//		Real lame1 = (*this->m_pMaterial)[Material::Property::Lame1];
//		Real lame2 = (*this->m_pMaterial)[Material::Property::Lame2];
//
//		Real V0 = this->m_intVolume;
//
//		Real mHx[12][12];
//
//#include "../codegen/FEMVol_CoNH_Hessian_Test.mcg"
//
//		for (int i = 0; i < 12; ++i)
//			for (int j = 0; j < 12; ++j)
//				m_mHessian(i, j) = mHx[i][j]; // Copy
	}

}