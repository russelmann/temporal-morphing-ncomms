//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================


#include <PhySim/Geometry/RigidBody.h>
#include <PhySim/Geometry/Node.h>

#include <PhySim/Utils/eulerRotation.h>
#include <PhySim/Utils/rigidBodyTransformEuler.h>
#include <PhySim/Utils/rigidBodyTransformQuat.h>
#include <PhySim/Utils/getRBNodeJacobianEuler.h>
#include <PhySim/Utils/getRBNodeJacobianQuat.h>
#include <PhySim/Utils/quaternionNormalized.h>
#include <PhySim/Utils/getQuatRotNodeJacobian.h>
#include <PhySim/Utils/getQuatRotXCoordHessian.h>
#include <PhySim/Utils/getQuatRotYCoordHessian.h>
#include <PhySim/Utils/getQuatRotZCoordHessian.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	RigidBody::RigidBody(Polytope* pPoly)
	{
		assert(pPoly != NULL);

		this->pPolytope = pPoly;
		this->m_vIntDoFs.resize(2);
		this->m_vIntDoFs[0] = new DoFSet(3);
		this->m_vIntDoFs[1] = new DoFSet(3);

		this->Init();
		this->UpdatePositions();
		this->UpdateJacobian();
	}

	RigidBody::~RigidBody()
	{
		delete this->pPolytope;
		delete this->m_vIntDoFs[0];
		delete this->m_vIntDoFs[1];
		this->m_vIntDoFs.clear();
	}

	void RigidBody::Init()
	{
		this->QuaterDoF()->SetPositions(Vector3d::Zero());
		this->CenterDoF()->SetPositions(Vector3d::Zero());

		// Initialize local coordinates

		// Center

		int N = (int) this->pPolytope->Nodes().size();
		Vector3d centroid = pPolytope->ComputeCentroid();
		this->CenterDoF()->SetPositions(centroid);

		// Rotation

		this->QuaterDoF()->SetPositions(Vector3d::Zero());

		this->m_mR.setIdentity();
		this->m_mT.resize(N, 3);

		this->m_mJValues = MatrixXd::Zero(3*N,6);
		this->m_mJPointers = MatrixXp::Zero(3*N,6);
		this->m_vmHValues.resize(3 * N);
		this->m_vmHPointers.resize(3 * N);
		for (int i = 0; i < 3 * N; ++i)
		{
			this->m_vmHValues[i] = MatrixXd::Zero(3 * N, 6);
			this->m_vmHPointers[i] = MatrixXp::Zero(3 * N, 6);
		}

		this->m_vOutDoFs.resize(N);

		for (int i = 0; i < N; ++i)
		{
			this->m_vOutDoFs[i] = this->pPolytope->Nodes()[i]->DoF();
			this->m_mT.row(i) = m_vOutDoFs[i]->GetPosition_0() - centroid;
		}
	}

	Vector3d RigidBody::GetCentroid() 
	{
		return this->CenterDoF()->GetPosition_x();
	}

	Matrix3d RigidBody::GetRotation() 
	{
		const Vector3d& rot = this->QuaterDoF()->GetPosition();

		// Rotate identity columns

		if (rot.isZero(1e-9))
		{
			return this->m_mR;
		}
		else
		{
			Matrix3d mZ = Matrix3d::Zero();
			Matrix3d mRT;
			rigidBodyTransformQuat(this->m_mR.col(0).data(), mZ.col(0).data(), rot.data(), mRT.col(0).data());
			rigidBodyTransformQuat(this->m_mR.col(1).data(), mZ.col(1).data(), rot.data(), mRT.col(1).data());
			rigidBodyTransformQuat(this->m_mR.col(2).data(), mZ.col(2).data(), rot.data(), mRT.col(2).data());

			return mRT;
		}
	}

	void RigidBody::RecomputeCentroid()
	{
		this->CenterDoF()->SetPosition_x(Poly()->ComputeCentroid(PhySim::DEF));
	}

	void RigidBody::RecomputeRotation()
	{
		Eigen::MatrixXd mVX;
		Eigen::MatrixXd mV0;
		this->Poly()->GetNodeMatrix(mV0, PhySim::Space::MAT);
		this->Poly()->GetNodeMatrix(mVX, PhySim::Space::DEF);

		Eigen::VectorXd vt;
		Eigen::MatrixXd mR;
		PhySim::computeBestRigidTransform(mV0.transpose(), 
										  mVX.transpose(), 
										  mR, vt);
		this->QuaterDoF()->SetPosition_x(VectorXd::Zero(3));
		this->m_mR = mR;
	}

	void RigidBody::UpdatePositions(Space s)
	{
		const Vector3d& rot0 = this->QuaterDoF()->GetPosition(s);

		// Is rotation over a threshold?

		if (rot0.norm() > 1e-3)
		{
			// Update accumulated rotation

			this->m_mR = this->GetRotation();

			// Reset current rotation to zero

			this->QuaterDoF()->SetPosition(Vector3d::Zero(), s);
		}

		// Update node positions

		int N = (int) this->pPolytope->Nodes().size();

		const Vector3d& rot = this->QuaterDoF()->GetPosition(s);
		const Vector3d& cen = this->CenterDoF()->GetPosition(s);

		for (int i = 0; i < N; ++i)
		{
			Vector3d vec = this->m_mR*this->m_mT.row(i).transpose();
			
			Vector3d pos;

			if (rot.isZero(1e-9))
			{
				pos = cen + vec;
			}
			else
			{
				rigidBodyTransformQuat(vec.data(),
									   cen.data(),
									   rot.data(),
									   pos.data());
			}

			this->m_vOutDoFs[i]->SetPosition(pos, s);
		}
	}

	void RigidBody::UpdateJacobian(Space s)
	{
		int N = (int) this->pPolytope->Nodes().size();

		const Vector3d& rot = this->QuaterDoF()->GetPosition(s);
		const Vector3d& cen = this->CenterDoF()->GetPosition(s);

		Vector3d rotC;
		if (rot.isZero(1e-9)) // Correct me!
			rotC = Vector3d(0.0, 1e-9, 0.0);
		else rotC = rot; 

		// Update node Jacobian

		for (int i = 0; i < N; ++i)
		{
			Vector3d vec = this->m_mR*this->m_mT.row(i).transpose();

			MatrixXd mDnDr(3, 3);
			getQuatRotNodeJacobian(vec.data(), rotC.data(), mDnDr.data());

			m_mJValues.block(3 * i, 0, 3, 3) = Matrix3d::Identity();
			m_mJValues.block(3 * i, 3, 3, 3) = mDnDr.transpose();
		}
	}

	void RigidBody::UpdateHessian(Space s)
	{
		int N = (int) this->pPolytope->Nodes().size();

		const Vector3d& rot = this->QuaterDoF()->GetPosition(s);
		const Vector3d& cen = this->CenterDoF()->GetPosition(s);

		Vector3d rotC;
		if (rot.isZero(1e-9)) // Correct me!
			rotC = Vector3d(0.0, 1e-9, 0.0);
		else rotC = rot;

		// Update coordinates Hessian

		for (int i = 0; i < N; ++i)
		{
			Vector3d vec = this->m_mR*this->m_mT.row(i).transpose();

			MatrixXd mH(3,3);

			m_vmHValues[3 * i + 0].setZero();
			getQuatRotXCoordHessian(vec.data(), rotC.data(), mH.data());
			m_vmHValues[3 * i + 0].block(3, 3, 3, 3) = mH.transpose();

			m_vmHValues[3 * i + 1].setZero();
			getQuatRotYCoordHessian(vec.data(), rotC.data(), mH.data());
			m_vmHValues[3 * i + 1].block(3, 3, 3, 3) = mH.transpose();

			m_vmHValues[3 * i + 2].setZero();
			getQuatRotZCoordHessian(vec.data(), rotC.data(), mH.data());
			m_vmHValues[3 * i + 2].block(3, 3, 3, 3) = mH.transpose();
		}
	}

	void RigidBody::AssembleGlobal_Jacobian(VectorTd& vtotalJ)
	{
		int localOffset_i = 0;
		int localOffset_j = 0;

		size_t numDoFSetOut = this->m_vOutDoFs.size();
		size_t numDoFSetInt = this->m_vIntDoFs.size();
		for (size_t i = 0; i < numDoFSetOut; ++i)
		{
			int freeId_i = this->m_vOutDoFs[i]->GetId_Free();
			if (freeId_i < 0)
			{
				localOffset_i += this->m_vOutDoFs[i]->GetNumDim();
				continue;
			}

			for (size_t j = 0; j < numDoFSetInt; ++j)
			{
				int freeId_j = this->m_vIntDoFs[j]->GetId_Free();
				if (freeId_j < 0)
				{
					localOffset_j += this->m_vIntDoFs[j]->GetNumDim();
					continue;
				}

				assembleSparseMatrix(vtotalJ,
					this->m_vOutDoFs[i]->GetOffset_Free(),
					this->m_vIntDoFs[j]->GetOffset_Free(),
					m_mJValues.block(localOffset_i, localOffset_j, this->m_vOutDoFs[i]->GetNumDim(), this->m_vIntDoFs[j]->GetNumDim()));

				localOffset_j += this->m_vIntDoFs[j]->GetNumDim();
			}

			localOffset_j = 0;
			localOffset_i += this->m_vOutDoFs[i]->GetNumDim();
		}
	}

	void RigidBody::AssembleGlobal_FastPreallocatedJacobian()
	{
		int localOffset_i = 0;
		int localOffset_j = 0;

		size_t numDoFSetOut = this->m_vOutDoFs.size();
		size_t numDoFSetInt = this->m_vIntDoFs.size();

		for (size_t i = 0; i < numDoFSetOut; ++i)
		{
			int freeId_i = this->m_vOutDoFs[i]->GetId_Free();
			if (freeId_i < 0)
			{
				localOffset_i += this->m_vOutDoFs[i]->GetNumDim();
				continue;
			}

			for (size_t j = 0; j < numDoFSetInt; ++j)
			{
				int freeId_j = this->m_vIntDoFs[j]->GetId_Free();
				if (freeId_j < 0)
				{
					localOffset_j += this->m_vIntDoFs[j]->GetNumDim();
					continue;
				}

				for (int ii = 0; ii < this->m_vOutDoFs[i]->GetNumDim(); ++ii)
				{
					for (int jj = 0; jj < this->m_vIntDoFs[j]->GetNumDim(); ++jj)
					{
						if (m_mJPointers(localOffset_i + ii, localOffset_j + jj) != NULL)
						{
							*m_mJPointers(localOffset_i + ii, localOffset_j + jj) += m_mJValues(localOffset_i + ii, localOffset_j + jj);
						}
					}
				}

				localOffset_j += this->m_vIntDoFs[j]->GetNumDim();
			}

			localOffset_j = 0;
			localOffset_i += this->m_vOutDoFs[i]->GetNumDim();
		}
	}

	void RigidBody::AllocateGlobal_Jacobian(const CoefMap& mp)
	{
		int localOffset_i = 0;
		int localOffset_j = 0;

		this->m_mJPointers.setZero();

		size_t numDoFSetOut = this->m_vOutDoFs.size();
		size_t numDoFSetInt = this->m_vIntDoFs.size();
		for (size_t i = 0; i < numDoFSetOut; ++i)
		{
			int freeId_i = this->m_vOutDoFs[i]->GetId_Free();
			if (freeId_i < 0)
			{
				localOffset_i += this->m_vOutDoFs[i]->GetNumDim();
				continue;
			}

			for (size_t j = 0; j < numDoFSetInt; ++j)
			{
				int freeId_j = this->m_vIntDoFs[j]->GetId_Free();
				if (freeId_j < 0)
				{
					localOffset_j += this->m_vIntDoFs[j]->GetNumDim();
					continue;
				}

				for (int ii = 0; ii < this->m_vOutDoFs[i]->GetNumDim(); ++ii)
				{
					for (int jj = 0; jj < this->m_vIntDoFs[j]->GetNumDim(); ++jj)
					{
						m_mJPointers(localOffset_i + ii, localOffset_j + jj) = mp.at(IntPair(this->m_vOutDoFs[i]->GetOffset_Free() + ii, this->m_vIntDoFs[j]->GetOffset_Free() + jj));
					}
				}

				localOffset_j += this->m_vIntDoFs[j]->GetNumDim();
			}

			localOffset_j = 0;
			localOffset_i += this->m_vOutDoFs[i]->GetNumDim();
		}
	}

	void RigidBody::AssembleGlobal_VectorHessianProduct(const VectorXd& vv, VectorTd& vtotalJ)
	{
		int localOffset_i = 0;

		size_t numDoFSetOut = this->m_vOutDoFs.size();
		size_t numDoFSetInt = this->m_vIntDoFs.size();
		for (size_t i = 0; i < numDoFSetOut; ++i)
		{
			int freeId_i = this->m_vOutDoFs[i]->GetId_Free();
			if (freeId_i < 0)
			{
				localOffset_i += this->m_vOutDoFs[i]->GetNumDim();
				continue;
			}

			for (size_t ii = 0; ii < this->m_vOutDoFs[i]->GetNumDim(); ++ii)
			{
				Real value = vv(this->m_vOutDoFs[i]->GetOffset_Free() + ii);

				const MatrixXd& Mvi = m_vmHValues[localOffset_i + ii];

				int localOffset_j = 0;
				int localOffset_k = 0;

				for (size_t j = 0; j < numDoFSetInt; ++j)
				{
					int freeId_j = this->m_vIntDoFs[j]->GetId_Free();
					if (freeId_j < 0)
					{
						localOffset_j += this->m_vIntDoFs[j]->GetNumDim();
						continue;
					}
					for (size_t k = 0; k < numDoFSetInt; ++k)
					{
						int freeId_k = this->m_vIntDoFs[k]->GetId_Free();
						if (freeId_k < 0)
						{
							localOffset_k += this->m_vIntDoFs[k]->GetNumDim();
							continue;
						}

						assembleSparseMatrix(vtotalJ,
							this->m_vIntDoFs[j]->GetOffset_Free(),
							this->m_vIntDoFs[k]->GetOffset_Free(),
							value*Mvi.block(localOffset_j, localOffset_k, this->m_vIntDoFs[j]->GetNumDim(), this->m_vIntDoFs[k]->GetNumDim()));

						localOffset_k += this->m_vIntDoFs[k]->GetNumDim();
					}

					localOffset_k = 0;
					localOffset_j += this->m_vIntDoFs[j]->GetNumDim();
				}
			}

			localOffset_i += this->m_vOutDoFs[i]->GetNumDim();
		}
	}

	void RigidBody::AssembleGlobal_FastPreallocatedVectorHessianProduct(const VectorXd& vv)
	{
		int localOffset_i = 0;

		size_t numDoFSetOut = this->m_vOutDoFs.size();
		size_t numDoFSetInt = this->m_vIntDoFs.size();
		for (size_t i = 0; i < numDoFSetOut; ++i)
		{
			int freeId_i = this->m_vOutDoFs[i]->GetId_Free();
			if (freeId_i < 0)
			{
				localOffset_i += this->m_vOutDoFs[i]->GetNumDim();
				continue;
			}

			for (size_t ii = 0; ii < this->m_vOutDoFs[i]->GetNumDim(); ++ii)
			{
				Real value = vv(this->m_vOutDoFs[i]->GetOffset_Free() + ii);

				const MatrixXd& Mvi = m_vmHValues[localOffset_i + ii];
				const MatrixXp& Mpi = m_vmHPointers[localOffset_i + ii];

				int localOffset_j = 0;
				int localOffset_k = 0;

				for (size_t j = 0; j < numDoFSetInt; ++j)
				{
					int freeId_j = this->m_vIntDoFs[j]->GetId_Free();
					if (freeId_j < 0)
					{
						localOffset_j += this->m_vIntDoFs[j]->GetNumDim();
						continue;
					}
					for (size_t k = 0; k < numDoFSetInt; ++k)
					{
						int freeId_k = this->m_vIntDoFs[k]->GetId_Free();
						if (freeId_k < 0)
						{
							localOffset_k += this->m_vIntDoFs[k]->GetNumDim();
							continue;
						}

						for (int jj = 0; jj < this->m_vIntDoFs[j]->GetNumDim(); ++jj)
						{
							for (int kk = 0; kk < this->m_vIntDoFs[k]->GetNumDim(); ++kk)
							{
								if (Mpi(localOffset_j + jj, localOffset_k + kk) != NULL)
								{
									*Mpi(localOffset_j + jj, localOffset_k + kk) += value*Mvi(localOffset_j + jj, localOffset_k + kk);
								}
							}
						}

						localOffset_k += this->m_vIntDoFs[k]->GetNumDim();
					}

					localOffset_k = 0;
					localOffset_j += this->m_vIntDoFs[j]->GetNumDim();
				}
			}

			localOffset_i += this->m_vOutDoFs[i]->GetNumDim();
		}
	}

	void RigidBody::AllocateGlobal_VectorHessianProduct(const CoefMap& mp)
	{
		size_t numDoFSetInt = this->m_vIntDoFs.size();

		int localOffset_j = 0;
		int localOffset_k = 0;

		for (size_t j = 0; j < numDoFSetInt; ++j)
		{
			int freeId_j = this->m_vIntDoFs[j]->GetId_Free();
			if (freeId_j < 0)
			{
				localOffset_j += this->m_vIntDoFs[j]->GetNumDim();
				continue;
			}
			for (size_t k = 0; k < numDoFSetInt; ++k)
			{
				int freeId_k = this->m_vIntDoFs[k]->GetId_Free();
				if (freeId_k < 0)
				{
					localOffset_k += this->m_vIntDoFs[k]->GetNumDim();
					continue;
				}

				for (int jj = 0; jj < this->m_vIntDoFs[j]->GetNumDim(); ++jj)
				{
					for (int kk = 0; kk < this->m_vIntDoFs[k]->GetNumDim(); ++kk)
					{
						for (int i = 0; i < (int) this->m_vmHValues.size(); ++i)
							m_vmHPointers[i](localOffset_j + jj, localOffset_k + kk) = mp.at(IntPair(this->m_vIntDoFs[j]->GetOffset_Free() + jj, this->m_vIntDoFs[k]->GetOffset_Free() + kk));
					}
				}

				localOffset_k += this->m_vIntDoFs[k]->GetNumDim();
			}

			localOffset_k = 0;
			localOffset_j += this->m_vIntDoFs[j]->GetNumDim();
		}
	}
}