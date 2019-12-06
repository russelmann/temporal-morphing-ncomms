//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Energies/EnergyElement_SpringLinear.h>

#include <PhySim/Geometry/Node.h>
#include <PhySim/Geometry/Edge.h>
#include <PhySim/Models/Model_Particles.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	EnergyElement_SpringLinear::EnergyElement_SpringLinear(Model_Particles* pModel, Edge* pEdge, Material* pMaterial) : EnergyElement(pModel, pMaterial)
	{
		this->m_pEdge = pEdge;

		this->m_vDoFs.resize(2);
		this->m_vDoFs[0] = pEdge->GetOrigin()->DoF();
		this->m_vDoFs[1] = pEdge->GetHead()->DoF();

		this->m_vgradient.resize(3);
		this->m_mHessian.resize(3,3);
		this->m_mHessianPFree.resize(6, 6);
		this->m_mHessianPFull.resize(6, 6);

		this->Init();
	}

	EnergyElement_SpringLinear::~EnergyElement_SpringLinear()
	{
		// Nothing to do...
	}

	void EnergyElement_SpringLinear::Init()
	{
		this->m_intVolume = this->m_restLength = this->ComputeRestLength();
	}

	Real EnergyElement_SpringLinear::ComputeRestLength()
	{
		return this->m_pEdge->ComputeVolume();
	}

	void EnergyElement_SpringLinear::ComputeAndStore_Energy()
	{
		// Get defo nodes

		const Vector3d& ax = this->m_vDoFs[0]->GetPosition_x();
		const Vector3d& bx = this->m_vDoFs[1]->GetPosition_x();

		// Strain

		double s = ((bx - ax).norm() / this->m_restLength - 1);

		// Stiffness

		double k = (*this->m_pMaterial)[Material::Property::StretchK];

		this->m_energy = 0.5*k*this->m_restLength*s*s;
	}

	void EnergyElement_SpringLinear::ComputeAndStore_Gradient()
	{
		// Get defo nodes

		const Vector3d& ax = this->m_vDoFs[0]->GetPosition_x();
		const Vector3d& bx = this->m_vDoFs[1]->GetPosition_x();

		// Compute lengths

		double Lx = (bx - ax).norm();

		// Strain

		double s = (Lx / this->m_restLength - 1);

		// Stiffness

		double k = (*this->m_pMaterial)[Material::Property::StretchK];

		// GradientFull

		Vector3d ub = (bx - ax) / Lx;
		Vector3d gb = k*s*ub;
		this->m_vgradient.block<3, 1>(0, 0) = gb;
	}

	void EnergyElement_SpringLinear::ComputeAndStore_Hessian()
	{
		// Get defo nodes

		const Vector3d& ax = this->m_vDoFs[0]->GetPosition_x();
		const Vector3d& bx = this->m_vDoFs[1]->GetPosition_x();

		// Compute lengths

		double L0 = this->m_restLength;
		double Lx = (bx - ax).norm();

		// Strain

		double s = (Lx / L0 - 1);

		// Stiffness

		double k = (*this->m_pMaterial)[Material::Property::StretchK];

		// HessianFull

		Vector3d ub = (bx - ax) / Lx;
		double L0inv = 1.0 / L0;
		double Lxinv = 1.0 / Lx;
		Matrix3d mI = Matrix3d::Identity();
		Matrix3d ububT = ub*ub.transpose();
		this->m_mHessian.block<3, 3>(0, 0) = k*(L0inv*ububT + s*((mI - ububT)*Lxinv));
	}

	void EnergyElement_SpringLinear::AssembleGlobal_Gradient(VectorXd& vtotalGradient, bool full)
	{
		if (!full)
		{
			if (!this->m_vDoFs[0]->IsFixed())
			{
				int freeId = this->m_vDoFs[0]->GetId_Free();
				vtotalGradient.block<3, 1>(3 * freeId, 0) += -this->m_vgradient;
			}

			if (!this->m_vDoFs[1]->IsFixed())
			{
				int freeId = this->m_vDoFs[1]->GetId_Free();
				vtotalGradient.block<3, 1>(3 * freeId, 0) += this->m_vgradient;
			}
		}
		else
		{
         assert(false);
		}
	}

	void EnergyElement_SpringLinear::AssembleGlobal_Hessian(VectorTd& vtotalHessian, bool full)
	{
		if (!full)
		{
			int freeId_0 = this->m_vDoFs[0]->GetId_Free();
			int freeId_1 = this->m_vDoFs[1]->GetId_Free();
			int offset0 = 3 * freeId_0;
			int offset1 = 3 * freeId_1;

			if (freeId_0 != -1)
			{
				assembleSparseMatrix(vtotalHessian, offset0, offset0, m_mHessian.block<3, 3>(0, 0));
			}

			if (freeId_1 != -1)
			{
				assembleSparseMatrix(vtotalHessian, offset1, offset1, m_mHessian.block<3, 3>(0, 0));
			}

			if (freeId_0 != -1 && freeId_1 != -1)
			{
				if (freeId_0 > freeId_1)
				{
					assembleSparseMatrix(vtotalHessian, offset0, offset1, -m_mHessian.block<3, 3>(0, 0));
				}
				else
				{
					assembleSparseMatrix(vtotalHessian, offset1, offset0, -m_mHessian.block<3, 3>(0, 0));
				}
			}
		}
		else
		{
         assert(false);
      }
	}

	void EnergyElement_SpringLinear::AssembleGlobal_FastPreallocatedHessian(bool full)
	{
		if (!full)
		{
			int freeId_0 = this->m_vDoFs[0]->GetId_Free();
			int freeId_1 = this->m_vDoFs[1]->GetId_Free();
			int offset0 = 3 * freeId_0;
			int offset1 = 3 * freeId_1;

			if (freeId_0 != -1)
			{
				for (int i = 0; i < 3; ++i)
				{
					for (int j = 0; j < 3; ++j)
					{
						*m_mHessianPFree(i, j) += m_mHessian(i, j);
					}
				}
			}

			if (freeId_1 != -1)
			{
				for (int i = 0; i < 3; ++i)
				{
					for (int j = 0; j < 3; ++j)
					{
						*m_mHessianPFree(3 + i, 3 + j) += m_mHessian(i, j);
					}
				}
			}

			if (freeId_0 != -1 && freeId_1 != -1)
			{
				if (freeId_0 > freeId_1)
				{
					for (int i = 0; i < 3; ++i)
					{
						for (int j = 0; j < 3; ++j)
						{
							*m_mHessianPFree(i, 3 + j) += -m_mHessian(i, j);
						}
					}
				}
				else
				{
					for (int i = 0; i < 3; ++i)
					{
						for (int j = 0; j < 3; ++j)
						{
							*m_mHessianPFree(3 + i, j) += -m_mHessian(i, j);
						}
					}
				}
			}
		}
		else
		{
         assert(false);
      }
	}

}
