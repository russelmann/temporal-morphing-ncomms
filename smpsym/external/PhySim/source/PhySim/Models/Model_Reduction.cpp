//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Models/Model_Reduction.h>

#include <PhySim/Geometry/Node.h>
#include <PhySim/Energies/EnergyElement.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	Model_Reduction::Model_Reduction() : Model()
	{
		m_pIntModel = NULL;

		this->m_timerComputeReduction = CustomTimer(10, "CAL_REC");
		this->m_timerAssembleReduction = CustomTimer(10, "ASS_REC");
	}

	void Model_Reduction::Init(Model* pModel)
	{
		Model::Init();

		m_pIntModel = pModel;
	}

	Model_Reduction::~Model_Reduction()
	{
		// Nothing to do here...
	}

	void Model_Reduction::Free()
	{
		Model::Free();

		if (this->m_pIntModel != NULL)
			this->m_pIntModel = NULL;
	}

	void Model_Reduction::SetFullDOFPosition(const VectorXd& vx)
	{
		Model::SetFullDOFPosition(vx);
		this->UpdateInternalKinematics();
		this->m_pIntModel->DirtyDeformed();
		this->DirtyDeformed();
	}

	void Model_Reduction::SetFreeDOFPosition(const VectorXd& vx)
	{
		Model::SetFreeDOFPosition(vx);
		this->UpdateInternalKinematics();
		this->m_pIntModel->DirtyDeformed();
		this->DirtyDeformed();
	}

	inline const MatrixSd& Model_Reduction::GetJacobian()
	{
		if (this->IsDirty_Jacobian())
			this->ComputeAndStore_Jacobian();

		return this->m_mJaco.m_msparseMatrix;
	}

	void Model_Reduction::PrepareForSimulation()
	{
		this->m_pIntModel->PrepareForSimulation();

		// Allocate free DOF

		int idCount = 0;

		this->m_numFreeDoF = 0;

		size_t numDoFSet = this->m_vDoFs.size();
		for (size_t i = 0; i < numDoFSet; ++i)
		{
			if (this->m_vDoFs[i]->IsFixed())
			{
				this->m_vDoFs[i]->SetId_Free(-1);
				this->m_vDoFs[i]->SetOffset_Free(-1);
			}
			else
			{
				this->m_vDoFs[i]->SetId_Free(idCount++);
				this->m_vDoFs[i]->SetOffset_Free(m_numFreeDoF);

				this->m_numFreeDoF += this->m_vDoFs[i]->GetNumDim();
			}
		}

		this->m_energy = 0;
		this->m_vgradFree.resize(m_numFreeDoF);
		this->m_mHessFree = FastMatrixSd(m_numFreeDoF, m_numFreeDoF);
		this->m_mMassFree = FastMatrixSd(m_numFreeDoF, m_numFreeDoF);
		this->m_mJaco = FastMatrixSd(this->m_pIntModel->GetNumFreeDOF(), m_numFreeDoF);

		this->ComputeAndStore_Jacobian();

		this->ComputeAndStore_Mass();
		this->ComputeAndStore_Energy();
		this->ComputeAndStore_Gradient();
		this->ComputeAndStore_Hessian();
	}

	void Model_Reduction::DirtyDeformed()
	{
		Model::DirtyDeformed();

		this->m_dirtyFlags = this->m_dirtyFlags | DirtyFlags::Reduction;
	}

	bool Model_Reduction::IsDirty_Energy() const
	{
		return Model::IsDirty_Energy() || this->m_pIntModel->IsDirty_Energy();
	}

	bool Model_Reduction::IsDirty_Gradient(bool full) const
	{
		return Model::IsDirty_Gradient(full) || this->m_pIntModel->IsDirty_Gradient(full);
	}

	bool Model_Reduction::IsDirty_Hessian(bool full) const 
	{
		return Model::IsDirty_Hessian(full) || this->m_pIntModel->IsDirty_Hessian(full);
	}

	bool Model_Reduction::IsDirty_Mass(bool full) const
	{
		return Model::IsDirty_Mass(full) || this->m_pIntModel->IsDirty_Mass(full);
	}

	bool Model_Reduction::IsDirty_Jacobian() const
	{
		return (this->m_dirtyFlags & (int) DirtyFlags::Reduction) != 0;
	}

	void Model_Reduction::ComputeAndStore_Energy()
	{
		this->m_energy = this->m_pIntModel->GetEnergy();
		this->m_dirtyFlags = (DirtyFlags)(this->m_dirtyFlags & ~DirtyFlags::Energy);
	}

	void Model_Reduction::ComputeAndStore_Gradient(bool full)
	{
		if (!full)
		{
			this->m_vgradFree = this->GetJacobian().transpose()*this->m_pIntModel->GetGradient();
			this->m_dirtyFlags = (DirtyFlags)(this->m_dirtyFlags & ~DirtyFlags::GradientFree);
		}
		else
		{
         assert(false);
		}
	}

	void Model_Reduction::ComputeAndStore_Hessian(bool full)
	{
		if (!full)
		{
			MatrixSd mIntH = m_pIntModel->GetHessian().selfadjointView<Lower>();
			this->m_mHessFree.m_msparseMatrix = GetJacobian().transpose()*mIntH*GetJacobian();
			this->m_dirtyFlags = (DirtyFlags)(this->m_dirtyFlags & ~DirtyFlags::HessianFree);

			if (!this->m_mHessFree.HasTripletData())
			{
				this->m_mHessFree.BuildTripletsFromMatrix();
			}

			if (!this->m_mHessFree.HasMappingData())
			{
				this->m_mHessFree.BuildMappingFromMatrix();
			}
		}
		else
		{
			assert(false);
		}
	}

	void Model_Reduction::ComputeAndStore_Mass(bool full)
	{
		if (!full)
		{
			this->m_pIntModel->ComputeAndStore_Mass();
			const MatrixSd& J = this->GetJacobian();
			this->m_mMassFree.m_msparseMatrix = J.transpose()*m_pIntModel->GetMass()*J;
			this->m_dirtyFlags = (DirtyFlags)(this->m_dirtyFlags & ~DirtyFlags::MassFree);
		}
		else
		{
			assert(false);
		}
	}


}
