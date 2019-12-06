//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Models/Model_Reduction_RBCloud.h>

#include <PhySim/Geometry/Node.h>
#include <PhySim/Geometry/Polytope.h>
#include <PhySim/Geometry/RigidBody.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	Model_Reduction_RBCloud::Model_Reduction_RBCloud()
	{
		// Nothing to do here...
	}

	Model_Reduction_RBCloud::~Model_Reduction_RBCloud()
	{
		// Nothing to do here...
	}

	void Model_Reduction_RBCloud::Free()
	{
		Model_Reduction::Free();

		int numRB = (int) this->m_vBody.size();
		for (int i = 0; i < numRB; ++i)
			delete this->m_vBody[i];
		this->m_vBody.clear();
	}

	IModel::StateP Model_Reduction_RBCloud::CreateState(Space s) const
	{
		StateP pS = StateP(new State());
		this->GetState(pS, s);
		return pS;
	}

	bool Model_Reduction_RBCloud::HasState(IModel::StateP pS, Space s) const
	{
		if (!Model_Reduction::HasState(pS, s))
			return false;

		State* pSInt = static_cast<State*>(pS.get());

		size_t numRB = this->m_vBody.size();

		for (size_t i = 0; i < numRB; ++i)
			if (pSInt->m_vmR[i] != this->m_vBody[i]->R())
				return false;

		return true;
	}

	void Model_Reduction_RBCloud::GetState(IModel::StateP pS, Space s) const
	{
		Model_Reduction::GetState(pS, s);

		State* pSInt = static_cast<State*>(pS.get());

		size_t numRB = this->m_vBody.size();

		pSInt->m_vmR.resize(numRB);
		for (size_t i = 0; i < numRB; ++i)
			pSInt->m_vmR[i] = this->m_vBody[i]->R();
	}

	void Model_Reduction_RBCloud::SetState(const IModel::StateP pS, Space s)
	{
		Model_Reduction::SetState(pS, s);

		const State* pSInt = static_cast<const State*>(pS.get());

		size_t numRB = this->m_vBody.size();
		assert(pSInt->m_vmR.size() == numRB);

		for (size_t i = 0; i < numRB; ++i)
			this->m_vBody[i]->R() = pSInt->m_vmR[i];
	}

	void Model_Reduction_RBCloud::Init(Model* pModel, const vector<Polytope*>& vPoly, Type type)
	{
		Model_Reduction::Init(pModel);

		this->m_type = type;

		const vector<DoFSet*> vIntDoF = this->m_pIntModel->GetDoFSets();
		int numIntDoFSet = (int) vIntDoF.size();
		bVector vrbDoFStencil(numIntDoFSet,false);

		this->m_numFullDoF = 0;
		this->m_numRigidDoF = 0;
		this->m_numOtherDoF = 0;

		// Add rigidbody DoF

		int numRB = (int)vPoly.size();
		this->m_vBody.resize(numRB);
		this->m_vDoFs.reserve(2*numRB);
		for (int i = 0; i < numRB; ++i)
		{
			for (size_t j = 0; j < vPoly[i]->Nodes().size(); ++j) // Add stencil
				vrbDoFStencil[vPoly[i]->Nodes()[j]->DoF()->GetId_Full()] = true;

			this->m_vBody[i] = new RigidBody(vPoly[i]);
			this->m_vBody[i]->CenterDoF()->SetId_Full(2 * i + 0);
			this->m_vBody[i]->QuaterDoF()->SetId_Full(2 * i + 1);
			this->m_vBody[i]->CenterDoF()->SetOffset_Full(this->m_numFullDoF + 0);
			this->m_vBody[i]->QuaterDoF()->SetOffset_Full(this->m_numFullDoF + 3);
			this->m_vDoFs.push_back(this->m_vBody[i]->CenterDoF());
			this->m_vDoFs.push_back(this->m_vBody[i]->QuaterDoF());
			this->m_numFullDoF += 6;
			this->m_numRigidDoF += 6;
		}

		// Add non-rigidbody DoF sets

		int countDoFSet = 0;

		for (int i = 0; i < numIntDoFSet; ++i)
		{
			if (vrbDoFStencil[i])
				continue; // Added

			DoFSet* pMappedDoF = new DoFSet(vIntDoF[i]->GetNumDim());
			pMappedDoF->SetPosition_x(vIntDoF[i]->GetPosition_x());
			pMappedDoF->SetPosition_0(vIntDoF[i]->GetPosition_0());
			pMappedDoF->SetOffset_Full(this->m_numFullDoF);
			pMappedDoF->SetId_Full(2*numRB + countDoFSet++);
			this->m_vDoFs.push_back(pMappedDoF);

			this->m_numFullDoF += pMappedDoF->GetNumDim();
			this->m_numOtherDoF += pMappedDoF->GetNumDim();
	
			this->m_mmapDoF.insert(pair<DoFSet*, DoFSet*>(vIntDoF[i], pMappedDoF));
		}

		this->DirtyUndeformed();

		logSimu("\n--");
		logSimu("\nINITIALIZED %s", this->GetName().c_str());
		logSimu("\nTotal internal DOF: %i", m_pIntModel->GetNumFullDOF());
		logSimu("\nTotal external DOF: %i", this->GetNumFullDOF());
		logSimu("\n--");
	}

	void Model_Reduction_RBCloud::UpdateInternalKinematics()
	{
		int numRB = (int) this->m_vBody.size();

		for (int i = 0; i < numRB; ++i)
			this->m_vBody[i]->UpdatePositions();

		map<DoFSet*, DoFSet*>::iterator itCur = this->m_mmapDoF.begin();
		map<DoFSet*, DoFSet*>::iterator itEnd = this->m_mmapDoF.end();
		for (; itCur != itEnd; ++itCur)
		{
			itCur->first->SetPosition(itCur->second->GetPosition(Space::DEF), Space::DEF);
		}

		this->m_pIntModel->DirtyDeformed();
	}

	void Model_Reduction_RBCloud::ComputeAndStore_Jacobian()
	{
		int numRB = (int) this->m_vBody.size();

		// Compute

		if (m_isProfiling) this->m_timerComputeReduction.start();
#pragma omp parallel for
		for (int i = 0; i < numRB; ++i)
		{
			this->m_vBody[i]->UpdateJacobian();
			this->m_vBody[i]->UpdateHessian();
		}
		if (m_isProfiling) this->m_timerComputeReduction.stopStoreLog();

		// Assemble

		if (m_isProfiling) this->m_timerAssembleReduction.start();
		if (m_fastAssembly && !this->m_mJaco.HasMappingData())
		{
			// Gather triplets

			this->m_mJaco.m_vvalueTriplets.clear();

			for (size_t i = 0; i < numRB; ++i)
				this->m_vBody[i]->AssembleGlobal_Jacobian(this->m_mJaco.m_vvalueTriplets);

			map<DoFSet*, DoFSet*>::iterator itCur = this->m_mmapDoF.begin();
			map<DoFSet*, DoFSet*>::iterator itEnd = this->m_mmapDoF.end();
			for (itCur = this->m_mmapDoF.begin(); itCur != itEnd; ++itCur)
			{
				DoFSet* intDoF = itCur->first;
				DoFSet* extDoF = itCur->second;
				int numDim = intDoF->GetNumDim();
				for (int j = 0; j < numDim; ++j)
					this->m_mJaco.m_vvalueTriplets.push_back(Triplet<Real>(intDoF->GetOffset_Free() + j, extDoF->GetOffset_Free() + j, 1.0));
			}

			// Allocate mapping

			this->m_mJaco.BuildMatrixFromTriplets();
			this->m_mJaco.BuildMappingFromMatrix();

			for (size_t i = 0; i < numRB; ++i)
				this->m_vBody[i]->AllocateGlobal_Jacobian(this->m_mJaco.m_coeffMap);

			for (itCur = this->m_mmapDoF.begin(); itCur != itEnd; ++itCur)
			{
				DoFSet* intDoF = itCur->first;
				DoFSet* extDoF = itCur->second;
				int numDim = intDoF->GetNumDim();
				for (int j = 0; j < numDim; ++j)
					this->m_vmapJac.push_back(this->m_mJaco.m_coeffMap.at(IntPair(intDoF->GetOffset_Free() + j, extDoF->GetOffset_Free() + j)));
			}
		}
		else
		{
			this->m_mJaco.m_msparseMatrix *= 0.0;

			for (size_t i = 0; i < numRB; ++i)
				this->m_vBody[i]->AssembleGlobal_FastPreallocatedJacobian();

			for (int i = 0; i < (int) m_vmapJac.size(); ++i) 
				*m_vmapJac[i] = 1.0; // Non-reduced DoFs I
		}
		if (m_isProfiling) this->m_timerAssembleReduction.stopStoreLog();

		// Clean

		this->m_dirtyFlags = (DirtyFlags)(this->m_dirtyFlags & ~DirtyFlags::Reduction);
	}

	void Model_Reduction_RBCloud::ComputeAndStore_Hessian(bool full)
	{
		int numRB = (int) this->m_vBody.size();

		if (!full)
		{
			if (!this->m_mHessFree.HasTripletData())
			{
				MatrixSd mIntH = this->m_pIntModel->GetHessian().selfadjointView<Lower>();
				m_mHessFree.m_msparseMatrix = GetJacobian().transpose()*mIntH*GetJacobian();

				this->m_mHessFree.BuildTripletsFromMatrix();

				for (size_t i = 0; i < numRB; ++i)
					this->m_vBody[i]->AssembleGlobal_VectorHessianProduct(m_pIntModel->GetGradient(), this->m_mHessFree.m_vvalueTriplets);

				this->m_mHessFree.BuildMatrixFromTriplets();
				this->m_mHessFree.BuildMappingFromMatrix();

				for (size_t i = 0; i < numRB; ++i)
					this->m_vBody[i]->AllocateGlobal_VectorHessianProduct(this->m_mHessFree.m_coeffMap);
			}
			else
			{
				MatrixSd mIntH = m_pIntModel->GetHessian().selfadjointView<Lower>();
				MatrixSd mExtH = this->GetJacobian().transpose()*mIntH*GetJacobian();

				VectorTd vExtH;
				eigenSparseMatrixToTriplets(mExtH, vExtH);

				this->m_mHessFree.m_msparseMatrix *= 0.0;

				for (size_t i = 0; i < vExtH.size(); ++i)
					*this->m_mHessFree.m_coeffMap.at(pair<int, int>(vExtH[i].row(), vExtH[i].col())) = vExtH[i].value();

				for (size_t i = 0; i < numRB; ++i)
					this->m_vBody[i]->AssembleGlobal_FastPreallocatedVectorHessianProduct(m_pIntModel->GetGradient());
			}

			this->m_dirtyFlags = (DirtyFlags)(this->m_dirtyFlags & ~DirtyFlags::HessianFree);
		}
		else
		{
			assert(false);
		}
	}

}