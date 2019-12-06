//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Models/Model.h>

#include <PhySim/Geometry/Node.h>
#include <PhySim/Geometry/Node.h>
#include <PhySim/Geometry/Edge.h>
#include <PhySim/Geometry/Face_Tri.h>
#include <PhySim/Geometry/Face_Quad.h>
#include <PhySim/Geometry/Cell_Tetra.h>
#include <PhySim/Geometry/Cell_Hexa.h>

#include <PhySim/Energies/Material.h>
#include <PhySim/Energies/EnergyElement.h>
#include <PhySim/Energies/MassElement_Lumped.h>
#include <PhySim/Energies/EnergyElement_Force.h>
#include <PhySim/Energies/EnergyElement_Gravity.h>
#include <PhySim/Energies/EnergyElement_External.h>

#include <PhySim/Models/BCondition.h>

#include <PhySim/Solvers/LinearSolver_EigenLDLT.h>
#include <PhySim/Solvers/LinearSolver_EigenCG.h>
#include <PhySim/Solvers/LinearSolver_BiCGSTAB.h>
#include <PhySim/Solvers/LinearSolver_CholmodLDLT.h>
#include <PhySim/Solvers/LinearSolver_SSparseSPQR.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	Model::Model()
	{
		this->m_presolveFixedBC = false;

		this->m_fastAssembly = true;
		this->m_isProfiling = true;
		this->m_timerComputeEnergy = CustomTimer(10, "CAL_ENE");
		this->m_timerComputeGradient = CustomTimer(10, "CAL_GRA");
		this->m_timerComputeHessian = CustomTimer(10, "CAL_HES");
		this->m_timerAssembleEnergy = CustomTimer(10, "ASS_ENE");
		this->m_timerAssembleGradient = CustomTimer(10, "ASS_GRA");
		this->m_timerAssembleHessian = CustomTimer(10, "ASS_HES");

		Free();

		this->m_pPreBCSolver = new LinearSolver_EigenLDLT();

		this->m_gradientFD = false;
		this->m_hessianFD = false;
	}

	Model::~Model()
	{
		Free();

		delete this->m_pPreBCSolver;
	}

	void Model::Free()
	{
		this->m_energy = -1;
		this->m_vgradFree.setZero(0);
		this->m_mHessFree.Clear();
		this->m_mMassFree.Clear();
		this->m_vDoFs.clear();
		this->m_numFreeDoF = -1;
		this->m_numFullDoF = -1;

		this->m_dirtyFlags = DirtyFlags::None;

		// Free energy elements

		size_t numEnerEle = m_venergyEle.size();
		for (int i = 0; i < numEnerEle; ++i)
			delete m_venergyEle[i];
		m_venergyEle.clear();

		// Free mass elements

		size_t numMassEle = m_vmassEle.size();
		for (int i = 0; i < numMassEle; ++i)
			delete m_vmassEle[i];
		m_vmassEle.clear();
		m_vmat_mass.clear();

		// Clear external elements

		this->m_vexternEle.clear();
	}

	void Model::Init()
	{
		Free();

		// Nothing to do here...

		this->DirtyUndeformed();
	}

	inline const vector<DoFSet*>& Model::GetDoFSets() const
	{
		return this->m_vDoFs;
	}

	IModel::StateP Model::CreateState(Space s) const
	{
		StateP pS = StateP(new State());
		this->GetState(pS, s);
		return pS;
	}

	bool Model::HasState(IModel::StateP pS, Space s) const
	{
		State* pSInt = static_cast<State*>(pS.get());

		size_t numDoFSet = this->m_vDoFs.size();

		for (size_t i = 0; i < numDoFSet; ++i)
			if (!pSInt->m_vx.block(this->m_vDoFs[i]->GetOffset_Full(), 0, this->m_vDoFs[i]->GetNumDim(), 1).isApprox(this->m_vDoFs[i]->GetPosition(s), 1e-9))
				return false;

		return true;
	}

	void Model::GetState(IModel::StateP pS, Space s) const
	{
		State* pSInt = static_cast<State*>(pS.get());

		pSInt->m_vx.setZero(this->m_numFullDoF);
		size_t numDoFSet = this->m_vDoFs.size();
		for (size_t i = 0; i < numDoFSet; ++i)
		{
			pSInt->m_vx.block(this->m_vDoFs[i]->GetOffset_Full(), 0, this->m_vDoFs[i]->GetNumDim(), 1) = this->m_vDoFs[i]->GetPosition(s);
		}
	}

	void Model::SetState(const IModel::StateP pS, Space s)
	{
		const State* pSInt = static_cast<const State*>(pS.get());

		assert(pSInt->m_vx.size() == m_numFullDoF);

		size_t numDoFSet = this->m_vDoFs.size();
		for (size_t i = 0; i < numDoFSet; ++i)
		{
			this->m_vDoFs[i]->SetPosition(pSInt->m_vx.block(this->m_vDoFs[i]->GetOffset_Full(), 0, this->m_vDoFs[i]->GetNumDim(), 1), s);
		}

		if (s == Space::DEF)
			this->DirtyDeformed();

		if (s == Space::MAT)
			this->DirtyUndeformed();
	}

	inline int Model::GetNumFullDOF() const
	{
		return this->m_numFullDoF;
	}

	inline int Model::GetNumFreeDOF() const
	{
		return this->m_numFreeDoF;
	}

	void Model::GetFullDOFPosition(VectorXd& vx) const
	{
		vx.resize(this->m_numFullDoF);

		size_t numDoFSet = this->m_vDoFs.size();
		for (size_t i = 0; i < numDoFSet; ++i)
		{
			vx.block(this->m_vDoFs[i]->GetOffset_Full(), 0, this->m_vDoFs[i]->GetNumDim(), 1) = this->m_vDoFs[i]->GetPosition_x();
		}
	}

	void Model::GetFreeDOFPosition(VectorXd& vx) const
	{
		vx.resize(this->m_numFreeDoF);

		size_t numDoFSet = this->m_vDoFs.size();
		for (size_t i = 0; i < numDoFSet; ++i)
		{
			if (this->m_vDoFs[i]->IsFixed())
				continue; // Ignore fixed

			vx.block(this->m_vDoFs[i]->GetOffset_Free(), 0, this->m_vDoFs[i]->GetNumDim(), 1) = this->m_vDoFs[i]->GetPosition_x();
		}
	}

	void Model::GetFreeDOFVelocity(VectorXd& vv) const
	{
		vv.resize(this->m_numFreeDoF);
		vv.setZero();

		size_t numDoFSet = this->m_vDoFs.size();
		for (size_t i = 0; i < numDoFSet; ++i)
		{
			if (this->m_vDoFs[i]->IsFixed())
				continue; // Ignore fixed

			vv.block(this->m_vDoFs[i]->GetOffset_Free(), 0, this->m_vDoFs[i]->GetNumDim(), 1) = this->m_vDoFs[i]->GetVelocity_x();
		}
	}

	void Model::SetFullDOFPosition(const VectorXd& vx)
	{
		assert(vx.size() == this->m_numFullDoF);

		size_t numDoFSet = this->m_vDoFs.size();
		for (size_t i = 0; i < numDoFSet; ++i)
		{
			this->m_vDoFs[i]->SetPosition_x(vx.block(this->m_vDoFs[i]->GetOffset_Full(), 0, this->m_vDoFs[i]->GetNumDim(), 1));
		}

		this->DirtyDeformed();
	}

	void Model::SetFreeDOFPosition(const VectorXd& vx)
	{
		assert(vx.size() == this->m_numFreeDoF);

		size_t numDoFSet = this->m_vDoFs.size();
		for (size_t i = 0; i < numDoFSet; ++i)
		{
			if (this->m_vDoFs[i]->IsFixed())
				continue; // Ignore fixed

			this->m_vDoFs[i]->SetPosition_x(vx.block(this->m_vDoFs[i]->GetOffset_Free(), 0, this->m_vDoFs[i]->GetNumDim(), 1));
		}

		this->DirtyDeformed();
	}

	void Model::SetFreeDOFVelocity(const VectorXd& vv)
	{
		assert(vv.size() == this->m_numFreeDoF);

		size_t numDoFSet = this->m_vDoFs.size();
		for (size_t i = 0; i < numDoFSet; ++i)
		{
			if (this->m_vDoFs[i]->IsFixed())
				continue; // Ignore fixed

			this->m_vDoFs[i]->SetVelocity_x(vv.block(this->m_vDoFs[i]->GetOffset_Free(), 0, this->m_vDoFs[i]->GetNumDim(), 1));
		}

		this->DirtyDeformed();
	}

	inline const Real& Model::GetEnergy()
	{
		if (this->IsDirty_Energy())
			this->ComputeAndStore_Energy();

		return this->m_energy;
	}

	inline const VectorXd& Model::GetGradient(bool full)
	{
		if (full)
		{
			if (this->IsDirty_Gradient(full))
			{
				this->ComputeAndStore_Hessian(full);
			}
		}
		else
		{
			if (this->IsDirty_Gradient(full))
			{
				if (!this->m_gradientFD)
					this->ComputeAndStore_Gradient(full);
				else this->ComputeAndStore_Gradient_FD(full);
			}
		}

		if (!full)
			return this->m_vgradFree;
		else return this->m_vgradFull;
	}

	inline const MatrixSd& Model::GetHessian(bool full)
	{
		if (full)
		{
			if (this->IsDirty_Hessian(full))
			{
				this->ComputeAndStore_Hessian(full);
			}
		}
		else
		{
			if (this->IsDirty_Hessian(full))
			{
				if (!this->m_hessianFD)
					this->ComputeAndStore_Hessian(full);
				else this->ComputeAndStore_Hessian_FD(full);
			}
		}

		if (!full)
			return this->m_mHessFree.m_msparseMatrix;
		else return this->m_mHessFull.m_msparseMatrix;
	}

	inline const MatrixSd& Model::GetMass(bool full)
	{
		if (this->IsDirty_Mass(full))
			this->ComputeAndStore_Mass(full);

		if (!full)
			return this->m_mMassFree.m_msparseMatrix;
		else return this->m_mMassFull.m_msparseMatrix;
	}

	void Model::PrepareForSimulation()
	{
		// Init boundary

		for (size_t i = 0; i < this->m_vDoFs.size(); ++i)
			this->m_vDoFs[i]->Unfix(); // Initialize unfixed

		this->m_energy = 0;
		this->m_vgradFull.resize(m_numFullDoF);
		this->m_mHessFull = FastMatrixSd(m_numFullDoF, m_numFullDoF);
		this->m_mMassFull = FastMatrixSd(m_numFullDoF, m_numFullDoF);

		//this->ComputeAndStore_Mass(true);
		//this->ComputeAndStore_Hessian(true);
		//this->ComputeAndStore_Gradient();
		//this->ComputeAndStore_Energy();

		if (this->m_presolveFixedBC)
		{
			size_t numDoFSet = this->m_vDoFs.size();

			// Gather fixed DoF

			bVector vfixedDoFSet(numDoFSet, false);
			m_numFixedDoF = 0;
			m_numUnfixDoF = 0;
			int offsetUnfix = 0;
			int offsetFixed = 0;
			for (size_t i = 0; i < m_vBC.size(); ++i)
			{
				if (m_vBC[i]->Setup().m_type != BCType::Fixed)
					continue; // Ignore any other B. Condition

				for (size_t j = 0; j < m_vBC[i]->Setup().m_vDoF.size(); ++j)
				{
					vfixedDoFSet[m_vBC[i]->Setup().m_vDoF[j]->GetId_Full()] = true;
					m_numFixedDoF += this->m_vBC[i]->Setup().m_vDoF[j]->GetNumDim();
				}
			}
			offsetFixed = m_numUnfixDoF = this->m_numFullDoF - m_numFixedDoF;

			// Build permutation

			VectorTd vP; vP.reserve(m_numFullDoF);
			for (size_t i = 0; i < numDoFSet; ++i)
			{
				DoFSet* pDoF = this->m_vDoFs[i];

				if (vfixedDoFSet[i])
				{
					for (size_t j = 0; j < pDoF->GetNumDim(); ++j)
						vP.push_back(Triplet<Real>(offsetFixed + (int)j, pDoF->GetOffset_Full() + (int)j, 1.0));
					offsetFixed += pDoF->GetNumDim();
				}
				else
				{
					for (size_t j = 0; j < pDoF->GetNumDim(); ++j)
						vP.push_back(Triplet<Real>(offsetUnfix + (int)j, pDoF->GetOffset_Full() + (int)j, 1.0));
					offsetUnfix += pDoF->GetNumDim();
				}
			}

			m_mFixedPerm = MatrixSd(m_numFullDoF, m_numFullDoF);
			m_mFixedPerm.setFromTriplets(vP.begin(), vP.end());
			m_mFixedPerm.makeCompressed();

			MatrixSd mHFull = this->GetHessian(true).selfadjointView<Lower>();
			MatrixSd mHPerm = m_mFixedPerm*mHFull*m_mFixedPerm.transpose();
			MatrixSd mA = mHPerm.block(0, 0, m_numUnfixDoF, m_numUnfixDoF);
			m_pPreBCSolver->Init(mA, LinearSolverOptions());
		}

		for (size_t i = 0; i < this->m_vBC.size(); ++i)
			this->m_vBC[i]->Init();

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

		this->ComputeAndStore_Mass(false);
		this->ComputeAndStore_Hessian(false);
		this->ComputeAndStore_Gradient();
		this->ComputeAndStore_Energy();
	}

	void Model::DirtyUndeformed()
	{
		this->m_dirtyFlags = this->m_dirtyFlags | DirtyFlags::MassFull;
		this->m_dirtyFlags = this->m_dirtyFlags | DirtyFlags::MassFree;

		this->DirtyDeformed();
	}

	void Model::DirtyDeformed()
	{
		this->m_dirtyFlags = this->m_dirtyFlags | DirtyFlags::Energy;
		this->m_dirtyFlags = this->m_dirtyFlags | DirtyFlags::GradientFull;
		this->m_dirtyFlags = this->m_dirtyFlags | DirtyFlags::GradientFree;
		this->m_dirtyFlags = this->m_dirtyFlags | DirtyFlags::HessianFull;
		this->m_dirtyFlags = this->m_dirtyFlags | DirtyFlags::HessianFree;
	}

	bool Model::IsDirty_Energy() const
	{
		return (this->m_dirtyFlags & (int)DirtyFlags::Energy) != 0;
	}

	bool Model::IsDirty_Gradient(bool full) const
	{
		if (!full)
			return (this->m_dirtyFlags & (int)DirtyFlags::GradientFree) != 0;
		return (this->m_dirtyFlags & (int)DirtyFlags::GradientFull) != 0;
	}

	bool Model::IsDirty_Hessian(bool full) const
	{
		if (!full)
			return (this->m_dirtyFlags & (int)DirtyFlags::HessianFree) != 0;
		else return (this->m_dirtyFlags & (int)DirtyFlags::HessianFull) != 0;
	}

	bool Model::IsDirty_Mass(bool full) const
	{
		if (!full)
			return (this->m_dirtyFlags & (int)DirtyFlags::MassFree) != 0;
		else return (this->m_dirtyFlags & (int)DirtyFlags::MassFull) != 0;
	}

	void Model::ComputeAndStore_Gradient_FD(bool full)
	{
		logSimu("\n[TRACE] Estimating %s model gradient using FD", this->GetName().c_str());

		Real eps = 1e-6;

		VectorXd vxF;
		this->GetFreeDOFPosition(vxF);

		IModel::StateP pS0 = this->CreateState();
		DirtyFlags prevFlags = this->m_dirtyFlags;
		this->GetState(pS0);

		// Compute estimation
		VectorXd vgF(m_numFreeDoF);

		for (int i = 0; i < m_numFreeDoF; ++i)
		{
			// +
			SetState(pS0);
			vxF(i) += eps;
			this->SetFreeDOFPosition(vxF);
			this->ComputeAndStore_Energy();
			Real ep = this->m_energy;

			// -
			SetState(pS0);
			vxF(i) -= 2 * eps;
			this->SetFreeDOFPosition(vxF);
			this->ComputeAndStore_Energy();
			Real em = this->m_energy;

			vgF(i) = (ep - em) / (2 * eps);

			vxF(i) += eps;
		}

		// Recover previous
		this->SetState(pS0);
		this->m_dirtyFlags = prevFlags;

		// Set free gradient
		this->m_vgradFree = vgF;

		// Clean
		this->m_dirtyFlags = (DirtyFlags)(this->m_dirtyFlags & ~DirtyFlags::GradientFree);

		// Store for debug
		//PhySim::logFile("finiteDiffGradient.csv", PhySim::matrixToString_CSV(m_vgradFree));
	}

	void Model::ComputeAndStore_Hessian_FD(bool full)
	{
		logSimu("\n[TRACE] Estimating %s model Hessian using FD", this->GetName().c_str());

		if (!this->m_mHessFree.HasMappingData())
			this->ComputeAndStore_Hessian(false);

		Real eps = 1e-6;

		VectorXd vxF;
		this->GetFreeDOFPosition(vxF);

		IModel::StateP pS0 = this->CreateState();
		DirtyFlags prevFlags = this->m_dirtyFlags;
		this->GetState(pS0);

		// Compute estimation
		MatrixXd mHF(m_numFreeDoF, m_numFreeDoF);

		for (int i = 0; i < m_numFreeDoF; ++i)
		{

			// +
			SetState(pS0);
			vxF(i) += eps;
			this->SetFreeDOFPosition(vxF);
			this->ComputeAndStore_Gradient();
			VectorXd vgp = this->m_vgradFree;

			// -
			SetState(pS0);
			vxF(i) -= 2 * eps;
			this->SetFreeDOFPosition(vxF);
			this->ComputeAndStore_Gradient();
			VectorXd vgm = this->m_vgradFree;

			mHF.col(i) = (vgp - vgm) / (2 * eps);

			vxF(i) += eps;
		}

		// Recover previous
		this->SetState(pS0);
		this->m_dirtyFlags = prevFlags;

		// Set free Hessian
		
		this->m_mHessFree.m_msparseMatrix *= 0.0;

		for (int i = 0; i < this->m_mHessFree.m_vpointTriplets.size(); ++i)
		{
			*(this->m_mHessFree.m_vpointTriplets[i].value()) = mHF(
				this->m_mHessFree.m_vpointTriplets[i].row(),
				this->m_mHessFree.m_vpointTriplets[i].col());
		}

		// Clean
		this->m_dirtyFlags = (DirtyFlags)(this->m_dirtyFlags & ~DirtyFlags::HessianFree);

		// Store for debug
		//MatrixSd mHfull = this->m_mHessFree.m_msparseMatrix.selfadjointView<Lower>();
		//PhySim::logFile("finiteDiffHessian.csv", PhySim::matrixToString_CSV(mHfull));
	}

	void Model::TestGlobalGradient()
	{
		Real eps = 1e-6;

		this->ComputeAndStore_Gradient_FD();
		VectorXd vgF = this->m_vgradFree;
		this->ComputeAndStore_Gradient();
		VectorXd vgA = this->m_vgradFree;

		// Test
		VectorXd vgD = vgA - vgF;
		Real testNorm = vgA.norm();
		Real realNorm = vgF.norm();
		Real diffNorm = vgD.norm();
		if (realNorm < 1e-6)
		{
			logSimu("\n[INVALID] Finite difference gradient near zero: %f", realNorm);
		}
		else
		{
			if (diffNorm / realNorm > 1e-6)
			{
				logSimu("\n[FAILURE] Global gradient test error: %f", diffNorm / realNorm);
				logFile("csvGradientTest_F.csv", vectorToString_CSV(vgF));
				logFile("csvGradientTest_A.csv", vectorToString_CSV(vgA));
			}
			else
			{
				logSimu("\n[SUCCESS] Global gradient test error: %f", diffNorm / realNorm);
			}
		}
	}

	void Model::TestGlobalHessian()
	{
		this->ComputeAndStore_Hessian_FD();
		MatrixSd mHFS = this->m_mHessFree.m_msparseMatrix.selfadjointView<Lower>();
		MatrixXd mHF = mHFS.toDense();
		this->ComputeAndStore_Hessian();
		MatrixSd mHAS = this->m_mHessFree.m_msparseMatrix.selfadjointView<Lower>();
		MatrixXd mHA = mHAS.toDense();

		// Test
		MatrixXd mHD = mHA - mHF;
		Real testNorm = mHA.norm();
		Real realNorm = mHF.norm();
		Real diffNorm = mHD.norm();
		if (realNorm < 1e-6)
		{
			logSimu("\n[INVALID] Finite difference Hessian near zero: %f", realNorm);
		}
		else
		{
			if (diffNorm / realNorm > 1e-6)
			{
				logSimu("\n[FAILURE] Global Hessian test error: %f", diffNorm / realNorm);
				logFile("csvHessianTest_F.csv", matrixToString_CSV(mHF));
				logFile("csvHessianTest_A.csv", matrixToString_CSV(mHA));
			}
			else
			{
				logSimu("\n[SUCCESS] Global Hessian test error: %f", diffNorm / realNorm);
			}
		}
	}

	void Model::AddExternalEnergyComponent(EnergyElement_External* pEle)
	{
		vector<EnergyElement_External*>::iterator it = find(
			this->m_vexternEle.begin(),
			this->m_vexternEle.end(),
			pEle);
		if (it == this->m_vexternEle.end())
		{
			this->m_vexternEle.reserve(this->m_vexternEle.size() + 1);
			this->m_venergyEle.reserve(this->m_venergyEle.size() + 1);
			this->m_vexternEle.push_back(pEle);
			this->m_venergyEle.push_back(pEle);
			this->m_vexternEle.back()->Init();

			this->DirtyDeformed();
		}
	}

	void Model::DelExternalEnergyComponent(EnergyElement_External* pEle)
	{
		vector<EnergyElement_External*>::iterator it = find(
			this->m_vexternEle.begin(),
			this->m_vexternEle.end(),
			pEle);
		if (it != this->m_vexternEle.end())
		{
			this->m_vexternEle.erase(it);
			this->m_venergyEle.erase(find(
				this->m_venergyEle.begin(),
				this->m_venergyEle.end(),
				pEle));

			this->DirtyDeformed();
		}
	}


	void Model::AddBoundaryCondition(BCondition* pBC)
	{
		vector<BCondition*>::iterator it = find(
			this->m_vBC.begin(),
			this->m_vBC.end(),
			pBC);

		assert(it == m_vBC.end());
		this->m_vBC.push_back(pBC);

		if (pBC->Element() != NULL)
		{
			this->m_venergyEle.push_back(pBC->Element());
		}
	}

	void Model::RemoveBoundaryCondition(BCondition* pBC)
	{
		vector<BCondition*>::iterator it = find(
			this->m_vBC.begin(),
			this->m_vBC.end(),
			pBC);

		assert(it != m_vBC.end());
		this->m_vBC.erase(it);

		if (pBC->Element() != NULL)
		{
			this->m_venergyEle.erase(find(
				this->m_venergyEle.begin(),
				this->m_venergyEle.end(),
				pBC->Element()));
		}
	}

	void Model::ClearBoundaryConditions()
	{
		for (size_t i = 0; i < this->m_vBC.size(); ++i)
		{
			if (this->m_vBC[i]->Element() != NULL)
			{
				this->m_venergyEle.erase(find(
					this->m_venergyEle.begin(),
					this->m_venergyEle.end(),
					m_vBC[i]->Element()));
			}
		}
		this->m_vBC.clear();
	}

	bool Model::BoundaryConditionsLoaded()
	{
		for (size_t i = 0; i < this->m_vBC.size(); ++i)
			if (!this->m_vBC[i]->IsFullyLoaded())
				return false;

		return true;
	}

	void Model::ResetBoundaryConditions()
	{
		bVector vapply(this->m_vBC.size(), false);

		for (size_t i = 0; i < this->m_vBC.size(); ++i)
			vapply[i] = this->m_vBC[i]->ResetLoading();

		for (size_t i = 0; i < this->m_vBC.size(); ++i)
		{
			if (vapply[i]) // Apply BC
				this->m_vBC[i]->Apply();
		}
	}

	void Model::FullBoundaryConditions()
	{
		bVector vapply(this->m_vBC.size(), false);

		bool applyFixedPoints = false;

		for (size_t i = 0; i < this->m_vBC.size(); ++i)
		{
			vapply[i] = this->m_vBC[i]->FullLoading();

			if (vapply[i] && m_vBC[i]->Setup().m_type == BCType::Fixed)
				applyFixedPoints = true; // At least one fixed BC applied
		}

		if (m_presolveFixedBC && applyFixedPoints)
			this->PresolveForFixedBoundary();

		for (size_t i = 0; i < this->m_vBC.size(); ++i)
		{
			if (vapply[i]) // Apply BC
				this->m_vBC[i]->Apply();
		}
	}

	void Model::StepBoundaryConditions()
	{
		bVector vapply(this->m_vBC.size(), false);

		bool applyFixedPoints = false;

		for (size_t i = 0; i < this->m_vBC.size(); ++i)
		{
			vapply[i] = this->m_vBC[i]->StepLoading();

			if (vapply[i] && m_vBC[i]->Setup().m_type == BCType::Fixed)
				applyFixedPoints = true; // At least one fixed BC applied
		}

		if (m_presolveFixedBC && applyFixedPoints)
			this->PresolveForFixedBoundary();

		for (size_t i = 0; i < this->m_vBC.size(); ++i)
		{
			if (vapply[i]) // Apply BC
				this->m_vBC[i]->Apply();
		}
	}

	void Model::PresolveForFixedBoundary()
	{
		size_t numDoFSet = this->m_vDoFs.size();

		VectorXd vxFull = VectorXd(m_numFullDoF);
		VectorXd dxFull = VectorXd(m_numFullDoF);
		vxFull.setZero();
		dxFull.setZero();
		this->GetFullDOFPosition(vxFull);

		// Gather fixed DoF displacement

		bVector vfixedDoFSet(numDoFSet, false);
		for (size_t i = 0; i < m_vBC.size(); ++i)
		{
			if (m_vBC[i]->Setup().m_type != BCType::Fixed)
				continue; // Ignore any other B. Condition

			vector<VectorXd> vEndVal; m_vBC[i]->GetCurValues(vEndVal);
			for (size_t j = 0; j < m_vBC[i]->Setup().m_vDoF.size(); ++j)
			{
				int numDim = this->m_vBC[i]->Setup().m_vDoF[j]->GetNumDim();
				int offset = this->m_vBC[i]->Setup().m_vDoF[j]->GetOffset_Full();
				VectorXd vxDoF = vxFull.block(offset, 0, numDim, 1);
				dxFull.block(offset, 0, numDim, 1) = vEndVal[j] - vxDoF;
			}
		}

		VectorXd vgFull = this->GetGradient(true);
		MatrixSd mHFull = this->GetHessian(true).selfadjointView<Lower>();
		MatrixSd mHPerm = m_mFixedPerm*mHFull*m_mFixedPerm.transpose();
		VectorXd vgPerm = m_mFixedPerm*vgFull;
		VectorXd dxPerm = m_mFixedPerm*dxFull;

		// Partition permuted mH, vg and dx

		int offsetFixed = m_numUnfixDoF;
		int offsetUnfix = 0;

		VectorXd vgU = vgPerm.block(offsetUnfix, 0, m_numUnfixDoF, 1); // Gradient at unfixed
		VectorXd dxF = dxPerm.block(offsetFixed, 0, m_numFixedDoF, 1); // Displacement at fixed
		MatrixSd mA = mHPerm.block(offsetUnfix, offsetUnfix, m_numUnfixDoF, m_numUnfixDoF);
		MatrixSd mC = mHPerm.block(offsetUnfix, offsetFixed, m_numUnfixDoF, m_numFixedDoF);
		VectorXd vb = -vgU - mC*dxF;

		VectorXd vx(vb.size());
		vx.setZero();
		this->m_pPreBCSolver->Solve(mA, vb, vx);
		dxPerm.block(0, 0, m_numUnfixDoF, 1) = vx;
		vxFull += m_mFixedPerm.transpose()*dxPerm;

		this->SetFullDOFPosition(vxFull);
	}

}