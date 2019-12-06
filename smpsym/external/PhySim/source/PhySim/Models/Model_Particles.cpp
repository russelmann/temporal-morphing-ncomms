//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Models/Model_Particles.h>

#include <PhySim/Geometry/Node.h>
#include <PhySim/Geometry/Node.h>
#include <PhySim/Geometry/Edge.h>
#include <PhySim/Geometry/Face_Tri.h>
#include <PhySim/Geometry/Face_Quad.h>
#include <PhySim/Geometry/Cell_Tetra.h>
#include <PhySim/Geometry/Cell_Hexa.h>
#include <PhySim/Geometry/NodeEBD.h>

#include <PhySim/Energies/EnergyElement.h>
#include <PhySim/Energies/MassElement_Lumped.h>
#include <PhySim/Energies/EnergyElement_Gravity.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	Model_Particles::Model_Particles() : Model()
	{
		this->m_pOptions = NULL;
	}

	Model_Particles::~Model_Particles()
	{
		delete this->m_pOptions;
	}

	Model_Particles::Options& Model_Particles::GetOptions()
	{
		if (this->m_pOptions == NULL) // Create if needed
			this->m_pOptions = new Model_Particles::Options();
		return *((Model_Particles::Options*) this->m_pOptions);
	}

	void Model_Particles::Free()
	{
		Model::Free();

		// Free particles vector

		size_t numNode = m_vnodes.size();
		for (int i = 0; i < numNode; ++i)
			delete m_vnodes[i];
		m_vnodes.clear();

		// Free elements vector

		size_t numEles = m_velems.size();
		for (int i = 0; i < numEles; ++i)
			delete m_velems[i];
		m_velems.clear();
	}

	void Model_Particles::Init()
	{
		Free();

		this->InitDiscretization();
		this->InitDegreesOfFreedom();
		this->InitEnergyElements();
		this->InitMassElements();

		this->InitializeRest();
		this->DirtyUndeformed();

		logSimu("\n--");
		logSimu("\nINITIALIZED %s", this->GetName().c_str());
		logSimu("\nTotal DOF: %i", this->m_numFullDoF);
		logSimu("\nNode number: %i", (int) this->m_vnodes.size());
		logSimu("\nDiscretization elements: %i", (int) this->m_velems.size());
		logSimu("\nEnergy elements number: %i", (int) this->m_venergyEle.size());
		logSimu("\nMass elements number: %i", (int) this->m_vmassEle.size());

		logSimu("\n--");
	}

	void Model_Particles::InitDiscretization()
	{
		const MatrixXd& mN = (this->m_pOptions->m_hasSubElem) ? m_pOptions->m_mNodesSub : m_pOptions->m_mNodes;
		const MatrixXi& mE = (this->m_pOptions->m_hasSubElem) ? m_pOptions->m_mElemsSub : m_pOptions->m_mElems;

		size_t numNodes = mN.rows();
		this->m_vnodes.resize(numNodes);

		// Create nodes

		for (size_t i = 0; i < numNodes; ++i)
		{
			this->m_vnodes[i] = new Node(mN.row(i));
		}

		// Create elements

		if (m_pOptions->m_discretization == Discretization::Nodes)
		{
			this->m_velems.resize(numNodes);
			for (size_t i = 0; i < numNodes; ++i)
				this->m_velems[i] = m_vnodes[i];
		}
		else
		{
			size_t numElems = mE.rows();
			this->m_velems.resize(numElems);

			for (size_t i = 0; i < numElems; ++i)
			{
				int numEleNodes = (int) mE.cols();

				vector<Node*> veleNodes(numEleNodes);
				for (int j = 0; j < numEleNodes; ++j)
					veleNodes[j] = this->m_vnodes[mE(i, j)];

				switch (m_pOptions->m_discretization)
				{
				case Discretization::Edges: this->m_velems[i] = new Edge(veleNodes); break;
				case Discretization::Triangles: this->m_velems[i] = new Face_Tri(veleNodes); break;
				case Discretization::Quadrangles: this->m_velems[i] = new Face_Quad(veleNodes); break;
				case Discretization::Tetrahedra4: this->m_velems[i] = new Cell_Tetra(veleNodes); break;
				case Discretization::Tetrahedra10: this->m_velems[i] = new Cell_Tetra(veleNodes); break;
				case Discretization::Hexahedra: this->m_velems[i] = new Cell_Hexa(veleNodes); break;
				default:
               assert(false); ;
				}
			}
		}
	}

	void Model_Particles::InitDegreesOfFreedom()
	{
		size_t numNodes = m_vnodes.size();
		this->m_vDoFs.resize(numNodes);
		this->m_numFullDoF = 0;

		for (size_t i = 0; i < numNodes; ++i)
		{
			this->m_vDoFs[i] = m_vnodes[i]->DoF();
			this->m_vDoFs[i]->SetId_Full((int)i);
			this->m_vDoFs[i]->SetOffset_Full((int)m_numFullDoF);
			this->m_numFullDoF += 3;
		}
	}

	void Model_Particles::InitMassElements()
	{
		assert(m_pOptions->m_material.HasProperty(Material::Property::Density));

		size_t numElems = m_velems.size();
		this->m_vmat_mass.resize(numElems);
		this->m_vmassEle.resize(numElems);

		Real density = m_pOptions->m_material[Material::Property::Density];

		for (size_t i = 0; i < numElems; ++i)
		{
			this->m_vmat_mass[i].AddProperty(Material::Property::Density, density);
			this->m_vmassEle[i] = new MassElement_Lumped(this, m_velems[i], &m_vmat_mass[i]);
		}
	}

	void Model_Particles::InitEnergyElements()
	{
		//EnergyElement_Gravity* pEle = new EnergyElement_Gravity(this);
		//pEle->SetGravityAcceleration(m_pOptions->m_vgravity);
		//this->m_venergyEle.resize(1);
		//this->m_venergyEle[0] = pEle;
		//this->m_pEle_gravity = pEle;
	}

	void Model_Particles::SetSubelementPositions(const MatrixXd& mV, Space s)
	{
		assert(mV.rows() == m_pOptions->m_mNodes.rows());

		for (int i = 0; i < mV.rows(); ++i)
		{
			m_vnodes[i]->DoF()->SetPosition(mV.row(i), s);
		}

		for (int i = 0; i < (int) m_velems.size(); ++i)
		{
			this->m_velems[i]->SetSubelementPositions(s);
		}

		if (s == Space::DEF)
			this->DirtyDeformed();

		if (s == Space::MAT)
			this->DirtyUndeformed();
	}

	inline const vector<Node*>& Model_Particles::GetNodes() const
	{
		return this->m_vnodes;
	}

	inline const vector<Polytope*>& Model_Particles::GetElements() const
	{
		return this->m_velems;
	}

	void Model_Particles::InitializeRest()
	{
		int numEle = (int) this->m_venergyEle.size();

#pragma omp parallel for
		for (int i = 0; i < numEle; ++i)
			this->m_venergyEle[i]->Init();
	}

	void Model_Particles::ComputeAndStore_Energy()
	{
		int numEle = (int) this->m_venergyEle.size();

		// Compute

		if (m_isProfiling) this->m_timerComputeEnergy.start();
#pragma omp parallel for
		for (int i = 0; i < numEle; ++i)
			this->m_venergyEle[i]->ComputeAndStore_Energy();
		if (m_isProfiling) this->m_timerComputeEnergy.stopStoreLog();

		// Sum

		if (m_isProfiling) this->m_timerAssembleEnergy.start();
		this->m_energy = 0;
		for (size_t i = 0; i < numEle; ++i)
			this->m_energy += this->m_venergyEle[i]->GetElementEnergy();
		if (m_isProfiling) this->m_timerAssembleEnergy.stopStoreLog();

		// Clean

		this->m_dirtyFlags = (DirtyFlags)(this->m_dirtyFlags & ~DirtyFlags::Energy);
	}

	void Model_Particles::ComputeAndStore_Gradient(bool full)
	{
		int numEle = (int) this->m_venergyEle.size();

		// Compute

		if (m_isProfiling) this->m_timerComputeGradient.start();
#pragma omp parallel for
		for (int i = 0; i < numEle; ++i)
			this->m_venergyEle[i]->ComputeAndStore_Gradient();
		if (m_isProfiling) this->m_timerComputeGradient.stopStoreLog();


		if (!full)
		{
			// Assemble

			if (m_isProfiling) this->m_timerAssembleGradient.start();
			this->m_vgradFree.setZero();
			for (size_t i = 0; i < numEle; ++i)
				this->m_venergyEle[i]->AssembleGlobal_Gradient(this->m_vgradFree, full);
			if (m_isProfiling) this->m_timerAssembleGradient.stopStoreLog();

			// Clean

			this->m_dirtyFlags = (DirtyFlags)(this->m_dirtyFlags & ~DirtyFlags::GradientFree);
		}
		else
		{
			// Assemble

			if (m_isProfiling) this->m_timerAssembleGradient.start();
			this->m_vgradFull.setZero();
			for (size_t i = 0; i < numEle; ++i)
				this->m_venergyEle[i]->AssembleGlobal_Gradient(this->m_vgradFull, full);
			if (m_isProfiling) this->m_timerAssembleGradient.stopStoreLog();

			// Clean

			this->m_dirtyFlags = (DirtyFlags)(this->m_dirtyFlags & ~DirtyFlags::GradientFull);
		}
	}

	void Model_Particles::ComputeAndStore_Hessian(bool full)
	{
		int numEle = (int) this->m_venergyEle.size();

		// Compute

		if (m_isProfiling) this->m_timerComputeHessian.start();
//#pragma omp parallel for
		for (int i = 0; i < numEle; ++i)
			this->m_venergyEle[i]->ComputeAndStore_Hessian();
		if (m_isProfiling) this->m_timerComputeHessian.stopStoreLog();

		if (!full)
		{
			// Assemble

			if (m_isProfiling) this->m_timerAssembleHessian.start();
			if (m_fastAssembly && !this->m_mHessFree.HasMappingData())
			{
				// Gather triplets

				this->m_mHessFree.m_vvalueTriplets.clear();

				for (size_t i = 0; i < numEle; ++i)
					this->m_venergyEle[i]->AssembleGlobal_Hessian(this->m_mHessFree.m_vvalueTriplets, full);

				// Allocate mapping

				this->m_mHessFree.BuildMatrixFromTriplets();
				this->m_mHessFree.BuildMappingFromMatrix();

				for (size_t i = 0; i < numEle; ++i)
					this->m_venergyEle[i]->AllocateGlobal_Hessian(this->m_mHessFree.m_coeffMap, full);
			}
			else
			{
				this->m_mHessFree.m_msparseMatrix *= 0.0;
				for (size_t i = 0; i < numEle; ++i)
					this->m_venergyEle[i]->AssembleGlobal_FastPreallocatedHessian(full);
			}
			if (m_isProfiling) this->m_timerAssembleHessian.stopStoreLog();

			// Clean

			this->m_dirtyFlags = (DirtyFlags)(this->m_dirtyFlags & ~DirtyFlags::HessianFree);
		}
		else
		{
			// Assemble

			if (m_isProfiling) this->m_timerAssembleHessian.start();
			if (m_fastAssembly && !this->m_mHessFull.HasMappingData())
			{
				// Gather triplets

				this->m_mHessFull.m_vvalueTriplets.clear();

				for (size_t i = 0; i < numEle; ++i)
					this->m_venergyEle[i]->AssembleGlobal_Hessian(this->m_mHessFull.m_vvalueTriplets, full);

				// Allocate mapping

				this->m_mHessFull.BuildMatrixFromTriplets();
				this->m_mHessFull.BuildMappingFromMatrix();

				for (size_t i = 0; i < numEle; ++i)
					this->m_venergyEle[i]->AllocateGlobal_Hessian(this->m_mHessFull.m_coeffMap, full);
			}
			else
			{
				this->m_mHessFull.m_msparseMatrix *= 0.0;
				for (size_t i = 0; i < numEle; ++i)
					this->m_venergyEle[i]->AssembleGlobal_FastPreallocatedHessian(full);
			}
			if (m_isProfiling) this->m_timerAssembleHessian.stopStoreLog();

			// Clean

			this->m_dirtyFlags = (DirtyFlags)(this->m_dirtyFlags & ~DirtyFlags::HessianFull);
		}
	}

	void Model_Particles::ComputeAndStore_Mass(bool full)
	{
		int numEle = (int) this->m_vmassEle.size();

		// Compute mass and lump

#pragma omp parallel for
		for (int i = 0; i < numEle; ++i)
			this->m_vmassEle[i]->ComputeAndStore_Mass();

		if (!full)
		{
			// Assemble mass matrix 

			for (size_t i = 0; i < numEle; ++i)
				this->m_vmassEle[i]->AssembleGlobal_MassMatrix(this->m_mMassFree.m_vvalueTriplets, full);

			this->m_mMassFree.BuildMatrixFromTriplets();

			this->m_dirtyFlags = (DirtyFlags)(this->m_dirtyFlags & ~DirtyFlags::MassFree);
		}
		else
		{
			// Assemble mass matrix 

			for (size_t i = 0; i < numEle; ++i)
				this->m_vmassEle[i]->AssembleGlobal_MassMatrix(this->m_mMassFull.m_vvalueTriplets, full);

			this->m_mMassFull.BuildMatrixFromTriplets();

			this->m_dirtyFlags = (DirtyFlags)(this->m_dirtyFlags & ~DirtyFlags::MassFull);
		
			double mass = 0;
			for (size_t i = 0; i < this->m_vmassEle.size(); ++i)
				mass += this->m_vmassEle[i]->GetElementMass();

			double volume = 0;
			for (int i = 0; i < this->m_velems.size(); ++i)
				volume += this->m_velems[i]->ComputeVolume();

			logSimu("\nTotal model volume: %f", volume);
			logSimu("\nTotal model mass: %f", mass);
		}
	}

	vector<DoFSet*> Model_Particles::SelectDoF(const Vector3d& vboxMin, const Vector3d& vboxMax, Space s) const
	{
		vector<DoFSet*> vDoFs;

		for (size_t i = 0; i < this->m_vnodes.size(); ++i)
		{
			const Vector3d& pos = this->m_vnodes[i]->DoF()->GetPosition(s);
			if (pos.x() > vboxMin.x() && pos.x() < vboxMax.x() &&
				pos.y() > vboxMin.y() && pos.y() < vboxMax.y() &&
				pos.z() > vboxMin.z() && pos.z() < vboxMax.z())
				vDoFs.push_back(this->m_vnodes[i]->DoF());
		}

		return vDoFs;
	}

	vector<NodeEBD> Model_Particles::ComputeEmbedding(const vector<Vector3d>& vp, Space s)
	{
		vector<NodeEBD> vnode;
		vnode.resize(vp.size());

		for (size_t i = 0; i < vp.size(); ++i)
		{
			bool embedded = false;

			for (size_t j = 0; j < this->m_velems.size(); ++j)
			{
				NodeEBD node = this->m_velems[j]->ComputeEmbedding(vp[i], s);

				if (node.Valid())
				{
					vnode[i] = node;
					embedded = true;
					break;
				}
			}
		}

		return vnode;
	}

}
