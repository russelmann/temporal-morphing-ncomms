//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Energies/MassElement_Lumped.h>

#include <PhySim/Geometry/Node.h>

#include <PhySim/Models/Model_Particles.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;


	MassElement_Lumped::MassElement_Lumped(Model_Particles* pModel, Polytope* pPolytope, Material* pMaterial)
	{
		this->m_pModel = pModel;
		this->m_pMaterial = pMaterial;
		this->m_pPolytope = pPolytope;
	}

	MassElement_Lumped::~MassElement_Lumped(void)
	{
	}

	void MassElement_Lumped::Init() 
	{
	}

	//void MassElement_Lumped::AddLumpedMassToNodes()
	//{
	//	size_t numNode = this->pPolytope->Nodes().size();
	//	Real massPerNode = this->mass / numNode;

	//	for (size_t i = 0; i < numNode; ++i)
	//	{
	//		this->pPolytope->Nodes()[i]->DoF()->AddMass(massPerNode);
	//	}
	//}

	void MassElement_Lumped::ComputeAndStore_Mass()
	{
		this->m_mass = this->m_pPolytope->ComputeVolume()*(*this->m_pMaterial)[Material::Property::Density];
	}

	void MassElement_Lumped::AssembleGlobal_MassLumped(VectorXd& vMass, bool full)
	{
		size_t numNode = this->m_pPolytope->Nodes().size();
		Real massPerNode = this->m_mass / numNode;

		if (!full)
		{
			for (size_t i = 0; i < numNode; ++i)
			{
				DoFSet* DoF = this->m_pPolytope->Nodes()[i]->DoF();

				// Is this a simulated DoF?
				if (DoF->GetId_Free() < 0)
				{
					continue;
				}

				vMass.block(DoF->GetOffset_Free(), 0, DoF->GetNumDim(), 1) += VectorXd::Ones(DoF->GetNumDim())*massPerNode;
			}
		}
		else
		{
			for (size_t i = 0; i < numNode; ++i)
			{
				DoFSet* DoF = this->m_pPolytope->Nodes()[i]->DoF();

				vMass.block(DoF->GetOffset_Full(), 0, DoF->GetNumDim(), 1) += VectorXd::Ones(DoF->GetNumDim())*massPerNode;
			}
		}
	}
	
	void MassElement_Lumped::AssembleGlobal_MassMatrix(VectorTd& vMass, bool full)
	{
		size_t numNode = this->m_pPolytope->Nodes().size();
		Real massPerNode = this->m_mass / numNode;

		if (!full)
		{
			for (size_t i = 0; i < numNode; ++i)
			{
				DoFSet* DoF = this->m_pPolytope->Nodes()[i]->DoF();

				// Is this a simulated DoF?
				if (DoF->GetId_Free() < 0)
				{
					continue;
				}

				for (int ii = 0; ii < DoF->GetNumDim(); ++ii)
					vMass.push_back(Triplet<Real>(
						DoF->GetOffset_Free() + ii,
						DoF->GetOffset_Free() + ii,
						massPerNode));
			}
		}
		else
		{
			for (size_t i = 0; i < numNode; ++i)
			{
				DoFSet* DoF = this->m_pPolytope->Nodes()[i]->DoF();

				for (int ii = 0; ii < DoF->GetNumDim(); ++ii)
					vMass.push_back(Triplet<Real>(
						DoF->GetOffset_Full() + ii,
						DoF->GetOffset_Full() + ii,
						massPerNode));
			}
		}
	}
}