//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Models/Model_MassSpring_SMP.h>

#include <PhySim/Energies/EnergyElement_FaceNodeColl_Step.h>
#include <PhySim/Energies/EnergyElement_FaceNodeColl_Quad.h>
#include <PhySim/Energies/EnergyElement_FaceNodeColl_LogB.h>

#include <PhySim/Energies/EnergyElement_FEM_MembraneCoNH.h>
#include <PhySim/Energies/EnergyElement_FEM_MembraneInNH.h>
#include <PhySim/Energies/EnergyElement_FEM_MembraneStVK.h>

#include <PhySim/Geometry/Node.h>
#include <PhySim/Geometry/Edge.h>

#include <PhySim/Geometry/Face_Tri.h>

#include <PhySim/Energies/Material.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	Model_MassSpring_SMP::Model_MassSpring_SMP() : Model_MassSpring()
	{
		// Nothing to do here...
	}

	Model_MassSpring_SMP::~Model_MassSpring_SMP()
	{
		// Nothing to do here...
	}

	void Model_MassSpring_SMP::Free()
	{
		Model_MassSpring::Free();

		this->m_venergyEle_coll.clear();
		this->m_venergyEle_mem.clear();
		this->m_vmat_coll.clear();
		this->m_vmat_mem.clear();

		for (size_t i = 0; i < this->m_velems_mem.size(); ++i)
			delete m_velems_mem[i]; // Delete membrane faces
		this->m_velems_mem.clear();
	}

	Model_MassSpring_SMP::Options& Model_MassSpring_SMP::GetOptions()
	{
		if (this->m_pOptions == NULL) // Create if needed
			this->m_pOptions = new Model_MassSpring_SMP::Options();
		return *((Model_MassSpring_SMP::Options*) this->m_pOptions);
	}

	void Model_MassSpring_SMP::InitEnergyElements()
	{
		Model_MassSpring::InitEnergyElements();

		Options* pOptions = static_cast<Options*>(this->m_pOptions);

		size_t curEles = m_venergyEle.size();

		// Add collision elements -------------------------

		size_t numColl = pOptions->m_mCollPairs.rows();
		assert(pOptions->m_vCollDists.rows() == numColl);

		if (numColl > 0)
		{
			this->m_venergyEle.reserve(curEles + numColl);
			this->m_venergyEle_coll.resize(numColl);
			this->m_vmat_coll.resize(numColl);
			curEles += numColl;

			for (int i = 0; i < numColl; ++i)
			{
				Real D = pOptions->m_vCollDists[i];
				this->m_vmat_coll[i].AddProperty(Material::Property::CollK, pOptions->m_material[Material::Property::CollK]);
				this->m_vmat_coll[i].AddProperty(Material::Property::CollT, pOptions->m_material[Material::Property::CollT] + D);

				vector<Node*> vcollNodes;
				for (int j = 0; j < 4; ++j)
					vcollNodes.push_back(m_vnodes[pOptions->m_mCollPairs(i, j)]);

				EnergyElement_FaceNodeColl_Quad* pEle = new EnergyElement_FaceNodeColl_Quad(this, vcollNodes, &this->m_vmat_coll[i]);
				this->m_venergyEle.push_back(pEle);
				this->m_venergyEle_coll[i] = pEle;
			}
		}

		// Add membrane elements -------------------------

		// Create elements

		size_t numMemEles = pOptions->m_mmemElems.rows();

		if (numMemEles > 0)
		{
			// Check material

			switch (pOptions->m_memMaterialModel)
			{
			case MaterialModel::StVK:
			case MaterialModel::CoNH:
			case MaterialModel::InNH:
				assert(pOptions->m_material.HasProperty(Material::Property::Lame1));
				assert(pOptions->m_material.HasProperty(Material::Property::Lame2));
				break;
			default:
            assert(false); ;
			}

			switch (pOptions->m_memDiscretization)
			{
			case Discretization::Triangles:
			case Discretization::Quadrangles:
				assert(pOptions->m_material.HasProperty(Material::Property::Thickness));
				break;
			default:
            assert(false); ;
			}

			// Add elements

			this->m_venergyEle.reserve(curEles + numMemEles);
			this->m_venergyEle_mem.resize(numMemEles);
			this->m_velems_mem.reserve(numMemEles);
			this->m_vmat_mem.resize(numMemEles);

			for (size_t i = 0; i < numMemEles; ++i)
			{
				this->m_vmat_mem[i].AddProperty(Material::Property::Thickness, pOptions->m_material[Material::Property::Thickness]);

				switch (pOptions->m_memMaterialModel)
				{
				case MaterialModel::StVK:
				case MaterialModel::CoNH:
				case MaterialModel::InNH:
					this->m_vmat_mem[i].AddProperty(Material::Property::Lame1, pOptions->m_material[Material::Property::Lame1]);
					this->m_vmat_mem[i].AddProperty(Material::Property::Lame2, pOptions->m_material[Material::Property::Lame2]);
					break;
				default:
               assert(false); ;
				}

				EnergyElement_FEM* pEle;
				vector<Node*> vnodes(3);
				Face* pFace;

				switch (pOptions->m_memDiscretization)
				{
				case Discretization::Triangles:
					vnodes[0] = this->m_vnodes[pOptions->m_mmemElems(i, 0)];
					vnodes[1] = this->m_vnodes[pOptions->m_mmemElems(i, 1)];
					vnodes[2] = this->m_vnodes[pOptions->m_mmemElems(i, 2)];
					pFace = new Face_Tri(vnodes);

					switch (pOptions->m_memMaterialModel)
					{
					case MaterialModel::StVK: pEle = new EnergyElement_FEM_MembraneStVK(this, pFace, &this->m_vmat_mem[i]); break;
					case MaterialModel::CoNH: pEle = new EnergyElement_FEM_MembraneCoNH(this, pFace, &this->m_vmat_mem[i]); break;
					case MaterialModel::InNH: pEle = new EnergyElement_FEM_MembraneInNH(this, pFace, &this->m_vmat_mem[i]); break;
					default:
                  assert(false); ;
					}
					break;
				default:
               assert(false); ;
				}

				this->m_venergyEle.push_back(pEle);
				this->m_venergyEle_mem[i] = pEle;
			}
		}
	}

	void Model_MassSpring_SMP::GetMembraneRestStrains(vector<MatrixXd>& vE0) const
	{
		size_t numFem = this->m_venergyEle_mem.size();

		vE0.resize(numFem);

		for (int i = 0; i < numFem; ++i)
		{
			vE0[i] = this->m_venergyEle_mem[i]->GetRestStrain();
		}
	}

	void Model_MassSpring_SMP::SetMembraneRestStrains(const vector<MatrixXd>& vE0)
	{
		size_t numFem = this->m_venergyEle_mem.size();

		assert(vE0.size() == numFem);

		for (int i = 0; i < numFem; ++i)
		{
			this->m_venergyEle_mem[i]->SetRestStrain(vE0[i]);

			this->m_venergyEle_mem[i]->Init();
		}

		this->DirtyDeformed();
	}

}
