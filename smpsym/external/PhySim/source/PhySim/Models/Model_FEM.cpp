//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Models/Model_FEM.h>

#include <PhySim/Energies/MassElement_Lumped.h>
#include <PhySim/Energies/EnergyElement_FEM.h>
#include <PhySim/Energies/EnergyElement_FEM_VolumeStVK.h>
#include <PhySim/Energies/EnergyElement_FEM_VolumeCoNH.h>
#include <PhySim/Energies/EnergyElement_FEM_VolumeCoMR.h>
#include <PhySim/Energies/EnergyElement_FEM_MembraneStVK.h>
#include <PhySim/Energies/EnergyElement_FEM_MembraneCoNH.h>
#include <PhySim/Energies/EnergyElement_Gravity.h>

#include <PhySim/Geometry/Node.h>
#include <PhySim/Geometry/Edge.h>
#include <PhySim/Geometry/Face_Tri.h>
#include <PhySim/Geometry/Face_Quad.h>
#include <PhySim/Geometry/Cell_Tetra.h>
#include <PhySim/Geometry/Cell_Hexa.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	Model_FEM::Model_FEM() : Model_Particles()
	{
		// Nothing to do here...
	}

	Model_FEM::~Model_FEM()
	{
		// Nothing to do...
	}

	void Model_FEM::Free()
	{
		Model_Particles::Free();

		this->venergyEle_fem.clear();
		this->vmat_fem.clear();
	}

	Model_FEM::Options& Model_FEM::GetOptions()
	{
		if (this->m_pOptions == NULL) // Create if needed
			this->m_pOptions = new Model_FEM::Options();
		return *((Model_FEM::Options*) this->m_pOptions);
	}

	void Model_FEM::InitDiscretization()
	{
		Options* pOptions = static_cast<Options*>(this->m_pOptions);

		if (pOptions->m_discretization == Discretization::Tetrahedra10)
		{
			pOptions->m_hasSubElem = true;

			// Create additional nodes and elements

			const MatrixXd& mV = pOptions->m_mNodes;
			const MatrixXi& mE = pOptions->m_mElems;

			map<pair<int, int>, Vector3d> mnewNodes;

			for (int i = 0; i < mE.rows(); ++i)
			{
				for (int ii = 0; ii < 4; ++ii)
					for (int jj = ii + 1; jj < 4; ++jj)
					{
						int node0 = mE(i, ii);
						int node1 = mE(i, jj);
						if (node0 > node1)
						{
							int temp = node0;
							node0 = node1;
							node1 = temp;
						}
						if (mnewNodes.find(pair<int,int>(node0,node1)) == mnewNodes.end())
						{
							mnewNodes.insert(pair<pair<int, int>, Vector3d>(pair<int, int>(node0, node1), 0.5*(mV.row(node0) + mV.row(node1))));
						}
					}
			}

			int numNew = (int) mnewNodes.size();

			// Build new node matrix
			 
			pOptions->m_mNodesSub = MatrixXd(mV.rows() + numNew, 3);
			pOptions->m_mNodesSub.block(0, 0, mV.rows(), 3) = mV;

			map<pair<int, int>, int> mnewIndices;
			map<pair<int, int>, Vector3d>::iterator itCur = mnewNodes.begin();
			map<pair<int, int>, Vector3d>::iterator itEnd = mnewNodes.end();
			int countNode = 0;
			for (; itCur != itEnd; ++itCur)
			{
				int newIndex = (int) mV.rows() + countNode++;
				pOptions->m_mNodesSub.row(newIndex) = itCur->second;
				mnewIndices.insert(pair<pair<int,int>, int>(itCur->first, newIndex));
			}
				
			// Build new elements matrix

			map<pair<int, int>, int> subNodeOrder;
			subNodeOrder.insert(pair<pair<int, int>, int>(pair<int, int>(0, 1), 4));
			subNodeOrder.insert(pair<pair<int, int>, int>(pair<int, int>(1, 2), 5));
			subNodeOrder.insert(pair<pair<int, int>, int>(pair<int, int>(0, 2), 6));
			subNodeOrder.insert(pair<pair<int, int>, int>(pair<int, int>(0, 3), 7));
			subNodeOrder.insert(pair<pair<int, int>, int>(pair<int, int>(1, 3), 8));
			subNodeOrder.insert(pair<pair<int, int>, int>(pair<int, int>(2, 3), 9));

			pOptions->m_mElemsSub = MatrixXi(mE.rows(), 10); // Tet10
			pOptions->m_mElemsSub.block(0, 0, mE.rows(), 4) = mE;

			for (int i = 0; i < mE.rows(); ++i)
			{
				for (int ii = 0; ii < 4; ++ii)
					for (int jj = ii + 1; jj < 4; ++jj)
					{
						int node0 = mE(i, ii);
						int node1 = mE(i, jj);
						if (node0 > node1)
						{
							int temp = node0;
							node0 = node1;
							node1 = temp;
						}

						pOptions->m_mElemsSub(i, subNodeOrder[pair<int, int>(ii, jj)]) = mnewIndices[pair<int, int>(node0, node1)];
					}
			}
		}

		Model_Particles::InitDiscretization();
	}
	 
	void Model_FEM::InitEnergyElements()
	{
		Model_Particles::InitEnergyElements();

		Options* pOptions = static_cast<Options*>(this->m_pOptions);

		// Check material

		switch (pOptions->m_materialModel)
		{
		case MaterialModel::CoRot:
			assert(this->m_pOptions->m_material.HasProperty(Material::Property::Young));
			assert(this->m_pOptions->m_material.HasProperty(Material::Property::Poisson));
			break;
		case MaterialModel::StVK:
		case MaterialModel::CoNH: 
			assert(this->m_pOptions->m_material.HasProperty(Material::Property::Lame1));
			assert(this->m_pOptions->m_material.HasProperty(Material::Property::Lame2));
			break;
		case MaterialModel::CoMR:
			assert(this->m_pOptions->m_material.HasProperty(Material::Property::Mooney01));
			assert(this->m_pOptions->m_material.HasProperty(Material::Property::Mooney10));
			assert(this->m_pOptions->m_material.HasProperty(Material::Property::Lame2));
			break;
		default:
         assert(false); ;
		}

		size_t numElems = this->m_velems.size();
		size_t curElems = this->m_venergyEle.size();
		this->m_venergyEle.reserve(curElems + numElems);
		this->venergyEle_fem.resize(numElems);
		this->vmat_fem.resize(numElems);

		for (size_t i = 0; i < numElems; ++i)
		{
			switch (pOptions->m_materialModel)
			{
			case MaterialModel::CoRot:
				this->vmat_fem[i].AddProperty(Material::Property::Young, this->m_pOptions->m_material[Material::Property::Young]);
				this->vmat_fem[i].AddProperty(Material::Property::Poisson, this->m_pOptions->m_material[Material::Property::Poisson]);
				break;
			case MaterialModel::StVK:
			case MaterialModel::CoNH:
				this->vmat_fem[i].AddProperty(Material::Property::Lame1, this->m_pOptions->m_material[Material::Property::Lame1]);
				this->vmat_fem[i].AddProperty(Material::Property::Lame2, this->m_pOptions->m_material[Material::Property::Lame2]);
				break;
			case MaterialModel::CoMR:
				this->vmat_fem[i].AddProperty(Material::Property::Mooney01, this->m_pOptions->m_material[Material::Property::Mooney01]);
				this->vmat_fem[i].AddProperty(Material::Property::Mooney10, this->m_pOptions->m_material[Material::Property::Mooney10]);
				this->vmat_fem[i].AddProperty(Material::Property::Lame2, this->m_pOptions->m_material[Material::Property::Lame2]);
				break;
			default:
            assert(false); ;
			}

			EnergyElement_FEM* pEle = NULL;
			switch (pOptions->m_discretization)
			{
			case Discretization::Triangles:
			case Discretization::Quadrangles:
				switch (pOptions->m_materialModel)
				{
				case MaterialModel::StVK: pEle = new EnergyElement_FEM_MembraneStVK(this, (Face*) this->m_velems[i], &this->vmat_fem[i]); break;
				case MaterialModel::CoNH: pEle = new EnergyElement_FEM_MembraneCoNH(this, (Face*) this->m_velems[i], &this->vmat_fem[i]); break;
				default:
               assert(false); ;
				}
				break;
			case Discretization::Tetrahedra4:
				
				assert(pOptions->m_numQuadrature == 1);

				switch (pOptions->m_materialModel)
				{
				case MaterialModel::StVK: pEle = new EnergyElement_FEM_VolumeStVK(this, (Cell*) this->m_velems[i], &this->vmat_fem[i], 1); break;
				case MaterialModel::CoNH: pEle = new EnergyElement_FEM_VolumeCoNH(this, (Cell*) this->m_velems[i], &this->vmat_fem[i], 1); break;
				case MaterialModel::CoMR: pEle = new EnergyElement_FEM_VolumeCoMR(this, (Cell*) this->m_velems[i], &this->vmat_fem[i], 1); break;
				default:
               assert(false); ;
				}
				break;
			case Discretization::Tetrahedra10:

				assert(pOptions->m_numQuadrature == 4 ||
					   pOptions->m_numQuadrature == 5 ||
					   pOptions->m_numQuadrature == 8 ||
					   pOptions->m_numQuadrature == 11);

				switch (pOptions->m_materialModel)
				{
				case MaterialModel::StVK: pEle = new EnergyElement_FEM_VolumeStVK(this, (Cell*) this->m_velems[i], &this->vmat_fem[i], pOptions->m_numQuadrature); break;
				case MaterialModel::CoNH: pEle = new EnergyElement_FEM_VolumeCoNH(this, (Cell*) this->m_velems[i], &this->vmat_fem[i], pOptions->m_numQuadrature); break;
				case MaterialModel::CoMR: pEle = new EnergyElement_FEM_VolumeCoMR(this, (Cell*) this->m_velems[i], &this->vmat_fem[i], pOptions->m_numQuadrature); break;
				default:
               assert(false); ;
				}
				break;
			case Discretization::Hexahedra:
				switch (pOptions->m_materialModel)
				{
				case MaterialModel::StVK: pEle = new EnergyElement_FEM_VolumeStVK(this, (Cell*) this->m_velems[i], &this->vmat_fem[i], pOptions->m_numQuadrature); break;
				case MaterialModel::CoNH: pEle = new EnergyElement_FEM_VolumeCoNH(this, (Cell*) this->m_velems[i], &this->vmat_fem[i], pOptions->m_numQuadrature); break;
				default:
               assert(false); ;
				}
				break;
			default:
            assert(false); ;
			}
			this->m_venergyEle.push_back(pEle);
			this->venergyEle_fem[i] = pEle;
		}
	}

	void Model_FEM::GetRestStrains(vector<MatrixXd>& vE0) const
	{
		size_t numFem = this->venergyEle_fem.size();

		vE0.resize(numFem);

		for (int i = 0; i < numFem; ++i)
		{
			vE0[i] = this->venergyEle_fem[i]->GetRestStrain();
		}
	}

	void Model_FEM::SetRestStrains(const vector<MatrixXd>& vE0)
	{
		size_t numFem = this->venergyEle_fem.size();

		assert(vE0.size() == numFem);

		for (int i = 0; i < numFem; ++i)
		{
			this->venergyEle_fem[i]->SetRestStrain(vE0[i]);
		}

		this->DirtyDeformed();
	}

}

