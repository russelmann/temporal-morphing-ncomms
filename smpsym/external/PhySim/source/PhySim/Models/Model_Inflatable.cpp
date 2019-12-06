//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Models/Model_Inflatable.h>

#include <PhySim/Energies/EnergyElement_PressureConstant.h>

#include <PhySim/Geometry/Face_Tri.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	Model_Inflatable::Model_Inflatable() : Model_FEM()
	{
		// Nothing to do here...
	}

	Model_Inflatable::~Model_Inflatable()
	{
		// Nothing to do here...
	}

	void Model_Inflatable::Free()
	{
		Model_FEM::Free();

		for (size_t i = 0; i < this->m_vintFaces.size(); ++i)
			delete this->m_vintFaces[i];
		this->m_vintFaces.clear();

		for (size_t i = 0; i < this->m_vextFaces.size(); ++i)
			delete this->m_vextFaces[i];
		this->m_vextFaces.clear();

		this->m_venergyEle_intPressure.clear();
		this->m_venergyEle_extPressure.clear();
	}

	Model_Inflatable::Options& Model_Inflatable::GetOptions()
	{
		if (this->m_pOptions == NULL) // Create if needed
			this->m_pOptions = new Model_Inflatable::Options();
		return *((Model_Inflatable::Options*) this->m_pOptions);
	}

	void Model_Inflatable::InitEnergyElements()
	{
		// Check possible discretizations

		assert(m_pOptions->m_discretization == Discretization::Triangles ||
			   m_pOptions->m_discretization == Discretization::Tetrahedra4);

		// Initialize FEM elements

		Model_FEM::InitEnergyElements();

		// Add the pressure materials

		Options* pOptions = static_cast<Options*>(this->m_pOptions);
		this->m_mat_intPressure.AddProperty(Material::Property::Pressure, pOptions->m_internalPressure);
		this->m_mat_extPressure.AddProperty(Material::Property::Pressure, pOptions->m_externalPressure);

		// Add pressure elements

		int numInt = (int) pOptions->m_vintFaces.rows();
		int numExt = (int) pOptions->m_vextFaces.rows();

		this->m_venergyEle.reserve(this->m_venergyEle.size() + numInt + numExt);
		this->m_venergyEle_intPressure.reserve(numInt);
		this->m_venergyEle_extPressure.reserve(numExt);

		for (int i = 0; i < numInt; ++i)
		{
			vector<Node*> vnodes(3);
			vnodes[0] = this->m_vnodes[pOptions->m_vintFaces(i, 0)];
			vnodes[1] = this->m_vnodes[pOptions->m_vintFaces(i, 1)];
			vnodes[2] = this->m_vnodes[pOptions->m_vintFaces(i, 2)];
			Face_Tri* pFace = new Face_Tri(vnodes);
			this->m_venergyEle_intPressure.push_back(new EnergyElement_PressureConstant(this, pFace, &m_mat_intPressure));
			this->m_venergyEle.push_back(this->m_venergyEle_intPressure.back());
		}

		for (int i = 0; i < numExt; ++i)
		{
			vector<Node*> vnodes(3);
			vnodes[0] = this->m_vnodes[pOptions->m_vextFaces(i, 0)];
			vnodes[1] = this->m_vnodes[pOptions->m_vextFaces(i, 1)];
			vnodes[2] = this->m_vnodes[pOptions->m_vextFaces(i, 2)];
			Face_Tri* pFace = new Face_Tri(vnodes);
			this->m_venergyEle_extPressure.push_back(new EnergyElement_PressureConstant(this, pFace, &m_mat_extPressure));
			this->m_venergyEle.push_back(this->m_venergyEle_extPressure.back());
		}
	}

	void Model_Inflatable::MultInternalPressure(Real factor)
	{
		Real pressure = this->m_mat_intPressure[Material::Property::Pressure];
		this->m_mat_intPressure[Material::Property::Pressure] = factor*pressure;

		this->DirtyDeformed();
	}

	void Model_Inflatable::SumInternalPressure(Real value)
	{
		Real pressure = this->m_mat_intPressure[Material::Property::Pressure];
		this->m_mat_intPressure[Material::Property::Pressure] = value+pressure;

		this->DirtyDeformed();
	}

}

