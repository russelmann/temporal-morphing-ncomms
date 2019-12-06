//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Models/Model_ThinShells.h>

#include <PhySim/Energies/MassElement_Lumped.h>
#include <PhySim/Energies/EnergyElement_FEM.h>
#include <PhySim/Energies/EnergyElement_FEM_VolumeStVK.h>
#include <PhySim/Energies/EnergyElement_FEM_VolumeCoNH.h>
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

	Model_ThinShells::Model_ThinShells() : Model_FEM()
	{
		// Nothing to do here...
	}

	Model_ThinShells::~Model_ThinShells()
	{
		// Nothing to do here...
	}

	void Model_ThinShells::Free()
	{
		Model_FEM::Free();

		this->venergyEle_hinges.clear();
		this->vmat_hinges.clear();
	}

	Model_ThinShells::Options& Model_ThinShells::GetOptions()
	{
		if (this->m_pOptions == NULL) // Create if needed
			this->m_pOptions = new Model_ThinShells::Options();
		return *((Model_ThinShells::Options*) this->m_pOptions);
	}

	void Model_ThinShells::InitEnergyElements()
	{
		// Check possible discretizations

		assert(m_pOptions->m_discretization == Discretization::Triangles ||
			   m_pOptions->m_discretization == Discretization::Quadrangles);

		// Initialize FEM elements

		Model_FEM::InitEnergyElements();

		// Check bending material

		assert(this->m_pOptions->m_material.HasProperty(Material::Property::BendingK));

		// Add the bending elements

		Options* pOptions = static_cast<Options*>(this->m_pOptions);

		// TODO
	}

}

