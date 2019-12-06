//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#pragma once

#include <PhySim/CommonIncludes.h>
#include <PhySim/PhySimInterface.h>

#include <PhySim/Models/Model_FEM.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	/**
	* Model_ThinShells.
	*
	* TODO.
	*/
	class Model_ThinShells : public Model_FEM
	{
	public:

		/**
		* TODO.
		*/
		struct Options : public Model_FEM::Options
		{
			Options()
			{
				m_material.InitRealisticFromYoungPoisson(1e6, 0.45, 1e3);
				m_discretization = Discretization::Triangles;
				m_materialModel = MaterialModel::StVK;
			}
		};

	protected:

		vector<EnergyElement_ShellHinge*> venergyEle_hinges;
		vector<Material> vmat_hinges;

	public:

		virtual string GetName() const override { return "Thin-Shells System"; }

		virtual Options& GetOptions() override;

		/**
		* Constructor.
		*/
		Model_ThinShells();

		/**
		* Destructor.
		*/
		virtual ~Model_ThinShells();

		// Initialization
		virtual void Free() override;

		virtual const vector<Material>& GetMaterials_Hinges() const { return this->vmat_hinges; }
		virtual void SetMaterials_Hinges(const vector<Material>& vM) { this->vmat_hinges = vM; }

		virtual const vector<EnergyElement_ShellHinge*>& GetEnergyElements_ShellHinge() const { return this->venergyEle_hinges; }

	protected:

		virtual void InitEnergyElements() override;

	};
}

