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

#include <PhySim/Models/Model_Particles.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	/**
	* Model_FEM.
	*
	* TODO.
	*/
	class Model_FEM : public Model_Particles
	{
	public:

		/**
		* TODO.
		*/
		struct Options : public Model_Particles::Options
		{
			/**
			* Constitutive model: co-rotational, St.VK, Neo-Hookean, Ogden, etc.
			*/
			MaterialModel m_materialModel;

			int m_numQuadrature;

			Options()
			{
				m_material.InitRealisticFromYoungPoisson(1e6, 0.45, 1e3);
				m_discretization = Discretization::Tetrahedra4;
				m_materialModel = MaterialModel::StVK;
				m_numQuadrature = 1;
			}
		};

	protected:

		vector<EnergyElement_FEM*> venergyEle_fem;
		vector<Material> vmat_fem;

	public:

		virtual string GetName() const override { return "FEM System"; }

		virtual Options& GetOptions() override;

		/**
		* Constructor.
		*/
		Model_FEM();

		/**
		* Destructor.
		*/
		virtual ~Model_FEM();

		// Initialization
		virtual void Free() override;

		// Parameter: rest strain
		virtual void GetRestStrains(vector<MatrixXd>& vE0) const;
		virtual void SetRestStrains(const vector<MatrixXd>& vE0);

		virtual const vector<Material>& GetMaterials_FEM() const { return this->vmat_fem; }
		virtual void SetMaterials_FEM(const vector<Material>& vM) { this->vmat_fem = vM; }

		virtual const vector<EnergyElement_FEM*>& GetEnergyElements_FEM() const { return this->venergyEle_fem; }

	protected:

		virtual void InitDiscretization() override;
		virtual void InitEnergyElements() override;

	};
}

