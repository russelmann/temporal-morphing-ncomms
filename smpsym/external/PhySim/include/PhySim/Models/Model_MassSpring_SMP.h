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

#include <PhySim/Models/Model_MassSpring.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	/**
	* Model_MassSpring_SMP.
	*
	* TODO.
	*/
	class Model_MassSpring_SMP : public Model_MassSpring
	{
	public:

		/**
		* TODO
		*/
		struct Options : public Model_MassSpring::Options
		{
			/**
			* (Nf x 4) matrix with node indices corresponding to face-node pairs to collide.
			*/
			MatrixXi m_mCollPairs;

			/**
			* Vectors of distances for each of the colliding pairs to simulate bumpers volume.
			*/
			VectorXd m_vCollDists;

			/**
			* Membrane material model.
			*/
			MaterialModel m_memMaterialModel;
			Discretization m_memDiscretization;
			MatrixXi m_mmemElems;

			Options() : Model_MassSpring::Options()
			{
				m_material.InitRealisticFromYoungPoisson(1e-6, 0.45, 1);

				m_discretization = Discretization::Edges;
				m_material.AddProperty(Material::Property::CollK, 1e3);
				m_material.AddProperty(Material::Property::CollT, 0.1);
			}
		};

	protected:

		vector<EnergyElement_FaceNodeColl_Quad*> m_venergyEle_coll;
		vector<EnergyElement_FEM*> m_venergyEle_mem;
		vector<Face*> m_velems_mem;
		vector<Material> m_vmat_coll;
		vector<Material> m_vmat_mem;

	public:

		virtual string GetName() const { return "Mass-Spring SMP System"; }
	
		virtual Options& GetOptions();

		/**
		* Constructor.
		*/
		Model_MassSpring_SMP();

		/**
		* Destructor.
		*/
		virtual ~Model_MassSpring_SMP();

		// Initialization
		virtual void Free() override;

		// Parameter: membrane rest strain
		virtual void GetMembraneRestStrains(vector<MatrixXd>& vE0) const;
		virtual void SetMembraneRestStrains(const vector<MatrixXd>& vE0);

		virtual const vector<EnergyElement_FaceNodeColl_Quad*>& GetEnergyElements_Collision() const { return this->m_venergyEle_coll; }
		virtual const vector<EnergyElement_FEM*>& GetEnergyElements_Membrane() const { return this->m_venergyEle_mem; }

	protected:
		virtual void InitEnergyElements();

	};
}

