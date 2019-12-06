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

#include <PhySim/Models/Model.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	/**
	* Model_Particles.
	*
	* TODO.
	*/
	class Model_Particles : public Model
	{
	public:

		/**
		* TODO.
		*/
		struct Options
		{
			/**
			* (Nn x 3) matrix with the node coordinates.
			*/
			MatrixXd m_mNodes;

			/**
			* (Ne x Ni) matrix with discretization elements.
			*/
			MatrixXi m_mElems;

			/**
			* Whether or not this particle system has sub-elements.
			*/

			bool m_hasSubElem;

			/**
			* (Nn_sub x 3) matrix with the sub-element node coordinates. This
			* are only used by some particle systems where the elements are 
			* subdivided for detailed computations, e.g. quadratic Tets. FEM.
			*/
			MatrixXd m_mNodesSub;

			/**
			* (Nn_sub x 3) matrix with discretization sub-elements. This
			* are only used by some particle systems where the elements are 
			* subdivided for detailed computations, e.g. quadratic Tets. FEM.
			*/
			MatrixXi m_mElemsSub;

			/**
			* Material properties
			*/
			Material m_material;

			/*
			* Gravity acceleration
			*/
			Vector3d m_vgravity;

			/**
			* Type of discretization
			*/
			Discretization m_discretization;

			Options()
			{
				m_material.AddProperty(Material::Property::Density, 1);
				m_discretization = Discretization::Nodes;
				m_vgravity = Vector3d::Zero();
				m_hasSubElem = false;
			}
		};

	protected:

		Options* m_pOptions;
		vector<Node*> m_vnodes;
		vector<Polytope*> m_velems;

	public:
		
		virtual string GetName() const override { return "Particle System"; }

		virtual Options& GetOptions();

		virtual void Init() override;
		virtual void Free() override;

	public:

		/**
		* Constructor.
		*/
		Model_Particles();

		/**
		* Denstructor.
		*/
		virtual ~Model_Particles();


		inline virtual const vector<Node*>& GetNodes() const;
		inline virtual const vector<Polytope*>& GetElements() const;

		virtual void InitializeRest();

		virtual void ComputeAndStore_Energy() override;
		virtual void ComputeAndStore_Mass(bool full) override;
		virtual void ComputeAndStore_Gradient(bool full) override;
		virtual void ComputeAndStore_Hessian(bool full) override;

		virtual void SetSubelementPositions(const MatrixXd& mV, Space s = Space::DEF);

		virtual vector<DoFSet*> SelectDoF(const Vector3d& vboxMin, const Vector3d& vboxMax, Space s = Space::DEF) const override;
		
		virtual vector<NodeEBD> ComputeEmbedding(const vector<Vector3d>& vp, Space s = Space::DEF);

	protected:

		virtual void InitDiscretization();
		virtual void InitDegreesOfFreedom();
		virtual void InitEnergyElements();
		virtual void InitMassElements();

	};
}