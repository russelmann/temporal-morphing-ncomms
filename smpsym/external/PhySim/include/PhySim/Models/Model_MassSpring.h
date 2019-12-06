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
	* Model_MassSpring.
	*
	* TODO.
	*/
	class Model_MassSpring : public Model_Particles
	{
	public:

		enum SpringType
		{
			QuadK,
			UPoly8,
			MPoly8
		};

		/**
		* TODO.
		*/
		struct Springs
		{
			VectorXi m_vlinear;
			MatrixXi m_mhinges;
			MatrixXi m_mcrosses;
		
			Springs()
			{
				m_vlinear.setZero();
				m_mhinges.setZero();
				m_mcrosses.setZero();
			}
		};

		/**
		* TODO
		*/
		struct Options : public Model_Particles::Options
		{
			/**
			* Structure of springs where vlinear, mhinges and mcrosses contain indices refering to mElems rows.
			*/
			Springs m_springs;
			SpringType m_linearType;
			SpringType m_hingesType;
			SpringType m_crossType;

			Options()
			{
				m_discretization = Discretization::Edges;
				m_material.AddProperty(Material::Property::Density, 1);
				m_material.AddProperty(Material::Property::StretchK, 1);
				m_material.AddProperty(Material::Property::BendingK, 0);
				m_material.AddProperty(Material::Property::ShearK, 0);
				m_linearType = SpringType::QuadK;
				m_hingesType = SpringType::QuadK;
				m_crossType = SpringType::QuadK;
			}
		};

	protected:

		vector<EnergyElement_SpringLinear*> m_venergyEle_linear;
		vector<EnergyElement_SpringHinge*> m_venergyEle_hinge;
		vector<EnergyElement_SpringCross*> m_venergyEle_cross;
		vector<Material> m_vmat_linear;
		vector<Material> m_vmat_hinge;
		vector<Material> m_vmat_cross;

	public:

		virtual string GetName() const override { return "Mass-Spring System"; }

		virtual Options& GetOptions() override;

		/**
		* Constructor.
		*/
		Model_MassSpring();

		/**
		* Destructor.
		*/
		virtual ~Model_MassSpring();

		// Initialization
		virtual void Free();
	
		// Parameter: rest length
		virtual VectorXd GetRestLinearLengths() const;
		virtual void SetRestLinearLengths(const VectorXd&);

		// Parameter: rest hinge angles
		virtual VectorXd GetRestHingeAngles() const;
		virtual void SetRestHingeAngles(const VectorXd&);

		// Parameter: rest cross ratios
		virtual VectorXd GetRestCrossStrain() const;
		virtual void SetRestCrossStrain(const VectorXd&);

		// Parameter: stretch polynomial coefficients
		virtual void GetStretchPolyCoeff(vector<VectorXd>& vcoeff) const;
		virtual void SetStretchPolyCoeff(const vector<VectorXd>& vcoeff);

		// Parameter: stretch K
		virtual VectorXd GetStretchK() const;
		virtual void SetStretchK(const VectorXd&);

		// Parameter: bending K
		virtual VectorXd GetBendingK() const;
		virtual void SetBendingK(const VectorXd&);

		// Parameter: shear K
		virtual VectorXd GetShearK() const;
		virtual void SetShearK(const VectorXd&);

		virtual const vector<EnergyElement_SpringLinear*>& GetEnergyElements_SpringsLinear() const { return this->m_venergyEle_linear; }
		virtual const vector<EnergyElement_SpringHinge*>& GetEnergyElements_SpringsHinge() const { return this->m_venergyEle_hinge; }
		virtual const vector<EnergyElement_SpringCross*>& GetEnergyElements_SpringsCross() const { return this->m_venergyEle_cross; }

	protected:
		virtual void InitEnergyElements() override;

	};
}

