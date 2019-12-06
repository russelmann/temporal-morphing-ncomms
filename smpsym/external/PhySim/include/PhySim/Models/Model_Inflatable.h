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
	class Model_Inflatable : public Model_FEM
	{
	public:

		/**
		* TODO.
		*/
		struct Options : public Model_FEM::Options
		{
			MatrixXi m_vintFaces;
			MatrixXi m_vextFaces;
			Real m_internalPressure;
			Real m_externalPressure;

			Options()
			{
				m_material.InitRealisticFromYoungPoisson(1e6, 0.45, 1e3);
				m_discretization = Discretization::Tetrahedra4;
				m_materialModel = MaterialModel::StVK;

				m_internalPressure = 101325; // 1 atm
				m_externalPressure = 101325; // 1 atm
			}
		};

	protected:

		vector<Face_Tri*> m_vintFaces;
		vector<Face_Tri*> m_vextFaces;
		vector<EnergyElement_PressureConstant*> m_venergyEle_intPressure;
		vector<EnergyElement_PressureConstant*> m_venergyEle_extPressure;
		Material m_mat_intPressure;
		Material m_mat_extPressure;

	public:

		virtual string GetName() const override { return "Inflatable System"; }

		virtual Options& GetOptions() override;

		/**
		* Constructor.
		*/
		Model_Inflatable();

		/**
		* Destructor.
		*/
		virtual ~Model_Inflatable();

		// Initialization
		virtual void Free() override;

		virtual void MultInternalPressure(Real factor);
		virtual void SumInternalPressure(Real factor);

		virtual const Material& GetMaterials_IntPressure() const { return this->m_mat_intPressure; }
		virtual const Material& GetMaterials_ExtPressure() const { return this->m_mat_extPressure; }

		virtual const vector<EnergyElement_PressureConstant*>& GetEnergyElements_InternalPressure() const { return this->m_venergyEle_intPressure; }
		virtual const vector<EnergyElement_PressureConstant*>& GetEnergyElements_ExternalPressure() const { return this->m_venergyEle_extPressure; }

	protected:

		virtual void InitEnergyElements() override;

	};
}

