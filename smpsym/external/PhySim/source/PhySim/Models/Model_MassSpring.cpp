//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Models/Model_MassSpring.h>

#include <PhySim/Energies/MassElement_Lumped.h>
#include <PhySim/Energies/EnergyElement_SpringLinear.h>
#include <PhySim/Energies/EnergyElement_SpringLinear_UPoly.h>
#include <PhySim/Energies/EnergyElement_SpringLinear_MPoly.h>
#include <PhySim/Energies/EnergyElement_SpringHinge.h>
#include <PhySim/Energies/EnergyElement_SpringCross.h>
#include <PhySim/Energies/EnergyElement_Gravity.h>

#include <PhySim/Geometry/Node.h>
#include <PhySim/Geometry/Edge.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	Model_MassSpring::Model_MassSpring() : Model_Particles()
	{
		// Nothing to do here...
	}

	Model_MassSpring::~Model_MassSpring()
	{
		// Nothing to do here...
	}

	void Model_MassSpring::Free()
	{
		Model_Particles::Free();

		this->m_venergyEle_linear.clear();
		this->m_venergyEle_hinge.clear();
		this->m_venergyEle_cross.clear();
		m_vmat_linear.clear();
		m_vmat_hinge.clear();
		m_vmat_cross.clear();
	}

	Model_MassSpring::Options& Model_MassSpring::GetOptions()
	{
		if (this->m_pOptions == NULL) // Create if needed
			this->m_pOptions = new Model_MassSpring::Options();
		return *((Model_MassSpring::Options*) this->m_pOptions);
	}

	void Model_MassSpring::InitEnergyElements()
	{
		Model_Particles::InitEnergyElements();

		Options* m_pOptions = static_cast<Options*>(this->m_pOptions);

		const VectorXi& m_vlinear = m_pOptions->m_springs.m_vlinear;
		const MatrixXi& m_mhinges = m_pOptions->m_springs.m_mhinges;
		const MatrixXi& m_mcrosses = m_pOptions->m_springs.m_mcrosses;
		size_t numLinear = m_vlinear.size();
		size_t numCrosses = m_mcrosses.rows();
		size_t numHinges = m_mhinges.rows();

		this->m_venergyEle.reserve(this->m_venergyEle.size() + numLinear + numCrosses + numHinges);
		this->m_venergyEle_linear.resize(numLinear);
		this->m_venergyEle_cross.resize(numCrosses);
		this->m_venergyEle_hinge.resize(numHinges);

		if (numLinear > 0) { ; assert(m_pOptions->m_material.HasProperty(Material::Property::StretchK)); }
		if (numHinges > 0) { ; assert(m_pOptions->m_material.HasProperty(Material::Property::BendingK)); }
		if (numCrosses > 0) { ; assert(m_pOptions->m_material.HasProperty(Material::Property::ShearK)); }

		// Create energy elements (linear)

		m_vmat_linear.resize(numLinear);
		for (size_t i = 0; i < numLinear; ++i)
		{
			m_vmat_linear[i].AddProperty(Material::Property::Density, m_pOptions->m_material[Material::Property::Density]);

			EnergyElement_SpringLinear* pEle = NULL;

			if (m_pOptions->m_linearType == Model_MassSpring::SpringType::QuadK)
			{
				m_vmat_linear[i].AddProperty(Material::Property::StretchK, m_pOptions->m_material[Material::Property::StretchK]);
				pEle = new EnergyElement_SpringLinear(this, (Edge*) this->m_velems[m_vlinear(i)], &this->m_vmat_linear[i]);
			}

			if (m_pOptions->m_linearType == Model_MassSpring::SpringType::UPoly8)
			{
				pEle = new EnergyElement_SpringLinear_UPoly(this, (Edge*) this->m_velems[m_vlinear(i)], &this->m_vmat_linear[i]);
			}

			if (m_pOptions->m_linearType == Model_MassSpring::SpringType::MPoly8)
			{
				pEle = new EnergyElement_SpringLinear_MPoly(this, (Edge*) this->m_velems[m_vlinear(i)], &this->m_vmat_linear[i]);
			}

			this->m_venergyEle.push_back(pEle);
			this->m_venergyEle_linear[i] = pEle;
		}

		// Create energy elements (hinges)

		m_vmat_hinge.resize(numHinges);
		for (size_t i = 0; i < numHinges; ++i)
		{
			m_vmat_hinge[i].AddProperty(Material::Property::BendingK, m_pOptions->m_material[Material::Property::BendingK]);

			EnergyElement_SpringHinge* pEle = new EnergyElement_SpringHinge(this, (Edge*) this->m_velems[m_mhinges(i, 0)], (Edge*) this->m_velems[m_mhinges(i, 1)], &this->m_vmat_hinge[i]);
			this->m_venergyEle.push_back(pEle);
			this->m_venergyEle_hinge[i] = pEle;
		}

		// Create energy elements (cross)

		m_vmat_cross.resize(numCrosses);
		for (size_t i = 0; i < numCrosses; ++i)
		{
			m_vmat_cross[i].AddProperty(Material::Property::ShearK, m_pOptions->m_material[Material::Property::ShearK]);

			EnergyElement_SpringCross* pEle = new EnergyElement_SpringCross(this, (Edge*) this->m_velems[m_mcrosses(i, 0)], (Edge*) this->m_velems[m_mcrosses(i, 1)], &this->m_vmat_cross[i]);
			this->m_venergyEle.push_back(pEle);
			this->m_venergyEle_cross[i] = pEle;
		}
	}

	VectorXd Model_MassSpring::GetRestLinearLengths() const
	{
		size_t numLinear = this->m_venergyEle_linear.size();

		VectorXd vl(numLinear);

		for (size_t i = 0; i < numLinear; ++i)
		{
			vl[i] = this->m_venergyEle_linear[i]->GetRestLength();
		}

		return vl;
	}

	void Model_MassSpring::SetRestLinearLengths(const VectorXd& vl)
	{
		size_t numLinear = this->m_venergyEle_linear.size();

		assert(vl.size() == numLinear);

		for (size_t i = 0; i < numLinear; ++i)
		{
			this->m_venergyEle_linear[i]->SetRestLength(vl[i]);
		}

		this->DirtyDeformed();
	}

	VectorXd Model_MassSpring::GetRestHingeAngles() const
	{
		size_t numHinges = this->m_venergyEle_hinge.size();

		VectorXd vl(numHinges);

		for (size_t i = 0; i < numHinges; ++i)
		{
			vl[i] = this->m_venergyEle_hinge[i]->GetRestAngle();
		}

		return vl;
	}

	void Model_MassSpring::SetRestHingeAngles(const VectorXd& va)
	{
		size_t numHinges = this->m_venergyEle_hinge.size();

		assert(va.size() == numHinges);

		for (size_t i = 0; i < numHinges; ++i)
		{
			this->m_venergyEle_hinge[i]->SetRestAngle(va[i]);
		}

		this->DirtyDeformed();
	}

	VectorXd Model_MassSpring::GetRestCrossStrain() const
	{
		size_t numCrosses = this->m_venergyEle_cross.size();

		VectorXd vr(numCrosses);

		for (size_t i = 0; i < numCrosses; ++i)
		{
			vr[i] = this->m_venergyEle_cross[i]->GetRestStrain();
		}

		return vr;
	}

	void Model_MassSpring::SetRestCrossStrain(const VectorXd& vr)
	{
		size_t numCrosses = this->m_venergyEle_cross.size();

		assert(vr.size() == numCrosses);

		for (size_t i = 0; i < numCrosses; ++i)
		{
			this->m_venergyEle_cross[i]->GetRestStrain(vr[i]);
		}

		this->DirtyDeformed();
	}

	void Model_MassSpring::GetStretchPolyCoeff(vector<VectorXd>& vcoeff) const
	{
		size_t numLinear = this->m_venergyEle_linear.size();

		vcoeff.resize(numLinear);

		for (size_t i = 0; i < numLinear; ++i)
		{
			vcoeff[i] = static_cast<EnergyElement_SpringLinear_UPoly*>(this->m_venergyEle_linear[i])->GetPolyCoeff();
		}
	}

	VectorXd Model_MassSpring::GetStretchK() const
	{
		size_t numLinear = this->m_venergyEle_linear.size();

		VectorXd vs(numLinear);

		for (size_t i = 0; i < numLinear; ++i)
		{
			vs[i] = (*this->m_venergyEle_linear[i]->GetMaterial())[Material::Property::StretchK];
		}

		return vs;
	}

	VectorXd Model_MassSpring::GetBendingK() const
	{
		size_t numHinges = this->m_venergyEle_hinge.size();

		VectorXd vs(numHinges);

		for (size_t i = 0; i < numHinges; ++i)
		{
			vs[i] = (*this->m_venergyEle_hinge[i]->GetMaterial())[Material::Property::BendingK];
		}

		return vs;
	}

	VectorXd Model_MassSpring::GetShearK() const
	{
		size_t numCrosses = this->m_venergyEle_cross.size();

		VectorXd vs(numCrosses);

		for (size_t i = 0; i < numCrosses; ++i)
		{
			vs[i] = (*this->m_venergyEle_cross[i]->GetMaterial())[Material::Property::ShearK];
		}

		return vs;
	}

	void Model_MassSpring::SetStretchPolyCoeff(const vector<VectorXd>& vcoeff)
	{
		size_t numLinear = this->m_venergyEle_linear.size();

		assert(vcoeff.size() == numLinear);

		for (size_t i = 0; i < numLinear; ++i)
		{
			static_cast<EnergyElement_SpringLinear_UPoly*>(this->m_venergyEle_linear[i])->SetPolyCoeff(vcoeff[i]);
		}

		this->DirtyDeformed();
	}

	void Model_MassSpring::SetStretchK(const VectorXd& vs)
	{
		size_t numLinear = this->m_venergyEle_linear.size();

		assert(vs.size() == numLinear);

		for (size_t i = 0; i < numLinear; ++i)
		{
			(*this->m_venergyEle_linear[i]->GetMaterial())[Material::Property::StretchK] = vs[i];
		}

		this->DirtyDeformed();
	}

	void Model_MassSpring::SetBendingK(const VectorXd& vs)
	{
		size_t numHinges = this->m_venergyEle_hinge.size();

		assert(vs.size() == numHinges);

		for (size_t i = 0; i < numHinges; ++i)
		{
			(*this->m_venergyEle_hinge[i]->GetMaterial())[Material::Property::BendingK] = vs[i];
		}

		this->DirtyDeformed();
	}

	void Model_MassSpring::SetShearK(const VectorXd& vs)
	{
		size_t numCrosses = this->m_venergyEle_cross.size();

		assert(vs.size() == numCrosses);

		for (size_t i = 0; i < numCrosses; ++i)
		{
			(*this->m_venergyEle_cross[i]->GetMaterial())[Material::Property::ShearK] = vs[i];
		}

		this->DirtyDeformed();
	}
}