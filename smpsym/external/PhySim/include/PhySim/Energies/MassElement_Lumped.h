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

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class MassElement_Lumped : public IMassElement
	{
	protected:

		Model_Particles *m_pModel;
		Polytope* m_pPolytope;
		Material* m_pMaterial;
		Real m_mass;

	public:

		MassElement_Lumped(Model_Particles* pModel, Polytope* pPolytope, Material* pMaterial);

		virtual ~MassElement_Lumped(void);

		virtual void Init();

		virtual void ComputeAndStore_Mass();

		//virtual void AddLumpedMassToNodes();

		virtual void AssembleGlobal_MassLumped(VectorXd& vMass, bool full = false);
		virtual void AssembleGlobal_MassMatrix(VectorTd& vMass, bool full = false);

		virtual Real GetElementMass() const { return this->m_mass; };

		inline virtual Material* GetMaterial() { return this->m_pMaterial; }
		inline virtual const Material* GetMaterial() const { return this->m_pMaterial; }
		inline virtual void SetMaterial(Material* pMat) { this->m_pMaterial = pMat; }

	};
}
