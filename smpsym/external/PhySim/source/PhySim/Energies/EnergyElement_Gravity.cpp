//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Energies/EnergyElement_Gravity.h>

#include <PhySim/Geometry/Node.h>
#include <PhySim/Models/Model_Particles.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	EnergyElement_Gravity::EnergyElement_Gravity(Model_Particles* pModel) : EnergyElement(pModel, NULL)
	{
		this->m_vDoFs = pModel->GetDoFSets();
		this->m_vgravity = Vector3d(0, -9.8, 0);
	}

	EnergyElement_Gravity::~EnergyElement_Gravity()
	{
		// Nothing to do here...
	}

	void EnergyElement_Gravity::Init()
	{
		// Nothing to do here...
	}

	void EnergyElement_Gravity::ComputeAndStore_Energy()
	{
		this->ComputeAndStore_Gradient();

		VectorXd vx;
		this->m_pModel->GetFullDOFPosition(vx);
		this->m_energy = vx.dot(m_vgradient);
	}

	void EnergyElement_Gravity::ComputeAndStore_Gradient()
	{
		this->m_vgradient.resize(m_pModel->GetNumFullDOF());
		int numNodes = (int) m_pModel->GetNodes().size();

		this->m_vgradient.setZero();

		for (size_t i = 0; i < numNodes; ++i)
		{
			m_vgradient.block(this->m_vDoFs[i]->GetOffset_Full(), 0, 3, 1) = m_vgravity;
		}

		m_vgradient = this->m_pModel->GetMass(true)*m_vgradient;
	}
}