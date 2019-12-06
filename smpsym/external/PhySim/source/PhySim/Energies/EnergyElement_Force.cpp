//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Energies/EnergyElement_Force.h>

#include <PhySim/Geometry/Node.h>

#include <PhySim/Models/Model_Particles.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	EnergyElement_Force::EnergyElement_Force(Model_Particles* pModel, const vector<DoFSet*>& vDoFs) : EnergyElement(pModel, NULL)
	{
		int countDoF = 0;
		this->m_vforces.resize(vDoFs.size());
		for (size_t i = 0; i < vDoFs.size(); ++i)
		{
			int numD = vDoFs[i]->GetNumDim();
			this->m_vforces[i].resize(numD);
			this->m_vforces[i].setZero();
			countDoF = countDoF + numD;
		}

		this->m_vDoFs = vDoFs;

		this->m_vgradient.resize(countDoF);
	}

	EnergyElement_Force::~EnergyElement_Force()
	{
		// Nothing to do here...
	}

	void EnergyElement_Force::Init()
	{
		// Nothing to do here...
	}

	void EnergyElement_Force::ComputeAndStore_Energy()
	{
		this->m_energy = 0;
		for (size_t i = 0; i < m_vDoFs.size(); ++i)
			this->m_energy += this->m_vDoFs[i]->GetPosition_x().dot(this->m_vforces[i]);
	}

	void EnergyElement_Force::ComputeAndStore_Gradient()
	{
		int dofOffset = 0;
		for (size_t i = 0; i < m_vDoFs.size(); ++i)
		{
			m_vgradient.block(dofOffset, 0, m_vDoFs[i]->GetNumDim(), 1) = this->m_vforces[i];
			dofOffset += m_vDoFs[i]->GetNumDim();
		}
	}
}