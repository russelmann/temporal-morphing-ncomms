//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Energies/EnergyElement_External.h>

#include <PhySim/Geometry/Node.h>

#include <PhySim/Models/Model_Particles.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	EnergyElement_External::EnergyElement_External(Model_Particles* pModel, Material* pMaterial, const vector<DoFSet*>& vDoFs)
	{
		this->m_pModel = pModel;
		this->m_pMaterial = pMaterial;
		this->m_vDoFs = vDoFs;

		m_energy = 0;
		m_vgradient = VectorXd::Zero(0);
		m_vHessian.clear();
		m_vHessianP.clear();
	}

	EnergyElement_External::~EnergyElement_External()
	{
		// Nothing to do here...
	}

	void EnergyElement_External::AssembleGlobal_Gradient(VectorXd& vtotalGradient, bool full)
	{
		if (this->m_vgradient.size() == 0)
			return; // Ignore assembly

		if (!full)
		{
			assert(m_vgradient.size() == vtotalGradient.size());
			vtotalGradient = vtotalGradient + this->m_vgradient;
		}
		else
		{
         assert(false); ;
		}
	}

	void EnergyElement_External::AssembleGlobal_Hessian(VectorTd& vtotalHessian, bool full)
	{
		if (this->m_vHessian.size() == 0)
			return; // Ignore assembly

		if (!full)
		{
			// Extend the vector if necessary

			size_t numCoef = this->m_vHessian.size();
			if (vtotalHessian.capacity() < vtotalHessian.size() + numCoef)
				vtotalHessian.reserve(vtotalHessian.size() + numCoef);
	
			vtotalHessian.insert(vtotalHessian.end(), this->m_vHessian.begin(), this->m_vHessian.end());
		}
		else
		{
         assert(false); ;
		}
	}

	void EnergyElement_External::AllocateGlobal_Hessian(const CoefMap& mp, bool full)
	{
		if (this->m_vHessian.size() == 0)
			return; // Ignore assembly

		if (!full)
		{
			this->m_vHessianP.clear();

			size_t numCoef = m_vHessian.size();
			this->m_vHessianP.reserve(numCoef);
			for (int i = 0; i < numCoef; ++i)
			{
				this->m_vHessianP.push_back(
					Triplet<Real*>(
						this->m_vHessian[i].row(),
						this->m_vHessian[i].col(),
						mp.at(IntPair(this->m_vHessian[i].row(), this->m_vHessian[i].col())))
				);
			}
		}
		else
		{
         assert(false); ;
		}
	}

	void EnergyElement_External::AssembleGlobal_FastPreallocatedHessian(bool full)
	{
		if (this->m_vHessian.size() == 0)
			return; // Ignore assembly

		if (!full)
		{
			size_t numCoef = m_vHessian.size();
			assert(m_vHessianP.size() == numCoef);

			for (int i = 0; i < numCoef; ++i)
				(*this->m_vHessianP[i].value()) += this->m_vHessian[i].value();
		}
		else
		{
         assert(false); ;
		}
	}
}
