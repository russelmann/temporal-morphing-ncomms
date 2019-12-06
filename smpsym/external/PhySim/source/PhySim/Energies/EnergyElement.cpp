//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Energies/EnergyElement.h>

#include <PhySim/Geometry/Node.h>

#include <PhySim/Models/Model_Particles.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	EnergyElement::EnergyElement(Model_Particles* pModel, Material* pMaterial)
	{
		this->m_pModel = pModel;
		this->m_pMaterial = pMaterial;

		m_energy = 0;
		m_vgradient = VectorXd::Zero(0);
		m_mHessian = MatrixXd::Zero(0,0);
		m_mHessianPFree = MatrixXp::Zero(0,0);
		m_mHessianPFull = MatrixXp::Zero(0,0);

		this->m_vDoFs.clear();

		m_intVolume = 0;
	}

	EnergyElement::~EnergyElement()
	{
		// Nothing to do...
	}

	void EnergyElement::AssembleGlobal_Gradient(VectorXd& vtotalGradient, bool full)
	{
		if (this->m_vgradient.size() == 0)
			return; // Ignore assembly

		if (!full)
		{
			int localOffset = 0;

			size_t numDoFSet = GetStencilSize();
			for (size_t i = 0; i < numDoFSet; ++i)
			{
				int freeId = this->m_vDoFs[i]->GetId_Free();

				// Is simulated?
				if (freeId < 0)
				{
					localOffset += this->m_vDoFs[i]->GetNumDim();
					continue;
				}

				vtotalGradient.block(this->m_vDoFs[i]->GetOffset_Free(), 0, this->m_vDoFs[i]->GetNumDim(), 1) += m_vgradient.block(localOffset, 0, this->m_vDoFs[i]->GetNumDim(), 1);

				localOffset += this->m_vDoFs[i]->GetNumDim();
			}
		}
		else
		{
			int localOffset = 0;

			size_t numDoFSet = GetStencilSize();
			for (size_t i = 0; i < numDoFSet; ++i)
			{
				vtotalGradient.block(this->m_vDoFs[i]->GetOffset_Full(), 0, this->m_vDoFs[i]->GetNumDim(), 1) += m_vgradient.block(localOffset, 0, this->m_vDoFs[i]->GetNumDim(), 1);

				localOffset += this->m_vDoFs[i]->GetNumDim();
			}
		}
	}

	void EnergyElement::AssembleGlobal_Hessian(VectorTd& vtotalHessian, bool full)
	{
		if (this->m_mHessian.size() == 0)
			return; // Ignore assembly

		if (!full)
		{
			int localOffset_i = 0;
			int localOffset_j = 0;

			size_t numDoFSet = GetStencilSize();
			for (size_t i = 0; i < numDoFSet; ++i)
			{
				int freeId_i = this->m_vDoFs[i]->GetId_Free();
				if (freeId_i < 0)
				{
					localOffset_i += this->m_vDoFs[i]->GetNumDim();
					continue;
				}

				for (size_t j = 0; j < numDoFSet; ++j)
				{
					int freeId_j = this->m_vDoFs[j]->GetId_Free();
					if (freeId_j < 0)
					{
						localOffset_j += this->m_vDoFs[j]->GetNumDim();
						continue;
					}

					// Is upper triangular?
					if (freeId_i < freeId_j)
					{
						localOffset_j += this->m_vDoFs[j]->GetNumDim();
						continue;
					}

					assembleSparseMatrix(vtotalHessian,
						this->m_vDoFs[i]->GetOffset_Free(),
						this->m_vDoFs[j]->GetOffset_Free(),
						m_mHessian.block(localOffset_i, localOffset_j, this->m_vDoFs[i]->GetNumDim(), this->m_vDoFs[j]->GetNumDim()));

					localOffset_j += this->m_vDoFs[j]->GetNumDim();
				}

				localOffset_j = 0;
				localOffset_i += this->m_vDoFs[i]->GetNumDim();
			}
		}
		else
		{
			int localOffset_i = 0;
			int localOffset_j = 0;

			size_t numDoFSet = GetStencilSize();
			for (size_t i = 0; i < numDoFSet; ++i)
			{
				int fullId_i = this->m_vDoFs[i]->GetId_Full();

				for (size_t j = 0; j < numDoFSet; ++j)
				{
					int fullId_j = this->m_vDoFs[j]->GetId_Full();

					// Is upper triangular?
					if (fullId_i < fullId_j)
					{
						localOffset_j += this->m_vDoFs[j]->GetNumDim();
						continue;
					}

					assembleSparseMatrix(vtotalHessian,
						this->m_vDoFs[i]->GetOffset_Full(),
						this->m_vDoFs[j]->GetOffset_Full(),
						m_mHessian.block(localOffset_i, localOffset_j, this->m_vDoFs[i]->GetNumDim(), this->m_vDoFs[j]->GetNumDim()));

					localOffset_j += this->m_vDoFs[j]->GetNumDim();
				}

				localOffset_j = 0;
				localOffset_i += this->m_vDoFs[i]->GetNumDim();
			}
		}
	}

	void EnergyElement::AllocateGlobal_Hessian(const CoefMap& mp, bool full)
	{
		if (this->m_mHessian.size() == 0)
			return; // Ignore assembly


		if (!full)
		{
			int localOffset_i = 0;
			int localOffset_j = 0;

			size_t numDoFSet = GetStencilSize();

			this->m_mHessianPFree.setZero();

			for (size_t i = 0; i < numDoFSet; ++i)
			{
				int freeId_i = this->m_vDoFs[i]->GetId_Free();
				if (freeId_i < 0)
				{
					localOffset_i += this->m_vDoFs[i]->GetNumDim();
					continue;
				}

				for (size_t j = 0; j < numDoFSet; ++j)
				{
					int freeId_j = this->m_vDoFs[j]->GetId_Free();
					if (freeId_j < 0)
					{
						localOffset_j += this->m_vDoFs[j]->GetNumDim();
						continue;
					}

					// Is upper triangular?
					if (freeId_i < freeId_j)
					{
						localOffset_j += this->m_vDoFs[j]->GetNumDim();
						continue;
					}

					for (int ii = 0; ii < this->m_vDoFs[i]->GetNumDim(); ++ii)
					{
						for (int jj = 0; jj < this->m_vDoFs[j]->GetNumDim(); ++jj)
						{
							m_mHessianPFree(localOffset_i + ii, localOffset_j + jj) = mp.at(IntPair(this->m_vDoFs[i]->GetOffset_Free() + ii, this->m_vDoFs[j]->GetOffset_Free() + jj));
						}
					}

					localOffset_j += this->m_vDoFs[j]->GetNumDim();
				}

				localOffset_j = 0;
				localOffset_i += this->m_vDoFs[i]->GetNumDim();
			}
		}
		else
		{
			int localOffset_i = 0;
			int localOffset_j = 0;

			size_t numDoFSet = GetStencilSize();

			this->m_mHessianPFull.setZero();

			for (size_t i = 0; i < numDoFSet; ++i)
			{
				int fullId_i = this->m_vDoFs[i]->GetId_Full();

				for (size_t j = 0; j < numDoFSet; ++j)
				{
					int fullId_j = this->m_vDoFs[j]->GetId_Full();

					// Is upper triangular?
					if (fullId_i < fullId_j)
					{
						localOffset_j += this->m_vDoFs[j]->GetNumDim(); 
						continue;
					}

					for (int ii = 0; ii < this->m_vDoFs[i]->GetNumDim(); ++ii)
					{
						for (int jj = 0; jj < this->m_vDoFs[j]->GetNumDim(); ++jj)
						{
							m_mHessianPFull(localOffset_i + ii, localOffset_j + jj) = mp.at(IntPair(this->m_vDoFs[i]->GetOffset_Full() + ii, this->m_vDoFs[j]->GetOffset_Full() + jj));
						}
					}

					localOffset_j += this->m_vDoFs[j]->GetNumDim();
				}

				localOffset_j = 0;
				localOffset_i += this->m_vDoFs[i]->GetNumDim();
			}
		}
	}

	void EnergyElement::AssembleGlobal_FastPreallocatedHessian(bool full)
	{
		if (this->m_mHessian.size() == 0)
			return; // Ignore assembly

		if (!full)
		{
			size_t numDoFSet = GetStencilSize();

			int localOffset_i = 0;
			int localOffset_j = 0;

			for (size_t i = 0; i < numDoFSet; ++i)
			{
				int freeId_i = this->m_vDoFs[i]->GetId_Free();
				if (freeId_i < 0)
				{
					localOffset_i += this->m_vDoFs[i]->GetNumDim();
					continue;
				}

				for (size_t j = 0; j < numDoFSet; ++j)
				{
					int freeId_j = this->m_vDoFs[j]->GetId_Free();
					if (freeId_j < 0)
					{
						localOffset_j += this->m_vDoFs[j]->GetNumDim();
						continue;
					}

					for (int ii = 0; ii < this->m_vDoFs[i]->GetNumDim(); ++ii)
					{
						for (int jj = 0; jj < this->m_vDoFs[j]->GetNumDim(); ++jj)
						{
							if (m_mHessianPFree(localOffset_i + ii, localOffset_j + jj) != NULL)
							{
								*m_mHessianPFree(localOffset_i + ii, localOffset_j + jj) += m_mHessian(localOffset_i + ii, localOffset_j + jj);
							}
						}
					}

					localOffset_j += this->m_vDoFs[j]->GetNumDim();
				}

				localOffset_j = 0;
				localOffset_i += this->m_vDoFs[i]->GetNumDim();
			}
		}
		else
		{
			size_t numDoFSet = GetStencilSize();

			int localOffset_i = 0;
			int localOffset_j = 0;

			for (size_t i = 0; i < numDoFSet; ++i)
			{
				int fullId_i = this->m_vDoFs[i]->GetId_Full();

				for (size_t j = 0; j < numDoFSet; ++j)
				{
					int fullId_j = this->m_vDoFs[j]->GetId_Full();

					for (int ii = 0; ii < this->m_vDoFs[i]->GetNumDim(); ++ii)
					{
						for (int jj = 0; jj < this->m_vDoFs[j]->GetNumDim(); ++jj)
						{
							if (m_mHessianPFull(localOffset_i + ii, localOffset_j + jj) != NULL)
							{
								*m_mHessianPFull(localOffset_i + ii, localOffset_j + jj) += m_mHessian(localOffset_i + ii, localOffset_j + jj);
							}
						}
					}

					localOffset_j += this->m_vDoFs[j]->GetNumDim();
				}

				localOffset_j = 0;
				localOffset_i += this->m_vDoFs[i]->GetNumDim();
			}
		}
	}
}