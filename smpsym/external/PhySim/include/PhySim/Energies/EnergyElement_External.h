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

#include <PhySim/Energies/Material.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class Node;

	class Model_Particles;

	class EnergyElement_External : public IEnergyElement
	{
	protected:

		Model_Particles *m_pModel;
		vector<DoFSet*> m_vDoFs;
		Material* m_pMaterial;

		double m_energy;
		VectorXd m_vgradient;
		VectorTd m_vHessian;
		VectorTp m_vHessianP;

	public:

		EnergyElement_External(Model_Particles* pModel, Material* pMaterial, const vector<DoFSet*>& vDoFs);

		virtual ~EnergyElement_External(void);

		inline virtual size_t GetStencilSize() { return this->m_vDoFs.size(); }
		inline virtual const vector<DoFSet*>& GetDoFStencil() const { return this->m_vDoFs; }
		inline virtual void SetDoFStencil(const vector<DoFSet*> vns) { this->m_vDoFs = vns; }

		/**
		* This will be called when initializing the elements to the undeformed state.
		*/
		virtual void Init() override {}

		/**
		* Compute the external component energy and store it in the variable m_energy.
		*/
		virtual void ComputeAndStore_Energy() override { }

		/**
		* Compute the external gradient and store it in the variable m_vgradient. The 
		* used here should be GLOBAL, i.e. the vector m_vgradient must have size N, where
		* N = m_pM->GetNumFullDOF(). The subvector of the gradient corresponding to the 
		* i-th DoF has size m_vDoF[i]->GetNumDim() and is assembled from vector position k,
		* with k = m_vDoF[i]->GetOffset_Free().
		*
		* NOTE 1:		Some of the DoF might not be used if they are fixed. In such case
		*				the methods m_vDoF[i]->GetId_Free() and m_vDoF[i]->GetOffset_Free()
		*				returns -1. Do not assemble any subvector including a fixed DoF.
		*/
		virtual void ComputeAndStore_Gradient() override { }

		/**
		* Compute the external HessianFull and store it in the variable m_vHessian,
		* in the form of vector of triplets (VectorTd). The vector m_vHessian must have size
		* equal to the non-zero coefficients of the external HessianFull and the indexing used 
		* here should be GLOBAL. The submatrix of the HessianFull corresponding to the i-th and
		* j-th DoFs has size [m_vDoF[i]->GetNumDim() x m_vDoF[j]->GetNumDim()] and it-s 
		* assembled from index (k,l), where:
		*
		* k = m_vDoF[i]->GetOffset_Free().
		* l = m_vDoF[j]->GetOffset_Free().
		*
		* NOTE 1:		Some of the DoF might not be used if they are fixed. In such case
		*				the methods m_vDoF[i]->GetId_Free() and m_vDoF[i]->GetOffset_Free()
		*				return -1. Do not assemble any submatrix including a fixed DoF.
		*
		* NOTE 2:		Only the lower triangular part of the HessianFull is actually used.
		*				The rest of the matrix will be ignored. Consider not assembling
		*				any submatrix (i,j) pertaining to the upper triangular part, i.e.
		*				m_vDoF[i]->GetId_Free() < m_vDoF[j]->GetId_Free().
		*/
		virtual void ComputeAndStore_Hessian() override { }

		virtual void AssembleGlobal_Gradient(VectorXd& vtotalGradient, bool full = false) override;
		virtual void AssembleGlobal_Hessian(VectorTd& vtotalHessian, bool full = false) override;

		virtual void AllocateGlobal_Hessian(const CoefMap& mp, bool full = false) override;
		virtual void AssembleGlobal_FastPreallocatedHessian(bool full = false) override;

		inline Real GetElementEnergy() const override { return this->m_energy; }

		/**
		* TODO: Don't use it yet
		*/
		inline const VectorXd& GetElementGradient() const override { return VectorXd(); }

		/**
		* TODO: Don't use it yet
		*/
		inline const MatrixXd& GetElementHessian() const override { return MatrixXd(); }

		inline virtual Material* GetMaterial() { return this->m_pMaterial; }
		inline virtual const Material* GetMaterial() const { return this->m_pMaterial; }
		inline virtual void SetMaterial(Material* pMat) { this->m_pMaterial = pMat; }

	};
}
