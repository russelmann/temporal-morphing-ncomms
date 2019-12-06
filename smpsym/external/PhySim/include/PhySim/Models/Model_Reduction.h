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
	* Model_Reduction.
	*
	* TODO.
	*/
	class Model_Reduction : public Model
	{
	protected:
		Model* m_pIntModel;

		CustomTimer m_timerComputeReduction;
		CustomTimer m_timerAssembleReduction;

		FastMatrixSd m_mJaco;

	public:

		/**
		* Constructor.
		*/
		Model_Reduction();

		/**
		* Destructor.
		*/
		virtual ~Model_Reduction();

		// Initialization
		virtual void Init(Model* pModel);

		virtual void Free() override;

		virtual Model* GetInternalModel() { return this->m_pIntModel; }

		virtual void SetFullDOFPosition(const VectorXd& vx) override;
		virtual void SetFreeDOFPosition(const VectorXd& vx) override;

		inline virtual const MatrixSd& GetJacobian();

		// Simulation preparation
		virtual void PrepareForSimulation() override;

		// Dirty flags
		virtual void DirtyDeformed() override;

		virtual bool IsDirty_Energy() const override;
		virtual bool IsDirty_Gradient(bool full = false) const override;
		virtual bool IsDirty_Hessian(bool full = false) const override;
		virtual bool IsDirty_Mass(bool full = false) const override;

		virtual bool IsDirty_Jacobian() const;

		// Mechanical computations
		virtual void ComputeAndStore_Energy() override;
		virtual void ComputeAndStore_Mass(bool full = false) override;
		virtual void ComputeAndStore_Gradient(bool full = false) override;
		virtual void ComputeAndStore_Hessian(bool full = false) override;

		virtual void UpdateInternalKinematics() = 0;
		virtual void ComputeAndStore_Jacobian() = 0;

	protected:
		virtual void InitDegreesOfFreedom() override { }

	};
}

