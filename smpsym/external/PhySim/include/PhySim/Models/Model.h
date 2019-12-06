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

#include <PhySim/Utils/CustomTimer.h>

#include <PhySim/Energies/Material.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	/**
	* Module (Abstract).
	*
	* Base abstract class that can be inherited from for model implementation.
	* It provides the basic functionality regarding energy, gradient and HessianFull
	* management, DoF get/set, state serialization, global derivative testing and
	* profiling. Any class deriving from Module must implement the following methods:
	*
	* - InitDegreesOfFreedom: Vector m_vDoF and m_numFullDoF are initialized
	* - ComputeAndStore_Mass: MassFull matrix is computed and stored in m_mMassFree
	* - ComputeAndStore_Energy: Energy is computed and stored in m_energy
	* - ComputeAndStore_Gradient: GradientFull is computed and stored in m_vgradFree
	* - ComputeAndStore_Hessian: HessianFull is computed and stored in m_mHessFree
	*/
	class Model : public IModel
	{

	public:

		/**
		* Basic structure for safe state saving and loading. This just constains the current value of the DoF.
		* Any class deriving from Module that stores additional data to determine the state of the simulation
		* should also extend this structure as well as CreateState(), GetState() and SetState() methods. 
		* 
		* For example, when creating a physics model that uses EulerAngles angles e to represent rotations, one might 
		* want to avoid DoF to store high values (due to Gimbal lock). A possible method to workaround this issue
		* is to occasionally cache the current rotation in some other format (e.g. a rotation matrix, Rc) and 
		* compute the true transformation as a composition of both rotations R = R(e)*Rc. Such a cached matrix
		* is not a DoF but certenly affects the state of the simulation.
		*
		* Safely saving and loading the kinematic state of the simulation is essential to consistently track
		* changes in the energy and its derivatives. For instance, for estimating derivatives using finite
		* differences or to assess whether or not back-tracking line-search is improving the solution when
		* minimizing the potential energy.
		*/
		struct State : public IModel::State
		{
		public:
			VectorXd m_vx;

			State() {}
			virtual ~State() { }
			State(const State& s) { m_vx = s.m_vx; }
		};
		typedef shared_ptr<State> StateP;

	protected:

		vector<DoFSet*> m_vDoFs;

		int m_numFullDoF;
		int m_numFreeDoF;

		Real m_energy;
		VectorXd m_vgradFree;
		VectorXd m_vgradFull;
		FastMatrixSd m_mHessFree;
		FastMatrixSd m_mMassFree;
		FastMatrixSd m_mHessFull;
		FastMatrixSd m_mMassFull;
		bool m_fastAssembly;

		DirtyFlags m_dirtyFlags;

		bool m_isProfiling;
		CustomTimer m_timerComputeEnergy;
		CustomTimer m_timerComputeGradient;
		CustomTimer m_timerComputeHessian;
		CustomTimer m_timerAssembleEnergy;
		CustomTimer m_timerAssembleGradient;
		CustomTimer m_timerAssembleHessian;

		vector<Material> m_vmat_mass;
		vector<IMassElement*> m_vmassEle;
		vector<IEnergyElement*> m_venergyEle;
		EnergyElement_Gravity* m_pEle_gravity;

		vector<EnergyElement_External*> m_vexternEle;

		vector<BCondition*>		m_vBC;
		bool					m_presolveFixedBC;
		MatrixSd				m_mFixedPerm;
		int						m_numFixedDoF;
		int						m_numUnfixDoF;

		LinearSolver*			m_pPreBCSolver;

		bool					m_gradientFD;
		bool					m_hessianFD;

	public:

		/**
		* Returns the name of the simulation model for logging.
		*/
		virtual string GetName() const override { return "Base Model"; }

		/**
		* Constructor.
		*/
		Model();

		/**
		* Destructor.
		*/
		virtual ~Model();

		// Initialization

		/**
		* Initializes the model. This method first calls Free() and then
		* calls InitDegreesOfFreedom() to initialize the rest of the data.
		*/
		virtual void Init();

		/**
		* Frees the model data. This method frees the memory allocated to store 
		* the DoF of the simulation model. Any class deriving from  Module that 
		* allocates additional memory must override this method as it is called 
		* by Module destructor. Make sure the memory allocated by Module is also
		* freed by always calling parent Free() method.
		*/
		virtual void Free();

		// Profiling

		/**
		* Get if performance profiling is active.
		*/
		virtual bool GetIsProfiling() const { return this->m_isProfiling; }

		/**
		* Set if performance profiling is active.
		*/
		virtual void SetIsProfiling(bool ip) { this->m_isProfiling = ip; }

		virtual void SetPresolveFixed(bool set) { this->m_presolveFixedBC = set; }
		virtual bool GetPresolveFixed() const { return this->m_presolveFixedBC; }

		// State save/load

		virtual IModel::StateP CreateState(Space s = Space::DEF) const override;
		virtual bool HasState(IModel::StateP pS, Space s = Space::DEF) const override;
		virtual void GetState(IModel::StateP pS, Space s = Space::DEF) const override;
		virtual void SetState(const IModel::StateP pS, Space s = Space::DEF) override;

		// Kinematic DOF 

		/**
		* Returns a reference to the current DoF vector.
		*/
		inline virtual const vector<DoFSet*>& GetDoFSets() const;

		inline virtual int GetNumFullDOF() const;
		virtual void GetFullDOFPosition(VectorXd& vx) const;
		virtual void SetFullDOFPosition(const VectorXd& vx);

		inline virtual int GetNumFreeDOF() const override;
		virtual void GetFreeDOFPosition(VectorXd& vx) const override;
		virtual void SetFreeDOFPosition(const VectorXd& vx) override;
		virtual void GetFreeDOFVelocity(VectorXd& vv) const override;
		virtual void SetFreeDOFVelocity(const VectorXd& vv) override;

		// Mechanical state
		inline virtual const Real& GetEnergy() override;
		inline virtual const VectorXd& GetGradient(bool full = false) override;
		inline virtual const MatrixSd& GetHessian(bool full = false) override;
		inline virtual const MatrixSd& GetMass(bool full = false) override;

		// Mechanical computations
		virtual void ComputeAndStore_Energy() = 0;
		virtual void ComputeAndStore_Mass(bool full = false) = 0;
		virtual void ComputeAndStore_Gradient(bool full = false) = 0;
		virtual void ComputeAndStore_Hessian(bool full = false) = 0;

		virtual void ComputeAndStore_Gradient_FD(bool full = false);
		virtual void ComputeAndStore_Hessian_FD(bool full = false);

		// Simulation preparation
		virtual void PrepareForSimulation();

		// Dirty flags
		virtual void DirtyUndeformed();
		virtual void DirtyDeformed();
		virtual bool IsDirty_Energy() const;
		virtual bool IsDirty_Gradient(bool full = false) const;
		virtual bool IsDirty_Hessian(bool full = false) const;
		virtual bool IsDirty_Mass(bool full = false) const;

		// Finite difference testing
		virtual void TestGlobalGradient();
		virtual void TestGlobalHessian();

		virtual const vector<IMassElement*>& GetMassElements() const { return this->m_vmassEle; }
		virtual const vector<IEnergyElement*>& GetEnergyElements() const { return this->m_venergyEle; }

		virtual const vector<EnergyElement_External*> GetExternalEnergyComponents() const { return this->m_vexternEle; }
		virtual void AddExternalEnergyComponent(EnergyElement_External* pEle);
		virtual void DelExternalEnergyComponent(EnergyElement_External* pEle);

		virtual bool BoundaryConditionsLoaded() override;
		virtual void ResetBoundaryConditions() override;
		virtual void StepBoundaryConditions() override;
		virtual void FullBoundaryConditions() override;
		virtual void PresolveForFixedBoundary() override;

		virtual void AddBoundaryCondition(BCondition* pBC);
		virtual void RemoveBoundaryCondition(BCondition* pBC);
		virtual void ClearBoundaryConditions();
		virtual const vector<BCondition*>& GetBoundaryConditions() { return this->m_vBC; }

		virtual bool& FiniteDifferent_Gradient() { return this->m_gradientFD; }
		virtual bool& FiniteDifferent_Hessian() { return this->m_hessianFD; }

		virtual vector<DoFSet*> SelectDoF(const Vector3d& vboxMin, const Vector3d& vboxMax, Space s = Space::DEF) const { return vector<DoFSet*>(); }

	protected:
		virtual void InitDegreesOfFreedom() = 0;

	};
}



