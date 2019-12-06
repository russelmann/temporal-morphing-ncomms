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

#include <PhySim/Models/Model_Reduction.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	/**
	* Model_Reduction_RBCloud
	*
	* TODO.
	*/
	class Model_Reduction_RBCloud: public Model_Reduction
	{
	public:

		enum Type
		{
			Quaternion,
			EulerAngles
		};

		/**
		* TODO.
		*/
		struct State : public Model::State
		{
		public:
			vector<Matrix3d>	m_vmR;

			State() {}
			virtual ~State() { }
			State(const State& s) { m_vmR = s.m_vmR; }
		};
		typedef shared_ptr<State> StateP;

	public:
		map<DoFSet*, DoFSet*>	 m_mmapDoF;

	protected:
		Type					 m_type;

		vector<RigidBody*>		 m_vBody;

		vector<Real*>			 m_vmapJac;
		int						 m_numRigidDoF;
		int						 m_numOtherDoF;

	public:

		virtual string GetName() const { return "RB Cloud Reduction"; }

		/**
		* Constructor.
		*/
		Model_Reduction_RBCloud();

		/**
		* Destructor.
		*/
		virtual ~Model_Reduction_RBCloud();

		// Initialization
		virtual void Init(Model* pModel, const vector<Polytope*>& vRB, Type type = Type::Quaternion);
		virtual void Free();

		virtual const int& NumRigidDoF() const { return this->m_numRigidDoF; }
		virtual const int& NumOtherDoF() const { return this->m_numOtherDoF; }

		IModel::StateP Model_Reduction_RBCloud::CreateState(Space s) const;
		bool HasState(IModel::StateP pS, Space s) const;
		void GetState(IModel::StateP pS, Space s) const;
		void SetState(const IModel::StateP pS, Space s);

		virtual vector<RigidBody*> GetRigidBodies() { return this->m_vBody; }

		virtual void UpdateInternalKinematics();
		virtual void ComputeAndStore_Jacobian();
		virtual void ComputeAndStore_Hessian(bool full = false);

	};
}