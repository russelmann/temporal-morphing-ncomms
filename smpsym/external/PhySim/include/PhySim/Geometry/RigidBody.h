#pragma once

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

#include <PhySim/Geometry/Polytope.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class RigidBody
	{
	protected:

		Polytope* pPolytope;

		Matrix3d m_mR; // Constant accumulated rotation
		MatrixXd m_mT; // Local translation vectors t

		// Quaternion q such that R = R(q)*Rc
		// Centroid c such that x = c + R(e)*Rc*t
		vector<DoFSet*> m_vIntDoFs; 

		// The nodes of the rigidbody polytope
		vector<DoFSet*> m_vOutDoFs;

		// Node Jacobian w.r.t. DoF

		MatrixXd m_mJValues;
		MatrixXp m_mJPointers;

		// Coords Hessian w.r.t. DoF

		vector<MatrixXd> m_vmHValues;
		vector<MatrixXp> m_vmHPointers;

	public:

		RigidBody(Polytope* pPoly);

		virtual ~RigidBody();

		void Init();

		inline Matrix3d& R() { return m_mR; }
		inline MatrixXd& T() { return m_mT; }
		inline MatrixXd& J() { return m_mJValues; }

		inline Polytope* Poly() { return pPolytope; }
		inline DoFSet* CenterDoF() { return this->m_vIntDoFs[0]; };
		inline DoFSet* QuaterDoF() { return this->m_vIntDoFs[1]; };

		Vector3d GetCentroid();
		Matrix3d GetRotation();

		void RecomputeCentroid();
		void RecomputeRotation();

		inline int GetNumNodes() const { return (int)pPolytope->Nodes().size(); }

		void AssembleGlobal_Jacobian(VectorTd& vtotalJ);
		void AssembleGlobal_FastPreallocatedJacobian();
		void AllocateGlobal_Jacobian(const CoefMap& mp);

		void AssembleGlobal_VectorHessianProduct(const VectorXd& vv, VectorTd& vtotalJ);
		void AssembleGlobal_FastPreallocatedVectorHessianProduct(const VectorXd& vv);
		void AllocateGlobal_VectorHessianProduct(const CoefMap& mp);

		void UpdatePositions(Space s = Space::DEF);
		void UpdateJacobian(Space s = Space::DEF);
		void UpdateHessian(Space s = Space::DEF);

	};

}