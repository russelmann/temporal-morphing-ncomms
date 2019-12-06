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

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	// Energies

	class EnergyElement;
	class EnergyElement_External;
	class EnergyElement_FaceNodeColl_LogB;
	class EnergyElement_FaceNodeColl_Quad;
	class EnergyElement_FaceNodeColl_Step;
	class EnergyElement_FaceNodeColl_Quad;
	class EnergyElement_FaceNodeColl_Step;
	class EnergyElement_FEM;
	class EnergyElement_FEM_MembraneCoNH;
	class EnergyElement_FEM_MembraneStVK;
	class EnergyElement_FEM_VolumeCoNH;
	class EnergyElement_FEM_VolumeStVK;
	class EnergyElement_Force;
	class EnergyElement_Gravity;
	class EnergyElement_PressureConstant;
	class EnergyElement_SpringCross;
	class EnergyElement_SpringHinge;
	class EnergyElement_SpringLinear;
	class EnergyElement_SpringLinear_UPoly;
	class EnergyElement_SpringLinear_MPoly;
	class EnergyElement_ShellHinge;
	class MassElement_Lumped;
	class Material;

	// Geometry

	class Cell;
	class Cell_Hexa;
	class Cell_Poly;
	class Cell_Tetra;
	class DoFSet;
	class Edge;
	class Face;
	class Face_Poly;
	class Face_Quad;
	class Face_Tri;
	class Node;
	class NodeEBD;
	class Polytope;
	class RigidBody;
	class TetGenReader;

	// Models

	class BC_Fixed;
	class BC_Force;
	class BC_Gravity;
	class BCondition;
	class Model;
	class Model_FEM;
	class Model_Inflatable;
	class Model_MassSpring;
	class Model_MassSpring_SMP;
	class Model_Particles;
	class Model_Reduction;
	class Model_Reduction_RBCloud;
	class Model_Reduction_Spectral;
	class Model_ThinShells;

	// Solvers

	class LinearSolver;
	class LinearSolver_CudaSC;
	class LinearSolver_EigenCG;
	class LinearSolver_EigenLDLT;
	class OptimProblem;
	class OptimProblem_BasicStatic;
	class OptimProblem_BasicDynamic;
	class OptimProblem_KnitroBridge;
	class OptimSolver_SQP;
	class OptimSolver_USQP_LS;
	class OptimSolver_KnitroBridge;

	// Interfaces
	
	class IDoFSet
	{
	public:

		virtual int GetNumDim() const = 0;

		virtual int GetId_Full() const = 0;
		virtual void SetId_Full(int id) = 0;

		virtual int GetId_Free() const = 0;
		virtual void SetId_Free(int id) = 0;

		virtual int GetOffset_Full() const = 0;
		virtual void SetOffset_Full(int o) = 0;

		virtual int GetOffset_Free() const = 0;
		virtual void SetOffset_Free(int o) = 0;

		virtual bool IsFixed() const = 0;

		virtual void Fix() = 0;
		virtual void Unfix() = 0;

		virtual const VectorXd& GetPosition(Space s) const = 0;
		virtual void SetPosition(const VectorXd& v, Space s) = 0;

		virtual const VectorXd& GetVelocity(Space s) const = 0;
		virtual void SetVelocity(const VectorXd& v, Space s) = 0;

	};

	class IPolytope
	{
	public:

		// Topology

		virtual int Order() = 0;
		virtual vector<Node*>& Nodes() = 0;
		virtual vector<Edge*>& Edges() = 0;
		virtual vector<Face*>& Faces() = 0;
		virtual vector<Edge*>& AdjacentEdges() = 0;
		virtual vector<Face*>& AdjacentFaces() = 0;
		virtual vector<Cell*>& AdjacentCells() = 0;

		// Geometry

		virtual void GetNodeMatrix(MatrixXd& mN, Space s = (Space) Space::DEF) const = 0;
		virtual void SetNodeMatrix(const MatrixXd& mN, Space s = (Space) Space::DEF) = 0;

		virtual Real ComputeVolume(Space s = Space::MAT) const = 0;
		virtual Vector3d ComputeCentroid(Space s = Space::MAT) const = 0;
		virtual Matrix3d ComputeRotation(Space f = Space::MAT, 
									     Space t = Space::DEF) const = 0;

		// Interpolation

		virtual bool IsValidParametric(const VectorXd& vp) const = 0;
		virtual NodeEBD ComputeEmbedding(const Vector3d& vp, Space s = Space::MAT) = 0;
		virtual Vector3d InterpolatePosition(const VectorXd& vp, Space s = Space::MAT) const = 0;
		virtual void ComputeShapeFunction(const VectorXd& vp, VectorXd& vN, Space s = Space::MAT) const = 0;
		virtual void ComputeShapeDerivative(const VectorXd& vp, MatrixXd& vN, Space s = Space::MAT) const = 0;
		virtual void ComputeNat2IsoTransform(VectorXd& vb, MatrixXd& mA, Space s = Space::MAT) const = 0;
		virtual void ComputeIso2NatTransform(VectorXd& vb, MatrixXd& mA, Space s = Space::MAT) const = 0;

		// Subelement

		virtual void SetSubelementPositions(Space s = Space::DEF) = 0;

		virtual DoFSet* DoF() = 0;
	};

	class IModel
	{
	public:
		struct State
		{
			State() {}
			virtual ~State() { }
			State(const State& s) { }
		};
		typedef shared_ptr<State> StateP;

	public:

		// Description
		virtual string GetName() const = 0;

		// Kinematic state
		virtual int GetNumFreeDOF() const = 0;
		virtual void GetFreeDOFPosition(VectorXd& vx) const = 0;
		virtual void SetFreeDOFPosition(const VectorXd& vx) = 0;
		virtual void GetFreeDOFVelocity(VectorXd& vv) const = 0;
		virtual void SetFreeDOFVelocity(const VectorXd& vv) = 0;

		virtual StateP CreateState(Space s = Space::DEF) const = 0;
		virtual bool HasState(StateP pS, Space s = Space::DEF) const = 0;
		virtual void GetState(StateP pS, Space s = Space::DEF) const = 0;
		virtual void SetState(const StateP pS, Space s = Space::DEF) = 0;

		// Mechanical state
		virtual const Real& GetEnergy() = 0;
		virtual const VectorXd& GetGradient(bool full = false) = 0;
		virtual const MatrixSd& GetHessian(bool full = false) = 0;

		// MassFull state
		virtual const MatrixSd& GetMass(bool full = false) = 0;

		// Boundary conditions
		virtual bool BoundaryConditionsLoaded() = 0;
		virtual void ResetBoundaryConditions() = 0;
		virtual void StepBoundaryConditions() = 0;
		virtual void FullBoundaryConditions() = 0;
		virtual void PresolveForFixedBoundary() = 0;

	};

	class IEnergyElement
	{
	public:

		virtual void Init() = 0;

		virtual void ComputeAndStore_Energy() = 0;
		virtual void ComputeAndStore_Gradient() = 0;
		virtual void ComputeAndStore_Hessian() = 0;

		virtual void AssembleGlobal_Gradient(VectorXd& vtotalVector, bool full = false) = 0;
		virtual void AssembleGlobal_Hessian(VectorTd& vtotalTriplets, bool full = false) = 0;

		virtual void AllocateGlobal_Hessian(const CoefMap& mp, bool full = false) = 0;
		virtual void AssembleGlobal_FastPreallocatedHessian(bool full = false) = 0;

		virtual const Real& GetIntegrationVolume() const = 0;
		virtual void SetIntegrationVolume(const Real& iv) = 0;

		virtual Real GetElementEnergy() const = 0;
		virtual const VectorXd& GetElementGradient() const = 0;
		virtual const MatrixXd& GetElementHessian() const = 0;
	};

	class IMassElement
	{
	public:

		virtual void Init() = 0;

		virtual void ComputeAndStore_Mass() = 0;

		virtual void AssembleGlobal_MassLumped(VectorXd& vMass, bool full = false) = 0;
		virtual void AssembleGlobal_MassMatrix(VectorTd& mMass, bool full = false) = 0;

		virtual Real GetElementMass() const = 0;
	};

	class IOptimProblem
	{
	public:

		virtual const string& GetName() const = 0;

		virtual int GetNumVariables() const = 0;
		virtual int GetNumConstraints() const = 0;

		virtual void GetVariables(VectorXd& vx) const = 0;
		virtual void SetVariables(const VectorXd& vx) = 0;

		virtual void GetUpperBound(VectorXd& vub) const = 0;
		virtual void GetLowerBound(VectorXd& vlb) const = 0;

		virtual bool GetEnergy(Real& e) = 0;
		virtual bool GetGradient(VectorXd& vg) = 0;
		virtual bool GetHessian(MatrixSd& mH) = 0;

		virtual bool GetConstraints(VectorXd& vc) = 0;
		virtual bool GetJacobian(MatrixSd& mJ) = 0;

		virtual bool IsFullyConstrained() = 0;

		// Callbacks

		virtual bool StartSolveCallback() = 0;
		virtual bool StartStepCallback() = 0;

		virtual bool EndSolveCallback() = 0;
		virtual bool EndStepCallback() = 0;

		virtual bool PreComputeStepCallback() = 0;
		virtual bool PosComputeStepCallback() = 0;
		
		virtual bool PrePerformStepCallback() = 0;
		virtual bool PosPerformStepCallback() = 0;

	};

	class IOptimSolver
	{
	public:

		virtual IOptimProblem& Problem() = 0;
		virtual SolveResult SolveStep() = 0;
		virtual SolveResult SolveFull() = 0;
		virtual void RegisterCallback_Step(void(*amazingCallback)(IOptimSolver*, void*), void* pReceiver) = 0;
		virtual void RegisterCallback_Full(void(*amazingCallback)(IOptimSolver*, void*), void* pReceiver) = 0;
		virtual bool& ForceSolveStop() = 0;

		virtual Real ComputeOptimality() = 0;
		virtual Real ComputeFeasibility() = 0;
		virtual bool IsOptimal() = 0;
		virtual bool IsFeasible() = 0;
	};

	class ILinearSolver
	{
	public:

		virtual void Init(const MatrixSd& mA, const LinearSolverOptions& options) = 0;
		virtual SolveResult Solve(MatrixSd& mA, const VectorXd& vb, VectorXd& vx) = 0;
		virtual SolveResult Solve(MatrixSd& mA, const MatrixXd& mB, MatrixXd& mX) = 0;
		virtual SolveResult Solve(MatrixSd& mA, const MatrixSd& mB, MatrixSd& mX) = 0;
	};

}

