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

	//class Model_Reduction_Spectral : public Model_Reduction
	//{
	//public:

	//	// Construction
	//	Model_Reduction_Spectral();

	//	// Destruction
	//	virtual ~Model_Reduction_Spectral();

	//	// Initialization
	//	virtual void Init(Module* pModel, Real proportion);

	//	inline virtual string GetName() const { return "Spectral Reduction"; }

	//	virtual IModule::State* CreateState() const;
	//	virtual void GetState(IModule::State* pS, Space s = Space::DEF) const;
	//	virtual void SetState(const IModule::State* pS, Space s = Space::DEF);

	//	virtual void UpdateInternalKinematics();
	//	virtual void ComputeAndStore_Jacobian();

	//	virtual void Free();

	//};
}