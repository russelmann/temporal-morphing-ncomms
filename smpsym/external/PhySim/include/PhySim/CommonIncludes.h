//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#pragma once

// Standard

#include <memory>
#include <fstream>
#include <sstream>
#include <iostream>
#include <set>

#include <stdio.h>
//#include <direct.h>

// Eigen

#include <Eigen/Dense>
#include <Eigen/Sparse>

#define GetCurrentDir _getcwd

//#ifdef WINDOWS
//#include <direct.h>
//#define GetCurrentDir _getcwd
//#else
//#include <unistd.h>
//#define GetCurrentDir getcwd
//#endif

// Typedefs

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	// Standard types
	typedef double Real;
	typedef vector<int> iVector;
	typedef vector<float> fVector;
	typedef vector<double> dVector;
	typedef vector<bool> bVector;
	typedef pair<int, int> IntPair;
	typedef map<IntPair, Real*> CoefMap;

	// Eigen types
	typedef vector<Triplet<Real> > VectorTd;
	typedef vector<Triplet<Real*> > VectorTp;
	typedef SparseMatrix<Real> MatrixSd;

	typedef Matrix<Real*, Dynamic, Dynamic> MatrixXp;

	// Other types
	typedef struct Frame3d
	{
		Frame3d()
		{
			this->tan.setZero();
			this->nor.setZero();
			this->bin.setZero();
		}

		~Frame3d()
		{
			// Nothing here...
		}

		Frame3d(const Frame3d& toCopy)
		{
			this->tan = toCopy.tan;
			this->nor = toCopy.nor;
			this->bin = toCopy.bin;
		}

		bool operator == (const Frame3d& other) const
		{
			return 
				(tan - other.tan).squaredNorm() < 1e-9 &&
				(nor - other.nor).squaredNorm() < 1e-9 &&
				(bin - other.bin).squaredNorm() < 1e-9;
		}

		Frame3d reverse() const
		{
			Frame3d out = *this;
			out.tan = -out.tan;
			out.bin = -out.bin;
			return out;
		}

		Vector3d tan;
		Vector3d nor;
		Vector3d bin;
	}
	Frame3d;

	enum MaterialModel
	{
		CoRot,	// Co-rotational
		StVK,	// Saint Venant-Kirchhoff
		CoNH,	// Compressible Neo-Hookean
		CoMR,	// Compressible Mooney-Rivlin
		InNH,	// Incompressible Neo-Hookean
		InMR,	// Incompressible Mooney-Rivlin
		Ogden	// Ogden
	};

	enum Discretization
	{
		Nodes,
		Edges,
		Triangles,			// Linear triangle
		Quadrangles,		// Bilinear quadrangle
		Tetrahedra4,		// Linear tetrahedron
		Tetrahedra10,		// Quadratic tetrahedron
		Hexahedra			// Trilinear hexahedron
	};

	enum LSSolverType
	{
		LS_EigenCG,
		LS_EigenLDLT,
		LS_EigenLU,
		LS_BiCGSTAB,
		LS_CholmodLDLT,
		LS_SSparseSPQR,
		LS_CUDALU,
		LS_CUDAQR,
		LS_CUDASC
	};

	enum QPSolverType
	{
		QP_Steepest,
		QP_Newton,
		QP_Gauss,
		QP_BFGS_D,
		QP_BFGS_I,
		QP_LBFGS
	};

	enum StepSelType
	{
		SS_LineSearch,
		SS_TrustRegion
	};

	enum SolveResult
	{
		Success,
		Failure,	// Generic failure
		Stopped,	// Solver stop forced by user
		MaxIter,	// Maximum iterations reached
		NonDesc,	// Non-descendent step found
		Singular,	// Singular matrix (rank deficient)
		NonSPD		// Non symmetric positive definite
	};

	enum DirtyFlags
	{
		None = 0,
		Energy = 1,
		GradientFull = 2,
		HessianFull = 4,
		MassFull = 8,
		GradientFree = 16,
		HessianFree = 32,
		MassFree = 64,
		Reduction = 128
	};

	typedef struct OptimSolverOptions
	{
		Real maxError;
		int maxIters;
		LSSolverType lsSolverType;
		QPSolverType qpSolverType;
		StepSelType stepSelType;
		int lineSearchIters;
		Real lineSearchFactor;
		Real maxStepSize;
		bool profileTime;

		OptimSolverOptions()
		{
			maxError = 1e-9;
			maxIters = 100;
			stepSelType = StepSelType::SS_LineSearch;
			lsSolverType = LSSolverType::LS_EigenLDLT;
			qpSolverType = QPSolverType::QP_Newton;
			lineSearchIters = 10;
			lineSearchFactor = 0.5;
			maxStepSize = 1.0;
			profileTime = true;
		}
	}
	OptimSolverOptions;

	typedef struct LinearSolverOptions
	{
		Real maxError;
		int maxIters;
		int regIters;
		Real regSign;
		LSSolverType type;
		bool profileTime;

		LinearSolverOptions()
		{
			type = LSSolverType::LS_EigenLDLT;
			maxError = 1e-3;
			maxIters = 1000;
			regIters = 3;
			regSign = 1;
			profileTime = true;
		}
	}
	LinearSolverOptions;

	inline DirtyFlags operator|(DirtyFlags a, DirtyFlags b)
	{
		return static_cast<DirtyFlags>(static_cast<int>(a) | static_cast<int>(b));
	}

	inline DirtyFlags operator&(DirtyFlags a, DirtyFlags b)
	{
		return static_cast<DirtyFlags>(static_cast<int>(a) & static_cast<int>(b));
	}

	inline DirtyFlags operator~(DirtyFlags a)
	{
		return static_cast<DirtyFlags>(~static_cast<int>(a));
	}

	static const int NUMSPACES = 2;
	enum Space
	{
		MAT = 0,
		DEF = 1,
	};

}

// Utils

#include <PhySim/Utils/Utils.h>
