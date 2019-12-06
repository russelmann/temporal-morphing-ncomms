//=============================================================================
//
//   Co-Rotational Finite Element Method source code
//   Based on Physics-Based Animation book, chapter 10
//
//   Based on original implementation by
//							Miguel A. Otaduy,    URJC Madrid
//							Alvaro G. Perez,     URJC Madrid
//							and Javier S. Zurdo, URJC Madrid
//
//	Authors:
//			Krisztian Amit Birkas, McGill Uni.
//			Jesus Perez Rodriguez, IST Austria
//
//=============================================================================

#pragma once

#include <PhySim/CommonIncludes.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	struct LT_Triplet
	{
		bool operator()(const Triplet<double> &a, const Triplet<double> &b)
		{
			return (a.row() < b.row()) || (a.row() == b.row() && a.col() < b.col());
		}
	};

	// 3D MESHES

	bool readTriMesh_Obj(const string& path, MatrixXd& mV, MatrixXi& mF);
	bool writeTriMesh_Obj(const string& path, const MatrixXd& mV, const MatrixXi& mF);

	// FILE SYSTEM

	bool listDirectoryFiles(const string& dirPath, vector<string>& vfiles);

	// STRINGS

	void splitString(const string& input, const string& delim, vector<string>& vs);

	// LOG/TRACE

	bool readFileLines(const string& path, vector<string>& lines);
	string vectorToString_CSV(const VectorXd& v);
	string matrixToString_CSV(const MatrixXd& M);
	void logFile(const string& path, const string& content, bool append = false);
	void logSimu(const char *format, ...);
	void logTime(const char *format, ...);
	void logSimuStream(ostream* ostream);
	void logTimeStream(ostream* ostream);
	bool readMatrix_CSV(const string& path, MatrixXd& mV);
	bool writeMatrix_CSV(const string& path, MatrixXd& mV);

	// EIGEN/STL

	void add3D(int i, const Vector3d& v, dVector& vd);
	void set3D(int i, const Vector3d& v, dVector& vd);
	Eigen::Vector3d get3D(int i, const dVector& vd);

	void getSubvector(int offset, int numEle, const dVector& vin, dVector& vout);
	void setSubvector(int offset, int numEle, const dVector& vin, dVector& vout);

	void getSubvector(int offset, int numEle, const VectorXd& vin, VectorXd& vout);
	void setSubvector(int offset, int numEle, const VectorXd& vin, VectorXd& vout);

	dVector toSTL(const Eigen::VectorXd& v);
	iVector toSTL(const Eigen::VectorXi& v);
	dVector toSTL(const Eigen::Vector3d& v);
	Eigen::VectorXd toEigen(const dVector& v);
	Eigen::VectorXi toEigen(const iVector& v);

	void assembleSparseMatrix(VectorTd& triplets, unsigned int row, unsigned int col, const Matrix3d& m33);

	void eigenTripletsToSparseMatrix(const VectorTd& triplets, MatrixSd& M);
	void eigenSparseMatrixToPointers(const MatrixSd& M, VectorTp& triplets);
	void eigenSparseMatrixToTriplets(const MatrixSd& M, VectorTd& triplets);

	void getDiagonalIndices(const VectorTd& vT, int rows, int cols, iVector& vDIdx);

	void transformTripletValuesToCUDA(const VectorTd& triplets, vector<double> & cuda);
	void transformCUDAToTripletValues(const vector<double>& cuda, VectorTd& triplets);
	void transformTripletValuesToCUDA(const VectorTd& triplets, vector<float>& cuda);
	void transformCUDAToTripletValues(const vector<float>& cuda, VectorTd& triplets);
	void transformEigenVectorToCUDA(const VectorXd& veigen, vector<double>& cuda);
	void transformEigenVectorToCUDA(const VectorXd& veigen, vector<float>& cuda);
	void transformCUDAToEigenVector(const vector<float>& cuda, VectorXd& veigen);
	void transformCUDAToEigenVector(const vector<double>& cuda, VectorXd& veigen);
	void transformEigenTripletsToCUDA(const VectorTd& triplets, vector<int>& rows, vector<int>& cols, vector<double>& vals);
	void transformEigenTripletsToCUDA(const VectorTd& triplets, vector<int>& rows, vector<int>& cols, vector<float>& vals);
	void transformCUDAToEigenTriplets(const vector<int>& rows, const vector<int>& cols, const vector<float>& vals, VectorTd& triplets);
	void transformCUDAToEigenTriplets(const vector<int>& rows, const vector<int>& cols, const vector<double>& vals, VectorTd& triplets);


	// GEOMETRY

	void computeBoundingBox(MatrixXd& mV, Vector3d& vmin, Vector3d& vmax);

	Frame3d rotateFrame(const Frame3d& F, Real angle);
	Frame3d rotateFrame(const Frame3d& F, Matrix3d& R);

	void rotateFrame(const vector<Frame3d>& vF, const dVector& va, vector<Frame3d>& vFRot);

	Frame3d parallelTransport(const Frame3d& F, const Vector3d& t2);
	Vector3d parallelTransport(const Vector3d& u, const Vector3d& t1, const Vector3d& t2);
	Vector3d parallelTransportNormalized(const Vector3d& u, const Vector3d& t1, const Vector3d& t2);

	Vector3d curvatureBinormal(const Vector3d& e0, const Vector3d& e1);

	void orthonormalizeFrame(Frame3d& F);

	bool lineLineIntersect(const dVector& p1, const dVector& p2, const dVector& p3, const dVector& p4, double& mua, double& mub, double tol);

	Matrix3d transformEulerToMatrix(const Vector3d& ve);
	Vector3d transformMatrixToEuler(const Matrix3d& mR);

	// INTERPOLATIONS

	void hermiteInterpolation(const vector<Vector3d>& vin, vector<Vector3d>& vout, bool loop, int outnPoint, const Vector3d& tanIni = Vector3d(), const Vector3d& tanEnd = Vector3d());
	void hermiteInterpolation(const vector<Vector3d>& vin, vector<Vector3d>& vout, bool loop, double outTol, const Vector3d& tanIni = Vector3d(), const Vector3d& tanEnd = Vector3d());
	void hermiteInterpolation(const vector<Frame3d>& vin, const dVector& vsin, vector<Frame3d>& vout, bool loop, int outnp);
	void hermiteInterpolation(const vector<Vector2d>& vin, const dVector& vsin, vector<Vector2d>& vout, bool loop, int outnp);
	void hermiteInterpolation(const vector<Frame3d>& vin, const dVector& vsin, vector<Frame3d>& vout, bool loop, double outTol);
	void hermiteInterpolation(const vector<Vector2d>& vin, const dVector& vsin, vector<Vector2d>& vout, bool loop, double outTol);

	// OTHER MATH UTILS

	Real computePolarDecomposition(const Matrix3d& M, Matrix3d& Q, Matrix3d& S, Real tol);

	void computeBestRigidTransform(const MatrixXd& mCloud0, const MatrixXd& mCloud1, MatrixXd& mR, VectorXd& vt);

	void eulerRotation(const double *eulerAngles, const double *vectorOriginal, double vectorRotated[3]);

	// TYPES

	struct FastMatrixSd
	{
	public:
		CoefMap m_coeffMap;
		MatrixSd m_msparseMatrix;
		VectorTd m_vvalueTriplets;
		VectorTp m_vpointTriplets;

		FastMatrixSd(int _numRows = 0, int _numCols = 0)
		{
			m_coeffMap.clear();
			m_vvalueTriplets.clear();
			m_msparseMatrix.resize(_numRows, _numCols);
		}

		void BuildMatrixFromTriplets()
		{
			this->m_msparseMatrix.setFromTriplets(
				this->m_vvalueTriplets.begin(),
				this->m_vvalueTriplets.end());
			this->m_msparseMatrix.makeCompressed();
		}

		void BuildTripletsFromMatrix()
		{
			eigenSparseMatrixToTriplets(this->m_msparseMatrix, this->m_vvalueTriplets);
		}

		void BuildMappingFromMatrix()
		{
			eigenSparseMatrixToPointers(this->m_msparseMatrix, this->m_vpointTriplets);

			this->m_coeffMap.clear();

			size_t numCoeff = m_vpointTriplets.size();

			for (size_t i = 0; i < numCoeff; ++i)
			{
				const Triplet<Real*>& pointer = m_vpointTriplets[i];
				IntPair key = IntPair(pointer.row(), pointer.col());
				pair<IntPair, Real*> item(key, pointer.value());
				m_coeffMap.insert(item);
			}

			assert(m_coeffMap.size() == m_vpointTriplets.size());
		}

		void Clear()
		{
			m_coeffMap.clear();
			m_vvalueTriplets.clear();
			m_vpointTriplets.clear();
			m_msparseMatrix.setZero();
		}

		inline bool HasTripletData() { return this->m_vvalueTriplets.size() != 0; }
		inline bool HasMatrixData() { return this->m_msparseMatrix.size() != 0; }
		inline bool HasMappingData() { return this->m_coeffMap.size() != 0; }

		inline int GetNumNZUncompressed() { return (int)m_vvalueTriplets.size(); }
		inline int GetNumNZCompressed() { return (int)m_msparseMatrix.nonZeros(); }

	};
}

