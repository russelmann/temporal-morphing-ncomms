//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <stdarg.h> 

#include <PhySim/Utils/Utils.h>

#include <PhySim/Utils/rodriguesRotation.h>

#include <PhySim/Utils/dirent.h>
#include <PhySim/Utils/objload.h>

#pragma warning(disable : 4996) // Disable fopen warning

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	class Logger
	{
	public:

		unique_ptr<ostream> LogSimuStreamRef;
		unique_ptr<ostream> LogTimeStreamRef;
		Logger()
		{
			LogSimuStreamRef = unique_ptr<ostream>(new fstream("logSimu.txt", fstream::out | fstream::trunc));
			LogTimeStreamRef = unique_ptr<ostream>(new fstream("logTime.txt", fstream::out | fstream::trunc));
		}

		virtual ~Logger()
		{
			LogSimuStreamRef.reset();
			LogTimeStreamRef.reset();
		}
	};

	Logger LoggerInstance;


	// 3D MESHES

	bool readTriMesh_Obj(const string& path, MatrixXd& mV, MatrixXi& mF)
	{
		obj::Model inputModel = obj::loadModelFromFile(path);

		int numV = (int) inputModel.vertex.size() / 3;

		mV = MatrixXd(numV, 3);
		for (int i = 0; i < numV; ++i)
		{
			mV(i, 0) = inputModel.vertex[3 * i + 0];
			mV(i, 1) = inputModel.vertex[3 * i + 1];
			mV(i, 2) = inputModel.vertex[3 * i + 2];
		}

		int numF = (int) inputModel.faces.begin()->second.size()/3;
		mF.resize(numF, 3);
		for (int i = 0; i < numF; ++i)
		{
			mF(i, 0) = inputModel.faces.begin()->second[3*i + 0];
			mF(i, 1) = inputModel.faces.begin()->second[3*i + 1];
			mF(i, 2) = inputModel.faces.begin()->second[3*i + 2];
		}

		return true;
	}

	bool writeTriMesh_Obj(const string& path, const MatrixXd& mV, const MatrixXi& mF)
	{
		ostringstream objString;

		objString << "# PhySim OBJ export" << "\n";

		objString << "\n# Vertices" << "\n";

		for (int i = 0; i < mV.rows(); ++i)
		{
			objString << "v";
			objString << " " << mV.row(i).x();
			objString << " " << mV.row(i).y();
			objString << " " << mV.row(i).z();
			objString << "\n";
		}

		objString << "\n# Faces" << "\n";

		for (int i = 0; i < mF.rows(); ++i)
		{
			objString << "f";
			objString << " " << mF.row(i).x() + 1;
			objString << " " << mF.row(i).y() + 1;
			objString << " " << mF.row(i).z() + 1;
			objString << "\n";
		}
		
		ofstream objWriter;
		objWriter.open(path, ios::out);
		objWriter << objString.str();
		objWriter.close();

		return true;
	}

	// FILE SYSTEM
	
	bool listDirectoryFiles(const string& dirPath, vector<string>& vfiles)
	{
#ifdef _WIN32
		DIR* dir;

		if ((dir = opendir(dirPath.c_str())) != NULL) 
		{
			struct dirent* ent;

			while ((ent = readdir(dir)) != NULL) 
			{
				if (string(ent->d_name) == ".") continue;
				if (string(ent->d_name) == "..") continue;

				vfiles.push_back(ent->d_name);
			}

			closedir(dir);
			
			return true;
		}
		else
		{
			return false;
		}
#else
      return false;
#endif
	}

	// STRINGS

	void splitString(const string& input, const string& delim, vector<string>& vs)
	{
		std::string s = input;
		std::string D = delim;

		vs.clear();

		size_t pos = 0;
		std::string token;
		while ((pos = s.find(D)) != std::string::npos) 
		{
			vs.push_back(s.substr(0, pos));
			s.erase(0, pos + D.length());
		}
		vs.push_back(s);
	}

	// LOG/TRACE

	string vectorToString_CSV(const VectorXd& v)
	{
		ostringstream str;
		int np = (int) v.size();
		for (int j = 0; j < np - 1; ++j)
			str << v(j) << ",";
		str << v(np - 1);

		return str.str();
	}

	string matrixToString_CSV(const MatrixXd& M)
	{
		ostringstream str;
		int nr = (int) M.rows();
		int nc = (int) M.cols();
		for (int i = 0; i < nr; ++i)
		{
			for (int j = 0; j < nc - 1; ++j)
				str << M(i, j) << ",";
			str << M(i, nc - 1);

			str << "\n";
		}

		return str.str();
	}

	void logFile(const string& path, const string& content, bool append)
	{
		FILE *fp = NULL;
		if (!append) 
			fp = fopen(path.c_str(), "wt");
		else fp = fopen(path.c_str(), "wa");
		fprintf(fp, "%s", content.c_str());
		fflush(fp);
		fclose(fp);
	}

	bool readFileLines(const string& path, vector<string>& vlines)
	{
		string line;

		ifstream file(path);
		if (file.is_open())
		{
			while (getline(file, line))
				vlines.push_back(line);
			file.close();

			return true;
		}
		else
		{
			return false;
		}
	}

	bool readMatrix_CSV(const string& path, MatrixXd& mV)
	{
		vector<string> vlines;
		if (!readFileLines(path, vlines))
			return false; // Invalid file

		int numRows = (int)vlines.size();
		vector<dVector> vlineValue(numRows);
		for (int i = 0; i < numRows; ++i)
		{
			vector<string> vvalues;
			splitString(vlines[i], ",", vvalues);
			int numValues = (int) vvalues.size();

			vlineValue[i].resize(numValues);
			for (int j = 0; j < numValues; ++j)
				vlineValue[i][j] = atof(vvalues[j].c_str());
		}

		int numCols = (int) vlineValue[0].size();
		mV.resize(numRows, numCols);
		for (int i = 0; i < numRows; ++i)
		for (int j = 0; j < numCols; ++j)
			mV(i, j) = vlineValue[i][j];

		return true;
	}

	bool writeMatrix_CSV(const string& path, MatrixXd& mV)
	{
		logFile(path, matrixToString_CSV(mV), false);

		return true;
	}

	void logSimuStream(ostream* pOstream)
	{
		LoggerInstance.LogSimuStreamRef.reset();

		LoggerInstance.LogSimuStreamRef = unique_ptr<ostream>(pOstream);
	}

	void logTimeStream(ostream* pOstream)
	{
		LoggerInstance.LogTimeStreamRef.reset();

		LoggerInstance.LogTimeStreamRef = unique_ptr<ostream>(pOstream);
	}

	void logSimu(const char *format, ...)
	{
#if _WIN32
		static char msg[60000];
		va_list vl;

		va_start(vl, format);
		vsprintf_s(msg, format, vl);
		va_end(vl);

		(*LoggerInstance.LogSimuStreamRef) << msg;
		LoggerInstance.LogSimuStreamRef->flush();
#else
      static char msg[60000];
      va_list vl;
      
      va_start(vl, format);
      vsprintf(msg, format, vl);
      va_end(vl);
      
      (*LoggerInstance.LogSimuStreamRef) << msg;
      LoggerInstance.LogSimuStreamRef->flush();
#endif
	}

	void logTime(const char *format, ...)
	{
#if _WIN32
		static char msg[1024];
		va_list vl;

		va_start(vl, format);
		vsprintf_s(msg, format, vl);
		va_end(vl);

		(*LoggerInstance.LogTimeStreamRef) << msg;
		LoggerInstance.LogTimeStreamRef->flush();
#else
      static char msg[60000];
      va_list vl;
      
      va_start(vl, format);
      vsprintf(msg, format, vl);
      va_end(vl);
      
      (*LoggerInstance.LogSimuStreamRef) << msg;
      LoggerInstance.LogSimuStreamRef->flush();
#endif
	}

	// EIGEN/STL

	void add3D(int i, const Eigen::Vector3d& v, dVector& vd)
	{
		int offset = 3 * i;
		assert(offset < (int)vd.size());
		for (int i = 0; i < 3; ++i)
			vd[offset + i] += v(i);
	}

	void set3D(int i, const Eigen::Vector3d& v, dVector& vd)
	{
		int offset = 3 * i;
		assert(offset < (int)vd.size());
		for (int i = 0; i < 3; ++i)
			vd[offset + i] = v(i);
	}

	Eigen::Vector3d get3D(int i, const dVector& vd)
	{
		Eigen::Vector3d v;

		int offset = 3 * i;
		assert(offset < (int)vd.size());
		for (int i = 0; i < 3; ++i)
			v(i) = vd[offset + i];

		return v;
	}


	void getSubvector(int i, int nele, const dVector& vin, dVector& vout)
	{
		vout.resize(nele);
		for (int it = 0; it < nele; it++)
			vout[it] = vin[i + it];
	}

	void setSubvector(int i, int nele, const dVector& vin, dVector& vout)
	{
		for (int it = 0; it < nele; it++)
			vout[i + it] = vin[it];
	}

	void getSubvector(int i, int nele, const VectorXd& vin, VectorXd& vout)
	{
		vout.resize(nele);
		for (int it = 0; it < nele; it++)
			vout(it) = vin(i + it);
	}

	void setSubvector(int i, int nele, const VectorXd& vin, VectorXd& vout)
	{
		for (int it = 0; it < nele; it++)
			vout(i + it) = vin(it);
	}

	dVector toSTL(const Eigen::VectorXd& v)
	{
		int n = (int)v.size();
		dVector out(n);
		for (int i = 0; i < n; ++i)
			out[i] = v(i);
		return out;
	}

	iVector toSTL(const Eigen::VectorXi& v)
	{
		int n = (int)v.size();
		iVector out(n);
		for (int i = 0; i < n; ++i)
			out[i] = v(i);
		return out;
	}

	dVector toSTL(const Eigen::Vector3d& v)
	{
		dVector out(3);
		out[0] = v(0);
		out[1] = v(1);
		out[2] = v(2);
		return out;
	}

	Eigen::VectorXd toEigen(const dVector& v)
	{
		int n = (int)v.size();
		VectorXd out(n);
		for (int i = 0; i < n; ++i)
			out(i) = v[i];
		return out;
	}

	Eigen::VectorXi toEigen(const iVector& v)
	{
		int n = (int)v.size();
		VectorXi out(n);
		for (int i = 0; i < n; ++i)
			out(i) = v[i];
		return out;
	}
	
	void assembleSparseMatrix(VectorTd & triplets, unsigned int row, unsigned int col, const Matrix3d& m33)
	{
		triplets.push_back(Triplet<double>(row + 0, col + 0, m33(0, 0)));
		triplets.push_back(Triplet<double>(row + 0, col + 1, m33(0, 1)));
		triplets.push_back(Triplet<double>(row + 0, col + 2, m33(0, 2)));

		triplets.push_back(Triplet<double>(row + 1, col + 0, m33(1, 0)));
		triplets.push_back(Triplet<double>(row + 1, col + 1, m33(1, 1)));
		triplets.push_back(Triplet<double>(row + 1, col + 2, m33(1, 2)));

		triplets.push_back(Triplet<double>(row + 2, col + 0, m33(2, 0)));
		triplets.push_back(Triplet<double>(row + 2, col + 1, m33(2, 1)));
		triplets.push_back(Triplet<double>(row + 2, col + 2, m33(2, 2)));
	}

	void eigenTripletsToSparseMatrix(const VectorTd& triplets, MatrixSd& M)
	{
		M.setFromTriplets(triplets.begin(), triplets.end());
	}

	void eigenSparseMatrixToPointers(const MatrixSd& M, VectorTp& triplets)
	{
		int NZ = M.nonZeros();
		if (triplets.capacity() < NZ)
			triplets.reserve(NZ);
		triplets.clear();

		for (int k = 0; k < M.outerSize(); ++k)
			for (MatrixSd::InnerIterator it(M, k); it; ++it)
				triplets.push_back(Triplet<Real*>(it.row(), it.col(), &(it.valueRef())));
	}

	void eigenSparseMatrixToTriplets(const MatrixSd& M, VectorTd& triplets)
	{
		int NZ = M.nonZeros();
		if (triplets.capacity() < NZ)
			triplets.reserve(NZ);
		triplets.clear();

		for (int k = 0; k < M.outerSize(); ++k)
			for (MatrixSd::InnerIterator it(M, k); it; ++it)
				triplets.push_back(Triplet<Real>(it.row(), it.col(), it.value()));
	}

	void getDiagonalIndices(const VectorTd& vT, int rows, int cols, iVector& vDIdx)
	{
		vDIdx.reserve(rows);

		for (size_t i = 0; i < vT.size(); ++i)
			if (vT[i].row() == vT[i].col())
				vDIdx.push_back((int) i);
	}

	void transformEigenVectorToCUDA(const VectorXd& veigen, vector<float>& vcuda)
	{
		int numVal = (int)veigen.size();

		vcuda.clear();
		vcuda.reserve(numVal);
		for (int i = 0; i < numVal; ++i)
			vcuda.push_back((float)veigen(i));
	}

	void transformEigenVectorToCUDA(const VectorXd& veigen, vector<double>& vcuda)
	{
		int numVal = (int)veigen.size();

		vcuda.clear();
		vcuda.reserve(numVal);
		for (int i = 0; i < numVal; ++i)
			vcuda.push_back((float)veigen(i));
	}

	void transformCUDAToEigenVector(const vector<double>& vcuda, VectorXd& veigen)
	{
		int numVal = (int)vcuda.size();

		veigen.resize(numVal);
		for (int i = 0; i < numVal; ++i)
			veigen(i) = vcuda[i];
	}

	void transformCUDAToEigenVector(const vector<float>& vcuda, VectorXd& veigen)
	{
		int numVal = (int)vcuda.size();

		veigen.resize(numVal);
		for (int i = 0; i < numVal; ++i)
			veigen(i) = vcuda[i];
	}

	void transformTripletValuesToCUDA(const VectorTd& vT, vector<float>& vals)
	{
		int numNZ = (int)vT.size();
		assert((int)vals.size() == numNZ);
		for (int i = 0; i < numNZ; ++i)
			vals[i] = (float) vT[i].value();
	}

	void transformTripletValuesToCUDA(const VectorTd& vT, vector<double>& vals)
	{
		int numNZ = (int)vT.size();
		assert((int)vals.size() == numNZ);
		for (int i = 0; i < numNZ; ++i)
			vals[i] = (float)vT[i].value();
	}

	void transformCUDAToTripletValues(const vector<float>& vals, VectorTd& vT)
	{
		int numNZ = (int)vals.size();
		assert((int)vT.size() == numNZ);
		for (int i = 0; i < numNZ; ++i)
			vT[i] = Triplet<Real>(vT[i].row(), vT[i].col(), vals[i]);
	}

	void transformCUDAToTripletValues(const vector<double>& vals, VectorTd& vT)
	{
		int numNZ = (int)vals.size();
		assert((int)vT.size() == numNZ);
		for (int i = 0; i < numNZ; ++i)
			vT[i] = Triplet<Real>(vT[i].row(), vT[i].col(), vals[i]);
	}

	void transformEigenTripletsToCUDA(const VectorTd& triplets, vector<int>& rows, vector<int>& cols, vector<float>& vals)
	{
		int numNZ = (int)triplets.size();
		rows.reserve(numNZ);
		cols.reserve(numNZ);
		vals.reserve(numNZ);
		for (int i = 0; i < numNZ; ++i)
		{
			const Triplet<double>& t = triplets[i];
			rows.push_back(t.row());
			cols.push_back(t.col());
			vals.push_back((float)t.value());
		}
	}

	void transformCUDAToEigenTriplets(const vector<int>& rows, const vector<int>& cols, const vector<float>& vals, VectorTd& triplets)
	{
		int numNZ = (int)rows.size();
		triplets.reserve(numNZ);
		for (int i = 0; i < numNZ; ++i)
			triplets.push_back(Triplet<double>(rows[i], cols[i], vals[i]));
	}
	
	// GEOMETRY

	Matrix3d transformEulerToMatrix(const Vector3d& ve)
	{
		Matrix3d Rx;
		Matrix3d Ry;
		Matrix3d Rz;
		
		Rx(0, 0) = 1; Rx(0, 1) = 0; Rx(0, 2) = 0;
		Rx(1, 0) = 0; Rx(1, 1) = cos(ve[0]); Rx(1, 2) = -sin(ve[0]);
		Rx(2, 0) = 0; Rx(2, 1) = sin(ve[0]); Rx(2, 2) = cos(ve[0]);

		Ry(0, 0) = cos(ve[1]); Ry(0, 1) = 0; Ry(0, 2) = sin(ve[1]);
		Ry(1, 0) = 0; Ry(1, 1) = 1; Ry(1, 2) = 0;
		Ry(2, 0) = -sin(ve[1]); Ry(2, 1) = 0; Ry(2, 2) = cos(ve[1]);

		Rz(0, 0) = cos(ve[2]); Rz(0, 1) = -sin(ve[2]); Rz(0, 2) = 0;
		Rz(1, 0) = sin(ve[2]); Rz(1, 1) = cos(ve[2]); Rz(1, 2) = 0;
		Rz(2, 0) = 0; Rz(2, 1) = 0; Rz(2, 2) = 1;

		return Rz*Ry*Rx;
	}

	Vector3d transformMatrixToEuler(const Matrix3d& mR)
	{
		Real rotXangle = atan2(-mR.row(1).z(), mR.row(2).z());
		Real cosYangle = sqrt(pow(mR.row(0).x(), 2) + pow(mR.row(0).y(), 2));
		Real rotYangle = atan2(mR.row(0).z(), cosYangle);
		Real sinXangle = sin(rotXangle);
		Real cosXangle = cos(rotXangle);
		Real rotZangle = atan2(cosXangle * mR.row(1).x() + sinXangle * mR.row(2).x(), cosXangle * mR.row(1).y() + sinXangle * mR.row(2).y());
		return Vector3d(rotXangle, rotYangle, rotZangle);
	}

	void computeBoundingBox(MatrixXd& mV, Vector3d& vmin, Vector3d& vmax)
	{
		vmin = Vector3d::Ones()*HUGE_VAL;
		vmax = -Vector3d::Ones()*HUGE_VAL;

		for (int i = 0; i < mV.rows(); ++i)
		{
			if (mV.row(i).x() < vmin.x()) 
				vmin.x() = mV.row(i).x();

			if (mV.row(i).y() < vmin.y()) 
				vmin.y() = mV.row(i).y();

			if (mV.row(i).z() < vmin.z()) 
				vmin.z() = mV.row(i).z();

			if (mV.row(i).x() > vmax.x()) 
				vmax.x() = mV.row(i).x();

			if (mV.row(i).y() > vmax.y()) 
				vmax.y() = mV.row(i).y();

			if (mV.row(i).z() > vmax.z()) 
				vmax.z() = mV.row(i).z();
		}
	}
	
	Frame3d parallelTransport(const Frame3d& F, const Vector3d& t2)
	{
		Frame3d Fin = F;

		Frame3d Fout;
		Fout.tan = t2.normalized();
		Fout.nor = parallelTransport(F.nor, F.tan, Fout.tan);
		Fout.bin = Fout.tan.cross(Fout.nor);

		return Fout;
	}

	Vector3d parallelTransport(const Vector3d& u, const Vector3d& t1, const Vector3d& t2)
	{
		double d = t1.dot(t2);

		Vector3d b = t1.cross(t2);

		return u * d + b * (b.dot(u)) * 1 / (1 + d) + b.cross(u);
	};

	Vector3d parallelTransportNormalized(const Vector3d& u, const Vector3d& t1, const Vector3d& t2)
	{
		Vector3d te = t1.normalized();
		Vector3d tf = t2.normalized();

		double d = te.dot(tf);

		Vector3d b = te.cross(tf);

		return u * d + b * (b.dot(u)) * 1 / (1 + d) + b.cross(u);
	};

	Vector3d computeCurvatureBinormal(const Vector3d& e0, const Vector3d& e1)
	{
		Vector3d t0 = e0.normalized();
		Vector3d t1 = e1.normalized();
		return t0.cross(t1) * 2 * (1 / (1 + t0.dot(t1)));
	}

	void orthonormalizeFrame(Frame3d& F)
	{
		F.tan.normalize();
		F.bin = F.tan.cross(F.nor).normalized();
		F.nor = F.bin.cross(F.tan).normalized();
	}

	Frame3d rotateFrame(const Frame3d& frame, Matrix3d& R)
	{
		Frame3d newFrame;
		newFrame.bin = R * frame.bin;
		newFrame.nor = R * frame.nor;
		newFrame.tan = R * frame.tan;
		return newFrame;
	}

	Frame3d rotateFrame(const Frame3d& vF, Real angle)
	{
		Frame3d rot;
		rot.tan = vF.tan;
		Vector3d nor0 = vF.nor;
		Vector3d bin0 = vF.bin;
		rodriguesRotation(rot.tan.data(), angle, nor0.data(), rot.nor.data());
		rodriguesRotation(rot.tan.data(), angle, bin0.data(), rot.bin.data());

		return rot;
	}

	void rotateFrame(const vector<Frame3d>& vF, const dVector& va, vector<Frame3d>& vFRot)
	{
		int ne = (int)vF.size();
		assert(ne == va.size());
		assert(ne == vFRot.size());
		for (int i = 0; i < ne; i++)
		{
			vFRot[i] = rotateFrame(vF[i], va[i]);
		}
	}

	bool lineLineIntersect(
		const dVector& p1, const dVector& q1,
		const dVector& p2, const dVector& q2,
		double& s, double& t, double tol)
	{
		double EPS1 = 1e-6;
		double EPS2 = tol;

		double d1[3], d2[3], r[3], a, e, f;
		double c1[3], c2[3];

		d1[0] = q1[0] - p1[0];
		d1[1] = q1[1] - p1[1];
		d1[2] = q1[2] - p1[2];

		d2[0] = q2[0] - p2[0];
		d2[1] = q2[1] - p2[1];
		d2[2] = q2[2] - p2[2];

		r[0] = p1[0] - p2[0];
		r[1] = p1[1] - p2[1];
		r[2] = p1[2] - p2[2];

		a = d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2];
		e = d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2];
		f = d2[0] * r[0] + d2[1] * r[1] + d2[2] * r[2];

		// Check if either or both segments degenerate into points
		//
		if ((a <= EPS1) && (e <= EPS1))
		{
			s = t = 0.0f;
			c1[0] = p1[0]; c1[1] = p1[1]; c1[2] = p1[2];
			c2[0] = p2[0]; c2[1] = p2[1]; c2[2] = p2[2];
			return (((c1[0] - c2[0])*(c1[0] - c2[0]) + (c1[1] - c2[1])*(c1[1] - c2[1]) + (c1[2] - c2[2])*(c1[2] - c2[2]))) < EPS1;
		}

		if (a <= EPS1)
		{
			// First segment degenerates into a point
			//
			s = 0.0f;
			t = f / e;
			if (t<0.0f) t = 0.0f;
			if (t>1.0f) t = 1.0f;
		}
		else
		{
			double c = d1[0] * r[0] + d1[1] * r[1] + d1[2] * r[2];

			if (e <= EPS1)
			{
				// Second segment degenerates into a point
				//
				t = 0.0f;
				s = -c / a;
				if (s<0.0f) s = 0.0f;
				if (s>1.0f) s = 1.0f;
			}
			else
			{
				// Nondegenerate case
				//
				double b = d1[0] * d2[0] + d1[1] * d2[1] + d1[2] * d2[2];
				double denom = a*e - b*b;

				if (denom != 0.0f)
				{
					s = (b*f - c*e) / denom;
					if (s<0.0f) s = 0.0f;
					if (s>1.0f) s = 1.0f;
				}
				else s = 0.0f;

				double tnom = b*s + f;
				if (tnom < 0.0f)
				{
					t = 0.0f;
					s = -c / a;
					if (s<0.0f) s = 0.0f;
					if (s>1.0f) s = 1.0f;
				}
				else if (tnom > e)
				{
					t = 1.0f;
					s = (b - c) / a;
					if (s<0.0f) s = 0.0f;
					if (s>1.0f) s = 1.0f;
				}
				else
					t = tnom / e;
			}
		}

		c1[0] = p1[0] + d1[0] * s;
		c1[1] = p1[1] + d1[1] * s;
		c1[2] = p1[2] + d1[2] * s;

		c2[0] = p2[0] + d2[0] * t;
		c2[1] = p2[1] + d2[1] * t;
		c2[2] = p2[2] + d2[2] * t;

		double dist = sqrt((c1[0] - c2[0])*(c1[0] - c2[0]) + (c1[1] - c2[1])*(c1[1] - c2[1]) + (c1[2] - c2[2])*(c1[2] - c2[2]));

		return dist < EPS2 && s >= 0.0 && s <= 1.0 && t >= 0.0 && t <= 1.0;
	}

	// INTERPOLATION
	
	void hermiteInterpolation(const vector<Vector3d>& vin, vector<Vector3d>& vout, bool loop, int outnp, const Vector3d& tanIni, const Vector3d& tanEnd)
	{
		int innp = (int)vin.size();

		dVector vsin;
		for (int i = 0; i < innp - 1; ++i)
			vsin.push_back((vin[i + 1] - vin[i]).norm());

		int inns = (int)vsin.size();

		double sT = 0;
		for (int i = 0; i < inns; ++i)
			sT += vsin[i]; // Sum length

		int outns = outnp - 1; // Out number of domains
		double outs = sT / outns; // Out domain length

								  // Must match ends
		vout.resize(outnp);
		vout.back() = vin.back();
		vout.front() = vin.front();

		int sit = 0;
		int nxtp = 1;
		double sums = 0;
		double nxts = outs;
		for (int i = 1; i < outnp - 1; i++)
		{
			double s = vsin[sit];
			while (i*outs > sums + s)
			{
				sit++;
				sums += s;
				s = vsin[sit];
			}

			double t = (nxts - sums) / s;
			double t2 = t*t;
			double t3 = t*t*t;

			const Vector3d& p0 = vin[sit]; // In point 0
			const Vector3d& p1 = vin[sit + 1]; // In point 1
			Vector3d ts = (p1 - p0)*(1.0 / s); // In tan

											   // Get tangents
			Vector3d t0, t1;
			if (inns == 1)
			{
				if (tanIni.norm() > 0.0)
				{
					t0 = tanIni;
				}
				else t0 = ts;

				if (tanEnd.norm() > 0.0)
				{
					t1 = tanEnd;
				}
				else t1 = ts;
			}
			else if (sit == 0 && !loop) // First domain (no loop)
			{
				if (tanIni.norm() > 0.0)
				{
					t0 = tanIni;
				}
				else t0 = ts;

				const Vector3d& pn = vin[sit + 2];
				t1 = 0.5*(ts + (pn - p1)*(1.0 / s));
			}
			else if (sit == inns - 1 && !loop) // Last domain (no loop)
			{
				if (tanEnd.norm() > 0.0)
				{
					t1 = tanEnd;
				}
				else t1 = ts;

				const Vector3d& pp = vin[sit - 1];
				t0 = 0.5*(ts + (p0 - pp)*(1.0 / s));
			}
			else if (sit == 0 && loop) // First domain (loop)
			{
				const Vector3d& pp = vin[innp - 2];
				const Vector3d& pn = vin[sit + 2];
				t0 = 0.5*(ts + (p0 - pp)*(1.0 / s)); // Mean
				t1 = 0.5*(ts + (pn - p1)*(1.0 / s)); // Mean
			}
			else if (sit == inns - 1 && loop) // Last domain (loop)
			{
				const Vector3d& pp = vin[sit - 1];
				const Vector3d& pn = vin[0];
				t0 = 0.5*(ts + (p0 - pp)*(1.0 / s)); // Mean
				t1 = 0.5*(ts + (pn - p1)*(1.0 / s)); // Mean
			}
			else
			{
				const Vector3d& pp = vin[sit - 1];
				const Vector3d& pn = vin[sit + 2];
				t0 = 0.5*(ts + (p0 - pp)*(1.0 / s)); // Mean
				t1 = 0.5*(ts + (pn - p1)*(1.0 / s)); // Mean
			}

			// Hermite interpolation

			Vector3d p = p0*(2 * t3 - 3 * t2 + 1) +
				t0*(t3 - 2 * t2 + t)*s +
				p1*(-2 * t3 + 3 * t2) +
				t1*(t3 - t2)*s;

			vout[nxtp++] = p; // Point

			nxts += outs;
		}
	}

	void hermiteInterpolation(const vector<Vector2d>& vin, const dVector& vsin, vector<Vector2d>& vout, bool loop, int outnp)
	{
		int innp = (int)vin.size();
		int inns = (int)vsin.size();

		assert(inns == innp - 1);

		double sT = 0;
		for (int i = 0; i < inns; ++i)
			sT += vsin[i]; // Sum length

		int outns = outnp - 1; // Out number of domains
		double outs = sT / outns; // Out domain length

								  // Must match ends
		vout.resize(outnp);
		vout.back() = vin.back();
		vout.front() = vin.front();

		int sit = 0;
		int nxtp = 1;
		double sums = 0;
		double nxts = outs;
		for (int i = 1; i < outnp - 1; i++)
		{
			double s = vsin[sit];
			while (i*outs > sums + s)
			{
				sit++;
				sums += s;
				s = vsin[sit];
			}

			double t = (nxts - sums) / s;
			double t2 = t*t;
			double t3 = t*t*t;

			const Vector2d& p0 = vin[sit]; // In point 0
			const Vector2d& p1 = vin[sit + 1]; // In point 1
			Vector2d ts = (p1 - p0)*(1.0 / s); // In tan

											   // Get tangents
			Vector2d t0, t1;
			if (inns == 1)
			{
				t0 = ts;
				t1 = ts;
			}
			else if (sit == 0 && !loop) // First domain (no loop)
			{
				t0 = ts;
				const Vector2d& pn = vin[sit + 2];
				t1 = 0.5*(ts + (pn - p1)*(1.0 / s));
			}
			else if (sit == inns - 1 && !loop) // Last domain (no loop)
			{
				t1 = ts;
				const Vector2d& pp = vin[sit - 1];
				t0 = 0.5*(ts + (p0 - pp)*(1.0 / s));
			}
			else if (sit == 0 && loop) // First domain (loop)
			{
				const Vector2d& pp = vin[innp - 2];
				const Vector2d& pn = vin[sit + 2];
				t0 = 0.5*(ts + (p0 - pp)*(1.0 / s)); // Mean
				t1 = 0.5*(ts + (pn - p1)*(1.0 / s)); // Mean
			}
			else if (sit == inns - 1 && loop) // Last domain (loop)
			{
				const Vector2d& pp = vin[sit - 1];
				const Vector2d& pn = vin[0];
				t0 = 0.5*(ts + (p0 - pp)*(1.0 / s)); // Mean
				t1 = 0.5*(ts + (pn - p1)*(1.0 / s)); // Mean
			}
			else
			{
				const Vector2d& pp = vin[sit - 1];
				const Vector2d& pn = vin[sit + 2];
				t0 = 0.5*(ts + (p0 - pp)*(1.0 / s)); // Mean
				t1 = 0.5*(ts + (pn - p1)*(1.0 / s)); // Mean
			}

			// Hermite interpolation

			Vector2d p = p0*(2 * t3 - 3 * t2 + 1) +
				t0*(t3 - 2 * t2 + t)*s +
				p1*(-2 * t3 + 3 * t2) +
				t1*(t3 - t2)*s;

			vout[nxtp++] = p; // Point

			nxts += outs;
		}
	}

	void hermiteInterpolation(const vector<Frame3d>& vin, const dVector& vsin, vector<Frame3d>& vout, bool loop, int outnp)
	{
		int innp = (int)vin.size();
		int inns = (int)vsin.size();

		assert(inns == innp - 1);

		double sT = 0;
		for (int i = 0; i < inns; ++i)
			sT += vsin[i]; // Sum length

		int outns = outnp - 1; // Out number of domains
		double outs = sT / outns; // Out domain length

								  // Must match ends
		vout.resize(outnp);
		vout.back() = vin.back();
		vout.front() = vin.front();

		int sit = 0;
		int nxtp = 1;
		double sums = 0;
		double nxts = outs;
		for (int i = 1; i < outnp - 1; i++)
		{
			double s = vsin[sit];
			while (i*outs > sums + s)
			{
				sit++;
				sums += s;
				s = vsin[sit];
			}

			double t = (nxts - sums) / s;
			double t2 = t*t;
			double t3 = t*t*t;

			const Frame3d& p0 = vin[sit]; // In point 0
			const Frame3d& p1 = vin[sit + 1]; // In point 1
			Frame3d ts; // In tan
			ts.tan = (p1.tan - p0.tan)*(1.0 / s);
			ts.nor = (p1.nor - p0.nor)*(1.0 / s);
			ts.bin = (p1.bin - p0.bin)*(1.0 / s);

			// Get tangents
			Frame3d t0, t1;
			if (inns == 1)
			{
				t0 = ts;
				t1 = ts;
			}
			else if (sit == 0 && !loop) // First domain (no loop)
			{
				t0 = ts;
				const Frame3d& pn = vin[sit + 2];
				t1.tan = 0.5*(ts.tan + (pn.tan - p1.tan)*(1.0 / s));
				t1.nor = 0.5*(ts.nor + (pn.nor - p1.nor)*(1.0 / s));
				t1.bin = 0.5*(ts.bin + (pn.bin - p1.bin)*(1.0 / s));
			}
			else if (sit == inns - 1 && !loop) // Last domain (no loop)
			{
				t1 = ts;
				const Frame3d& pp = vin[sit - 1];
				t0.tan = 0.5*(ts.tan + (p0.tan - pp.tan)*(1.0 / s));
				t0.nor = 0.5*(ts.nor + (p0.nor - pp.nor)*(1.0 / s));
				t0.bin = 0.5*(ts.bin + (p0.bin - pp.bin)*(1.0 / s));
			}
			else if (sit == 0 && loop) // First domain (loop)
			{
				const Frame3d& pp = vin[innp - 2];
				const Frame3d& pn = vin[sit + 2];
				t0.tan = 0.5*(ts.tan + (p0.tan - pp.tan)*(1.0 / s)); // Mean
				t0.nor = 0.5*(ts.nor + (p0.nor - pp.nor)*(1.0 / s)); // Mean
				t0.bin = 0.5*(ts.bin + (p0.bin - pp.bin)*(1.0 / s)); // Mean

				t1.tan = 0.5*(ts.tan + (pn.tan - p1.tan)*(1.0 / s)); // Mean
				t1.nor = 0.5*(ts.nor + (pn.nor - p1.nor)*(1.0 / s)); // Mean
				t1.bin = 0.5*(ts.bin + (pn.bin - p1.bin)*(1.0 / s)); // Mean
			}
			else if (sit == inns - 1 && loop) // Last domain (loop)
			{
				const Frame3d& pp = vin[sit - 1];
				const Frame3d& pn = vin[0];
				t0.tan = 0.5*(ts.tan + (p0.tan - pp.tan)*(1.0 / s)); // Mean
				t0.nor = 0.5*(ts.nor + (p0.nor - pp.nor)*(1.0 / s)); // Mean
				t0.bin = 0.5*(ts.bin + (p0.bin - pp.bin)*(1.0 / s)); // Mean

				t1.tan = 0.5*(ts.tan + (pn.tan - p1.tan)*(1.0 / s)); // Mean
				t1.nor = 0.5*(ts.nor + (pn.nor - p1.nor)*(1.0 / s)); // Mean
				t1.bin = 0.5*(ts.bin + (pn.bin - p1.bin)*(1.0 / s)); // Mean
			}
			else
			{
				const Frame3d& pp = vin[sit - 1];
				const Frame3d& pn = vin[sit + 2];
				t0.tan = 0.5*(ts.tan + (p0.tan - pp.tan)*(1.0 / s)); // Mean
				t0.nor = 0.5*(ts.nor + (p0.nor - pp.nor)*(1.0 / s)); // Mean
				t0.bin = 0.5*(ts.bin + (p0.bin - pp.bin)*(1.0 / s)); // Mean

				t1.tan = 0.5*(ts.tan + (pn.tan - p1.tan)*(1.0 / s)); // Mean
				t1.nor = 0.5*(ts.nor + (pn.nor - p1.nor)*(1.0 / s)); // Mean
				t1.bin = 0.5*(ts.bin + (pn.bin - p1.bin)*(1.0 / s)); // Mean
			}

			// Hermite interpolation

			Frame3d p;

			p.tan = p0.tan*(2 * t3 - 3 * t2 + 1) +
				t0.tan*(t3 - 2 * t2 + t)*s +
				p1.tan*(-2 * t3 + 3 * t2) +
				t1.tan*(t3 - t2)*s;

			p.bin = p0.bin*(2 * t3 - 3 * t2 + 1) +
				t0.bin*(t3 - 2 * t2 + t)*s +
				p1.bin*(-2 * t3 + 3 * t2) +
				t1.bin*(t3 - t2)*s;

			p.nor = p0.nor*(2 * t3 - 3 * t2 + 1) +
				t0.nor*(t3 - 2 * t2 + t)*s +
				p1.nor*(-2 * t3 + 3 * t2) +
				t1.nor*(t3 - t2)*s;

			orthonormalizeFrame(p);

			vout[nxtp++] = p; // Point

			nxts += outs;
		}
	}

	void hermiteInterpolation(const vector<Frame3d>& vin, const dVector& vsin, vector<Frame3d>& vout, bool loop, double outTol)
	{
		assert(vin.size() - 1 == vsin.size());

		double length = 0;
		int nPoint = (int)vin.size();
		for (int i = 0; i < nPoint - 1; ++i)
			length += vsin[i]; // Add total

		int outnp = (int) ceil(length / outTol) + 1; // Output
		hermiteInterpolation(vin, vsin, vout, loop, outnp);
	}

	void hermiteInterpolation(const vector<Vector3d>& vin, vector<Vector3d>& vout, bool loop, double outTol, const Vector3d& tanIni, const Vector3d& tanEnd)
	{
		double length = 0;
		int nPoint = (int)vin.size();
		for (int i = 0; i < nPoint - 1; ++i)
			length += (vin[i + 1] - vin[i]).norm();
		int outnp = (int) ceil(length / outTol) + 1; // Output

		hermiteInterpolation(vin, vout, loop, outnp, tanIni, tanEnd);
	}

	void hermiteInterpolation(const vector<Vector2d>& vin, const dVector& vsin, vector<Vector2d>& vout, bool loop, double outTol)
	{
		assert(vin.size() - 1 == vsin.size());

		double length = 0;
		int nPoint = (int)vin.size();
		for (int i = 0; i < nPoint - 1; ++i)
			length += vsin[i]; // Add total

		int outnp = (int) ceil(length / outTol) + 1; // Output
		hermiteInterpolation(vin, vsin, vout, loop, outnp);
	}

	// OTHER MATH UTILS

	Real computePolarDecomposition(const Matrix3d& M, Matrix3d& Q, Matrix3d& S, Real tol)
	{
		double det = M.determinant();
		if (abs(det) < 1e-9)
		{
			return -1.0;
		}

		Matrix3d Mk;
		Matrix3d Ek;
		double M_oneNorm, M_infNorm, E_oneNorm;

		// Initialize
		Mk = M.transpose();
		M_oneNorm = Mk.lpNorm<1>();
		M_infNorm = Mk.lpNorm<Eigen::Infinity>();

		// Iterate
		do
		{
			Matrix3d MadjTk;

			Vector3d row0 = Mk.transpose().col(0);
			Vector3d row1 = Mk.transpose().col(1);
			Vector3d row2 = Mk.transpose().col(2);

			MadjTk.row(0) = row1.cross(row2);
			MadjTk.row(1) = row2.cross(row0);
			MadjTk.row(2) = row0.cross(row1);

			det = Mk(0, 0) * MadjTk(0, 0) + Mk(0, 1) * MadjTk(0, 1) + Mk(0, 2) * MadjTk(0, 2);
			if (abs(det) < 1e-9)
			{
				return -1.0;
			}

			double MadjT_one = MadjTk.lpNorm<1>();
			double MadjT_inf = MadjTk.lpNorm<Eigen::Infinity>();

			double gamma = sqrt(sqrt((MadjT_one*MadjT_inf) / (M_oneNorm*M_infNorm)) / fabs(det));

			double g2 = 0.5 / (gamma * det);
			double g1 = gamma * 0.5;

			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
				{
					Ek(i, j) = Mk(i, j); // Modify previous Mk
					Mk(i, j) = g1 * Mk(i, j) + g2 * MadjTk(i, j);
					Ek(i, j) -= Mk(i, j);
				}

			E_oneNorm = Ek.lpNorm<1>();
			M_oneNorm = Mk.lpNorm<1>();
			M_infNorm = Mk.lpNorm<Eigen::Infinity>();
		} while (E_oneNorm > M_oneNorm * tol);

		Q = Mk.transpose();
		S = Q.inverse()*M;

		return 1.0;
	}

	void computeBestRigidTransform(const MatrixXd& mCloud0, const MatrixXd& mCloud1, MatrixXd& mR, VectorXd& vt)
	{
		assert(mCloud0.rows() == mCloud1.rows());
		assert(mCloud0.cols() == mCloud1.cols());

		int N = mCloud0.rows();
		int M = mCloud0.cols();

		// Compute centroids

		VectorXd vc0 = mCloud0.row(0);
		VectorXd vc1 = mCloud1.row(0);
		for (int i = 1; i < N; ++i)
		{
			vc0 += mCloud0.row(i);
			vc1 += mCloud1.row(i);
		}
		vc0 /= (Real)N;
		vc1 /= (Real)N;

		MatrixXd mVectors0 = mCloud0.rowwise() - vc0.transpose();
		MatrixXd mVectors1 = mCloud1.rowwise() - vc1.transpose();
		MatrixXd mS = mVectors0.transpose()*mVectors1;

		JacobiSVD<MatrixXd> svd(mS, ComputeFullU | ComputeFullV);
		const MatrixXd& mU = svd.matrixU();
		const MatrixXd& mV = svd.matrixV();
		MatrixXd mI = MatrixXd::Identity(M, M);
		mI(M - 1, M - 1) = (mV*mU.transpose()).determinant();
		mR = mV*mI*mU.transpose();
		vt = vc1 - mR*vc0;
	}

	void eulerRotation(const double *eulerAngles, const double *vectorOriginal, double cgret[3])
	{
		double t1;
		double t10;
		double t11;
		double t13;
		double t14;
		double t17;
		double t2;
		double t22;
		double t26;
		double t3;
		double t4;
		double t6;
		double t8;
		double t9;
		double vectorRotated[3];
		vectorRotated[0] = 0;
		vectorRotated[1] = 0;
		vectorRotated[2] = 0;
		t1 = eulerAngles[2];
		t2 = cos(t1);
		t3 = eulerAngles[1];
		t4 = cos(t3);
		t6 = vectorOriginal[0];
		t8 = sin(t3);
		t9 = t2 * t8;
		t10 = eulerAngles[0];
		t11 = sin(t10);
		t13 = sin(t1);
		t14 = cos(t10);
		t17 = vectorOriginal[1];
		t22 = vectorOriginal[2];
		vectorRotated[0] = t2 * t4 * t6 + (t9 * t11 - t13 * t14) * t17 + (t9 * t14 + t13 * t11) * t22;
		t26 = t13 * t8;
		vectorRotated[1] = t13 * t4 * t6 + (t26 * t11 + t2 * t14) * t17 + (t26 * t14 - t2 * t11) * t22;
		vectorRotated[2] = -t8 * t6 + t4 * t11 * t17 + t4 * t14 * t22;
		cgret[0] = vectorRotated[0];
		cgret[1] = vectorRotated[1];
		cgret[2] = vectorRotated[2];
	}

}

