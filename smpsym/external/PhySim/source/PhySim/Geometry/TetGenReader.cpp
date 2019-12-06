//==========================================================
//
//	PhySim library. Generic library for physical simulation.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#include <PhySim/Models/Model_Inflatable.h>

#include <PhySim/Geometry/TetGenReader.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	struct TetGenNode {
		int index;
		double x;
		double y;
		double z;

		TetGenNode(int i, double x, double y, double z)
		{
			index = i;
			this->x = x;
			this->y = y;
			this->z = z;
		}

		bool Equals(TetGenNode a)
		{
			double distanceSquared =
				std::pow(x - a.x, 2) +
				std::pow(y - a.y, 2) +
				std::pow(z - a.z, 2);

			return distanceSquared < 0.0001; // 0.001 distance tolerated
		}

		Eigen::Vector3d GetPosition() const
		{
			return Eigen::Vector3d(this->x, this->y, this->z);
		}

		bool Equals(TetGenNode a, double precision)
		{
			double distanceSquared =
				std::pow(this->x - a.x, 2) +
				std::pow(this->y - a.y, 2) +
				std::pow(this->z - a.z, 2);

			return distanceSquared < precision;
		}
	};

	struct TetGenFace
	{
		int A;
		int B;
		int C;

		TetGenFace(int A, int B, int C)
		{
			this->A = A;
			this->B = B;
			this->C = C;
		}
	};

	struct TetGenEdge
	{
		int A;
		int B;

		TetGenEdge(int A, int B)
		{
			this->A = A;
			this->B = B;
		}
	};

	struct TetGenTetra
	{
		int A;
		int B;
		int C;
		int D;

		TetGenTetra(int A, int B, int C, int D)
		{
			this->A = A;
			this->B = B;
			this->C = C;
			this->D = D;
		}
	};


	// Data structures used only internally -- END

	// Function declarations:
	vector<string> ReadFile(string fileName) {
		string line;
		vector<string> lines;

		std::ifstream myfile(fileName);
		if (myfile.is_open())
		{
			while (getline(myfile, line))
			{
				lines.push_back(line);
			}
			myfile.close();
		}

		return lines;
	}

	bool TetGenReader::readPoly(string file, MatrixXd& mNodes, MatrixXi& mTetra, MatrixXi& mEdges, vector<MatrixXi>& vmFaces)
	{
		vector<TetGenNode> vnodes = readNodes(file + ".1.node");
		vector<TetGenEdge> vedges = readEdges(file + ".1.edge");
		vector<TetGenTetra> vtetra = readTetra(file + ".1.ele");

		vector<vector<TetGenFace> > vfaces;
		readFacesPoly(file + ".1.face", vfaces);

		convertIndexing(vnodes, vtetra, vedges, vfaces);

		int numNodes = (int)vnodes.size();
		int numEdges = (int)vedges.size();
		int numTetra = (int)vtetra.size();
		int numFaceSets = (int)vfaces.size();
		mNodes.resize(numNodes, 3);
		mEdges.resize(numEdges, 2);
		mTetra.resize(numTetra, 4);
		vmFaces.resize(numFaceSets);

		for (int i = 0; i < numNodes; ++i)
			mNodes.row(i) = vnodes[i].GetPosition();

		for (int i = 0; i < numEdges; ++i)
		{
			mEdges(i, 0) = vedges[i].A;
			mEdges(i, 1) = vedges[i].B;
		}
		for (int i = 0; i < numTetra; ++i)
		{
			mTetra(i, 0) = vtetra[i].A;
			mTetra(i, 1) = vtetra[i].B;
			mTetra(i, 2) = vtetra[i].C;
			mTetra(i, 3) = vtetra[i].D;
		}
		for (int i = 0; i < numFaceSets; ++i)
		{
			int numFaces = (int)vfaces[i].size();
			
			vmFaces[i].resize(numFaces, 3);
			for (int j = 0; j < numFaces; ++j)
			{
				vmFaces[i](j, 0) = vfaces[i][j].A;
				vmFaces[i](j, 1) = vfaces[i][j].B;
				vmFaces[i](j, 2) = vfaces[i][j].C;
			}
		}

		return true;
	}

	void TetGenReader::readFacesPoly(string path, vector<vector<TetGenFace> >& vfaces)
	{
		vector<string> vlines;
		map<int, int> msetMap;
		int	countSet = 0;

		readFileLines(path, vlines);

		// skip first line - not relevant to nodes
		for (size_t i = 1; i < vlines.size(); ++i)
		{
			std::istringstream stm(vlines[i]);

			int index, A, B, C, mark;

			if (stm >> index >> A >> B >> C >> mark)
			{
				if (msetMap.find(mark) == msetMap.end())
				{
					msetMap.insert(pair<int,int>(mark, countSet++));

					// Extend faces vector
					vfaces.resize(countSet);
				}

				vfaces[msetMap[mark]].push_back(TetGenFace(C, B, A));
			}
		}
	}


	void TetGenReader::convertIndexing(vector<TetGenNode> &vn, vector<TetGenTetra> &vt, vector<TetGenEdge> &ve, vector<vector<TetGenFace> > &vf) {
		int shift = -vn[0].index;

		for (size_t i = 0; i < vn.size(); ++i) {
			vn[i].index += shift;
		}

		for (size_t i = 0; i < vt.size(); ++i) {
			vt[i].A += shift;
			vt[i].B += shift;
			vt[i].C += shift;
			vt[i].D += shift;
		}

		for (size_t i = 0; i < ve.size(); ++i) {
			ve[i].A += shift;
			ve[i].B += shift;
		}

		for (size_t i = 0; i < vf.size(); ++i) {
			for (size_t j = 0; j < vf[i].size(); ++j)
			{
				vf[i][j].A += shift;
				vf[i][j].B += shift;
				vf[i][j].C += shift;
			}
		}
	}

	vector<TetGenNode> TetGenReader::readNodes(string path) {
		vector<string> lines;
		readFileLines(path, lines);
		vector<TetGenNode> nodes;

		// skip first line - not relevant to nodes
		for (size_t i = 1; i < lines.size(); ++i) {
			std::istringstream stm(lines[i]);

			// perform formatted input from the string stream
			double x, y, z;
			int index;
			if (stm >> index >> x >> y >> z)
			{
				nodes.push_back(TetGenNode(index, x, y, z));
				//std::cout << index << " " << x << " " << y << " " << z << std::endl;
			}
		}

		return nodes;
	}

	vector<TetGenFace> TetGenReader::readFaces(string path) {
		vector<string> lines;
		readFileLines(path, lines);
		vector<TetGenFace> faces;

		// skip first line - not relevant to nodes
		for (size_t i = 1; i < lines.size(); ++i) {
			std::istringstream stm(lines[i]);

			// perform formatted input from the string stream
			int index;
			int A, B, C, temp;
			if (stm >> index >> A >> B >> C >> temp)
			{
				faces.push_back(TetGenFace(A, B, C));
				//std::cout << index << " " << A << " " << B << " " << C << std::endl;
			}
		}

		return faces;
	}


	vector<TetGenTetra> TetGenReader::readTetra(string path) {
		vector<string> lines;
		readFileLines(path, lines);
		vector<TetGenTetra> tets;

		// skip first line - not relevant to nodes
		for (size_t i = 1; i < lines.size(); ++i) {
			std::istringstream stm(lines[i]);

			// perform formatted input from the string stream
			int index;
			int A, B, C, D;
			if (stm >> index >> A >> B >> C >> D)
			{
				tets.push_back(TetGenTetra(A, B, C, D));
				//std::cout << index << " " << A << " " << B << " " << C << " " << D << std::endl;
			}
		}

		return tets;
	}

	vector<TetGenEdge> TetGenReader::readEdges(string path) {
		vector<string> lines;
		readFileLines(path, lines);
		vector<TetGenEdge> edges;

		// skip first line - not relevant to nodes
		for (size_t i = 1; i < lines.size(); ++i) {
			std::istringstream stm(lines[i]);

			// perform formatted input from the string stream
			int index;
			int A, B, temp;
			if (stm >> index >> A >> B >> temp)
			{
				edges.push_back(TetGenEdge(A, B));
				//std::cout << index << " " << A << " " << B << " " << C << std::endl;
			}
		}

		return edges;
	}

	// Deprecated

	//void TetGenReader::convertIndexingToStartFrom0(vector<TetGenNode> &n, vector<TetGenFace> &f, vector<TetGenTetra> &t, vector<TetGenEdge> &e) {
	//	int shift = -n[0].index;
	//	for (size_t i = 0; i < n.size(); ++i) {
	//		n[i].index += shift;
	//	}
	//	for (size_t i = 0; i < f.size(); ++i) {
	//		f[i].A += shift;
	//		f[i].B += shift;
	//		f[i].C += shift;
	//	}
	//	for (size_t i = 0; i < t.size(); ++i) {
	//		t[i].A += shift;
	//		t[i].B += shift;
	//		t[i].C += shift;
	//		t[i].D += shift;
	//	}
	//	for (size_t i = 0; i < e.size(); ++i) {
	//		e[i].A += shift;
	//		e[i].B += shift;
	//	}
	//}

	//void TetGenReader::convertIndexingToStartFrom0(vector<TetGenNode> &n, vector<TetGenFace> &f, vector<TetGenFace> &f2, vector<TetGenFace> &f3, vector<TetGenTetra> &t, vector<TetGenEdge> &e) {
	//	int shift = -n[0].index;
	//	for (size_t i = 0; i < n.size(); ++i) {
	//		n[i].index += shift;
	//	}
	//	for (size_t i = 0; i < f.size(); ++i) {
	//		f[i].A += shift;
	//		f[i].B += shift;
	//		f[i].C += shift;
	//	}
	//	for (size_t i = 0; i < f2.size(); ++i) {
	//		f2[i].A += shift;
	//		f2[i].B += shift;
	//		f2[i].C += shift;
	//	}
	//	for (size_t i = 0; i < f3.size(); ++i) {
	//		f3[i].A += shift;
	//		f3[i].B += shift;
	//		f3[i].C += shift;
	//	}
	//	for (size_t i = 0; i < t.size(); ++i) {
	//		t[i].A += shift;
	//		t[i].B += shift;
	//		t[i].C += shift;
	//		t[i].D += shift;
	//	}
	//	for (size_t i = 0; i < e.size(); ++i) {
	//		e[i].A += shift;
	//		e[i].B += shift;
	//	}
	//}

	//void TetGenReader::rearrangeIndices(vector<TetGenNode> &n, vector<TetGenFace> &f, vector<TetGenTetra> &t, vector<TetGenEdge> &e, vector<TetGenFace> &cf) {
	//	int shift = -n[0].index;
	//	int max = 0;
	//	for (size_t i = 0; i < n.size(); ++i) {
	//		n[i].index += shift;
	//	}
	//	for (size_t i = 0; i < f.size(); ++i) {
	//		f[i].A += shift;
	//		f[i].B += shift;
	//		f[i].C += shift;
	//	}
	//	for (size_t i = 0; i < t.size(); ++i) {
	//		t[i].A += shift;
	//		t[i].B += shift;
	//		t[i].C += shift;
	//		t[i].D += shift;
	//	}
	//	for (size_t i = 0; i < cf.size(); ++i) {
	//		cf[i].A += shift;
	//		cf[i].B += shift;
	//		cf[i].C += shift;
	//	}
	//	for (size_t i = 0; i < e.size(); ++i) {
	//		e[i].A += shift;
	//		e[i].B += shift;
	//	}
	//}

	//// checks if node mapping is unique (no a->b, a->c, or a->c, b->c)
	//bool TetGenReader::isNodeMappingCorrect(vector<std::pair<int, int>> pairs) {
	//	vector<int> a(pairs.size());
	//	vector<int> b(pairs.size());

	//	for (size_t i = 0; i < pairs.size(); ++i) {
	//		a[i] = pairs[i].first;
	//		b[i] = pairs[i].second;
	//	}

	//	std::sort(a.begin(), a.end());
	//	for (size_t i = 1; i < a.size(); ++i) {
	//		if (a[i] == a[i - 1]) {
	//			return false;
	//		}
	//	}

	//	std::sort(b.begin(), b.end());
	//	for (size_t i = 1; i < b.size(); ++i) {
	//		if (b[i] == b[i - 1]) {
	//			return false;
	//		}
	//	}

	//	return true;
	//}


	//void TetGenReader::readFacesPoly(string fileName, vector<TetGenFace>& faces, vector<TetGenFace>& internalFaces, vector<TetGenFace>& externaFaces) 
	//{
	//	vector<string> lines;
	//	readFileLinesInto(fileName, lines);

	//	// skip first line - not relevant to nodes
	//	for (size_t i = 1; i < lines.size(); ++i) 
	//	{
	//		std::istringstream stm(lines[i]);

	//		int index, A, B, C, mark;

	//		if (stm >> index >> A >> B >> C >> mark)
	//		{
	//			faces.push_back(TetGenFace(C, B, A));

	//			if (mark == 1) 
	//			{
	//				externaFaces.push_back(TetGenFace(C, B, A));
	//			}
	//			else 
	//			{
	//				internalFaces.push_back(TetGenFace(C, B, A));
	//			}
	//		}
	//	}
	//}

	//vector<int> TetGenReader::matchNodesUsingConvexHull(vector<TetGenNode> n1, vector<TetGenNode> n2) {
	//	double minX, maxX, minY, maxY, minZ, maxZ;

	//	minX = maxX = n2[0].x;
	//	minY = maxY = n2[0].y;
	//	minZ = maxZ = n2[0].z;
	//	for (int i = 1; i < n2.size(); ++i) {
	//		minX = std::min(n2[i].x, minX);
	//		maxX = std::max(n2[i].x, maxX);
	//		minY = std::min(n2[i].y, minY);
	//		maxY = std::max(n2[i].y, maxY);
	//		minZ = std::min(n2[i].z, minZ);
	//		maxZ = std::max(n2[i].z, maxZ);
	//	}

	//	vector<int> matches(n1.size());
	//	for (size_t i = 0; i < n1.size(); ++i) {
	//		bool matched = false;

	//		if (minX <= n1[i].x && n1[i].x <= maxX &&
	//			minY <= n1[i].y && n1[i].y <= maxY &&
	//			minZ <= n1[i].z && n1[i].z <= maxZ) {

	//			matches[n1[i].index] = 0;
	//			matched = true;
	//		}

	//		if (!matched) {
	//			matches[n1[i].index] = -1;
	//		}
	//	}

	//	return matches;
	//}

	//vector<int> TetGenReader::matchNodes(vector<TetGenNode> n1, vector<TetGenNode> n2) {
	//	vector<int> matches(n1.size());
	//	for (size_t i = 0; i < n1.size(); ++i) {
	//		bool matched = false;
	//		for (size_t j = 0; j < n2.size(); ++j) {
	//			if (n1[i].Equals(n2[j])) {
	//				if (matched) {
	//					// if we match more than once for the same node, there may be a problem!
	//					// instead of this after function runs, we can call isNodeMappingCorrect <- slower
	//					std::cout << "Matching error: Consider increase node equals precision" << std::endl;
	//				}

	//				matches[n1[i].index] = n2[j].index;
	//				matched = true;
	//			}
	//		}
	//		if (!matched) {
	//			matches[n1[i].index] = -1;
	//		}
	//	}

	//	return matches;
	//}

	//bool TetGenReader::readPoly(string file, MatrixXd& mNodes, MatrixXi& mTetra, MatrixXi& mEdges, MatrixXi& mIntFaces, MatrixXi& mExtFaces)
	//{
	//	vector<TetGenNode> vnodes = readNodes(file + ".1.node");
	//	vector<TetGenEdge> vedges = readEdges(file + ".1.edge");
	//	vector<TetGenTetra> vtetra = readTetra(file + ".1.ele");

	//	vector<TetGenFace> vfaces;
	//	vector<TetGenFace> vintFaces;
	//	vector<TetGenFace> vextFaces;
	//	readFacesPoly(file + ".1.face", vfaces, vintFaces, vextFaces);

	//	convertIndexingToStartFrom0(vnodes, vfaces, vintFaces, vextFaces, vtetra, vedges);

	//	int numNodes = (int)vnodes.size();
	//	int numEdges = (int)vedges.size();
	//	int numTetra = (int)vtetra.size();
	//	int numIntFaces = (int)vintFaces.size();
	//	int numExtFaces = (int)vextFaces.size();
	//	mNodes.resize(numNodes, 3);
	//	mEdges.resize(numEdges, 2);
	//	mTetra.resize(numTetra, 4);
	//	mIntFaces.resize(numIntFaces, 3);
	//	mExtFaces.resize(numExtFaces, 3);

	//	for (int i = 0; i < numNodes; ++i)
	//		mNodes.row(i) = vnodes[i].GetPosition();

	//	for (int i = 0; i < numEdges; ++i)
	//	{
	//		mEdges(i, 0) = vedges[i].A;
	//		mEdges(i, 1) = vedges[i].B;
	//	}
	//	for (int i = 0; i < numTetra; ++i)
	//	{
	//		mTetra(i, 0) = vtetra[i].A;
	//		mTetra(i, 1) = vtetra[i].B;
	//		mTetra(i, 2) = vtetra[i].C;
	//		mTetra(i, 3) = vtetra[i].D;
	//	}
	//	for (int i = 0; i < numIntFaces; ++i)
	//	{
	//		mIntFaces(i, 0) = vintFaces[i].A;
	//		mIntFaces(i, 1) = vintFaces[i].B;
	//		mIntFaces(i, 2) = vintFaces[i].C;
	//	}
	//	for (int i = 0; i < numExtFaces; ++i)
	//	{
	//		mExtFaces(i, 0) = vextFaces[i].A;
	//		mExtFaces(i, 1) = vextFaces[i].B;
	//		mExtFaces(i, 2) = vextFaces[i].C;
	//	}

	//	return true;
	//}


}

