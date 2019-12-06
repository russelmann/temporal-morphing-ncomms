//==========================================================
//
//	PhyEye library. Generic library for physics visualization.
//
//	Authors:
//			Jesus Perez Rodriguez, IST Austria
//
//==========================================================

#pragma once

#include <PhySim/CommonIncludes.h>
#include <PhySim/PhySimInterface.h>

namespace PhySim
{
	using namespace std;
	using namespace Eigen;

	struct TetGenNode;
	struct TetGenFace;
	struct TetGenEdge;
	struct TetGenTetra;

	static class TetGenReader
	{
	public:
		static bool readPoly(string file, MatrixXd& mNodes, MatrixXi& mTetra, MatrixXi& mEdges, vector<MatrixXi>& vmFaces);
		
	private:
		// Readers
		static vector<TetGenNode> readNodes(string fileName);
		static vector<TetGenFace> readFaces(string fileName);
		static vector<TetGenEdge> readEdges(string fileName);
		static vector<TetGenTetra> readTetra(string fileName);
		static void readFacesPoly(string fileName, vector<vector<TetGenFace> >& vfaces);
		static void convertIndexing(vector<TetGenNode> &vn, vector<TetGenTetra> &vt, vector<TetGenEdge> &ve, vector<vector<TetGenFace> > &vf);

		//// Normalize data
		//static vector<int> matchNodesUsingConvexHull(vector<TetGenNode> n1, vector<TetGenNode> n2);
		//static vector<int> matchNodes(vector<TetGenNode> n1, vector<TetGenNode> n2);
		//static void rearrangeIndices(vector<TetGenNode> &n, vector<TetGenFace> &f, vector<TetGenTetra> &t, vector<TetGenEdge> &e, vector<TetGenFace> &cf);
	
		//static void convertIndexingToStartFrom0(vector<TetGenNode> &n, vector<TetGenFace> &f, vector<TetGenTetra> &t, vector<TetGenEdge> &e);
		//static void convertIndexingToStartFrom0(vector<TetGenNode> &n, vector<TetGenFace> &f, vector<TetGenFace> &f2, vector<TetGenFace> &f3, vector<TetGenTetra> &t, vector<TetGenEdge> &e);
	

		//// Data correctness check
		//static bool isNodeMappingCorrect(vector<std::pair<int, int>> pairs);

		//static bool readPoly(string file, MatrixXd& mNodes, MatrixXi& mTetra, MatrixXi& mEdges, MatrixXi& mIntFaces, MatrixXi& mExtFaces);

	};

}
