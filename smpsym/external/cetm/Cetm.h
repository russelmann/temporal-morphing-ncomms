#ifndef CETM_H
#define CETM_H

#include "cetm/Mesh.h"
#include "Utils.h"
#include <stack>

namespace cetm
{

	class Cetm
	{
	public:
		// constructor
		Cetm(Mesh& mesh0, int optScheme0) :
			mesh(mesh0),
			thetas(mesh0.vertices.size()),
			lengths(mesh0.edges.size()),
			angles(mesh0.halfEdges.size()),
			OptScheme(optScheme0)
		{ }

		// parameterize
		void parameterize();

		// get scaling factors
		bool get_scaling_factors(Eigen::VectorXd& u) const;

		// computes quasi conformal error
		double computeQcError();


	protected:
		// sets target theta values
		void setTargetThetas();

		// computes energy, gradient and hessian
		void computeEnergy(double& energy, const Eigen::VectorXd& u);
		void computeGradient(Eigen::VectorXd& gradient, const Eigen::VectorXd& u);
		void computeHessian(Eigen::SparseMatrix<double>& hessian, const Eigen::VectorXd& u);

		// sets edge lengths
		void setEdgeLengthsAndAngles();

		// computes scale factors
		bool computeScaleFactors();

		// determines position of unfixed face vertex
		void performFaceLayout(HalfEdgeCIter he, const Eigen::Vector2d& dir,
			std::unordered_map<int, bool>& visited, std::stack<EdgeCIter>& stack);

		// sets uvs
		void setUVs();

		// normalize
		void normalize();

		// member variables
		std::unordered_map<int, int> index;
		Eigen::VectorXd thetas;
		Eigen::VectorXd lengths;
		Eigen::VectorXd angles;
		Solver solver;
		int OptScheme;
		Mesh& mesh;
	};

}

#endif
