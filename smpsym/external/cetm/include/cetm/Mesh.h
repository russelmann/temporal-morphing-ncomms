#ifndef MESH_H
#define MESH_H

#include "cetm/Types.h"
#include "cetm/Vertex.h"
#include "cetm/Edge.h"
#include "cetm/Face.h"
#include "cetm/HalfEdge.h"

namespace cetm
{

	class Mesh {
	public:
		// default constructor
		Mesh() { }

		// copy constructor
		Mesh(const Mesh& mesh)
		{
			*this = mesh;
		}

		// build mesh from vertices and faces
		Mesh(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces);

		// read mesh from file
		bool read(const std::string& fileName);

		// write mesh to file
		bool write(const std::string& fileName) const;

		// computes CETM conformal parameterization
		double parameterize(int optScheme = NEWTON);

		// output scaing factors
		bool get_scaling_factors(Eigen::VectorXd& u) const;

		// delaunayize
		void delaunayize();

		// computes mean edge length
		double meanEdgeLength();

		// member variables
		std::vector<HalfEdge> halfEdges;
		std::vector<Vertex> vertices;
		std::vector<Edge> edges;
		std::vector<Face> faces;
		std::vector<HalfEdgeIter> boundaries;
		Eigen::VectorXd u;

	private:
		// center mesh about origin and rescale to unit radius
		void normalize();
	};

}
#endif
