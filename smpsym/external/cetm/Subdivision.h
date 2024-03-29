#ifndef SUBDIVISION_H
#define SUBDIVISION_H

#include "cetm/Mesh.h"
#include "cetm/MeshIO.h"

namespace cetm {

	class Subdivision {
	public:
		// implements loop subdivision
		void subdivide(Mesh& mesh) const;

	protected:
		// updates current vertices
		Eigen::Vector3d updateVertexPosition(const VertexCIter& v) const;

		// creates new coordinates
		Eigen::Vector3d createNewCoordinate(const HalfEdgeCIter& h) const;

		// creates new faces
		void createNewFaces(MeshData& data, const std::vector<int>& indices) const;
	};
}

#endif 
