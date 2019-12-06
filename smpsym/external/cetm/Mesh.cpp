#include "cetm/Mesh.h"
#include "cetm/MeshIO.h"
#include "Cetm.h"

using namespace cetm;

Mesh::Mesh(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces)
{
	MeshData mdata;
	for (int i = 0; i < vertices.rows(); ++i)
		mdata.positions.push_back(vertices.row(i));
	for (int i = 0; i < faces.rows(); ++i)
	{
		std::vector<Index> indf;
		indf.emplace_back(Index(faces(i, 0), -2, -2));
		indf.emplace_back(Index(faces(i, 1), -2, -2));
		indf.emplace_back(Index(faces(i, 2), -2, -2));
		mdata.indices.push_back(indf);
	}

	MeshIO::buildMesh(mdata, *this);
}

bool Mesh::read(const std::string& fileName)
{
    std::ifstream in(fileName.c_str());

    if (!in.is_open()) {
        std::cerr << "Error: Could not open file for reading" << std::endl;
        return false;
    }
    
    bool readSuccessful = false;
    if ((readSuccessful = MeshIO::read(in, *this))) {
        normalize();
    }
    
    in.close();
    return readSuccessful;
}

bool Mesh::write(const std::string& fileName) const
{
    std::ofstream out(fileName.c_str());
    
    if (!out.is_open()) {
        std::cerr << "Error: Could not open file for writing" << std::endl;
        return false;
    }
    
    MeshIO::write(out, *this);
    
    out.close();
    return false;
}

double Mesh::parameterize(int optScheme)
{
	Cetm param(*this, optScheme);
    param.parameterize();
    double qcError = param.computeQcError();
	param.get_scaling_factors(u);
    
    return qcError;
}

bool Mesh::get_scaling_factors(Eigen::VectorXd& u) const
{
	u = this->u;
	return true;
}

void Mesh::delaunayize()
{
    std::stack<EdgeIter> stack;
    std::unordered_map<int, bool> marked;
    
    // initialize
    for (EdgeIter e = edges.begin(); e != edges.end(); e++) {
        stack.push(e);
        marked[e->index] = true;
    }
    
    // delaunayize
    while (!stack.empty()) {
        EdgeIter e = stack.top(); stack.pop();
        marked[e->index] = false;
        
        // check delaunay property
        if (!e->isBoundary() && e->cotanWeigth() < 0.0) {
            e->flip();
            
            HalfEdgeCIter he = e->he;
            std::vector<EdgeIter> quadEdges = {he->next->edge,
                                               he->next->next->edge,
                                               he->flip->next->edge,
                                               he->flip->next->next->edge};
            // Push unmarked quad edges
            for (size_t i = 0; i < quadEdges.size(); i++) {
                EdgeIter qe = quadEdges[i];
                if (!marked[qe->index]) {
                    marked[qe->index] = true;
                    stack.push(qe);
                }
            }
        }
    }
}

double Mesh::meanEdgeLength()
{
    double sum = 0.0;
    for (EdgeCIter e = edges.begin(); e != edges.end(); e++) {
        sum += e->length();
    }
    
    return sum / edges.size();
}

void Mesh::normalize()
{
    // compute center of mass
    Eigen::Vector3d cm = Eigen::Vector3d::Zero();
    for (VertexCIter v = vertices.begin(); v != vertices.end(); v++) {
        cm += v->position;
    }
    cm /= (double)vertices.size();
    
    // translate to origin and determine radius
    double rMax = 0;
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->position -= cm;
        rMax = std::max(rMax, v->position.norm());
    }
    
    // rescale to unit sphere
    for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
        v->position /= rMax;
    }
}
