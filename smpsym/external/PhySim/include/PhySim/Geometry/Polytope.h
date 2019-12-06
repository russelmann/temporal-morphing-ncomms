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

namespace PhySim
{

	using namespace std;
	using namespace Eigen;

	class Polytope : public IPolytope
	{
	protected:

		// Topo

		int m_order;
		vector<Node*> m_vnodes;
		vector<Edge*> m_vedges;
		vector<Face*> m_vfaces;
		vector<Cell*> m_vcells;
		vector<Edge*> m_vadjacentEdges;
		vector<Face*> m_vadjacentFaces;
		vector<Cell*> m_vadjacentCells;

		// Data

		DoFSet* m_dofs;

	public:

		inline vector<Node*>& Nodes() { return this->m_vnodes; };
		inline vector<Edge*>& Edges() { return this->m_vedges; };
		inline vector<Face*>& Faces() { return this->m_vfaces; };
		inline vector<Cell*>& Cells() { return this->m_vcells; };
		inline vector<Edge*>& AdjacentEdges() { return this->m_vadjacentEdges; };
		inline vector<Face*>& AdjacentFaces() { return this->m_vadjacentFaces; };
		inline vector<Cell*>& AdjacentCells() { return this->m_vadjacentCells; };

		virtual Real ComputeVolume(Space s = Space::MAT) const;
		virtual Vector3d ComputeCentroid(Space s = Space::MAT) const;
		virtual Matrix3d ComputeRotation(Space f = Space::MAT, 
										 Space t = Space::DEF) const;

		inline DoFSet* DoF() { return this->m_dofs; };

		virtual void GetNodeMatrix(MatrixXd& mN, Space s = Space::MAT) const override;
		virtual void SetNodeMatrix(const MatrixXd& mN, Space s = Space::MAT) override;

		virtual void GetQuadrature(int num, vector<VectorXd>& vQPos, vector<Real>& vQWei) const;

		/**
		* Computes the value of the shape function for each of the 
		* nodes of the polytope at the specified point. The point
		* must be specified in iso-parametric coordinates.
		*
		* @param vp The point where to evaluate the shape function in iso-parametric.
		* @param vN (Nx1)-vector The resulting value of the shape function at the point.
		*/
		virtual void ComputeShapeFunction(const VectorXd& vp, VectorXd& vN, Space s = Space::MAT) const override;

		/**
		* Computes the derivative of the shape function w.r.t. iso-parametric coordinates 
		* for each of the nodes of the polytope at the specified point. The point must be 
		* also specified in iso-parametric coordinates.
		*
		* @param vp The point where to evaluate the shape derivative in iso-parametric.
		* @param vN (NxO)-matrix The resulting value of the shape drivative at the point.
		*/
		virtual void ComputeShapeDerivative(const VectorXd& vp, MatrixXd& mB, Space s = Space::MAT) const override;

		/**
		* Computes the affine transformation to go from natural coordinates
		* in the specified space to coordinates of the iso-parametric polytope.
		* The transformation is provided in the shape i = b + A*n.
		*
		* @param vb The independent vector of the affine transformation.
		* @param mA The corresponding matrix of the affine transformation.
		* @param s The space from which the coordinates are transformed.
		*/
		virtual void ComputeNat2IsoTransform(VectorXd& vb, MatrixXd& mA, Space s = Space::MAT) const override;

		/**
		* Computes the affine transformation to go from coordinates
		* of the iso-parametric polytope to the natural coordinates in
		* the specified space. The transformation is provided in the 
		* shape n = b + A*i.
		*
		* @param vb The independent vector of the affine transformation.
		* @param mA The corresponding matrix of the affine transformation.
		* @param s The space to which the coordinates are transformed.
		*/
		virtual void ComputeIso2NatTransform(VectorXd& vb, MatrixXd& mA, Space s = Space::MAT) const override;


		/**
		* Interpolates the position of the nodes in the chosen space
		* at the specified point. The point must be in iso-parametric
		* coordinates.
		* 
		* @param vp The point where to interpolate node position.
		* @param s The space in which node position is interpolated.
		*/
		virtual Vector3d InterpolatePosition(const VectorXd& vp, Space s = Space::MAT) const override;

		/**
		* Compute the embedded node corresponding to the specified 
		* 3D position. The point must be in natural coordinates in 
		* the specified space.
		*
		* @param vp Coordinates of the point in natural coordinates.
		* @param s The space in which node position is interpolated.
		*/
		virtual NodeEBD ComputeEmbedding(const Vector3d& vp, Space s = Space::MAT) override;

		/**
		* Check if the specified iso-parametric coordinates are ok.
		*
		* @param vp The iso-parametric coordinates.
		*/
		virtual bool IsValidParametric(const VectorXd& vp) const override;

		/**
		* Set the positions of the subelement nodes according to the
		* particular discretization. For instance, 10-node quadratic
		* tetrahedron have 6 subelement nodes, one at each tetrahedron.
		*
		* @param s The space where positions are set.
		*/
		virtual void SetSubelementPositions(Space s = Space::DEF);

	};

}