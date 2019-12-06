// This file is part of smpsym, an inverse design and simulation tool for
// research paper Guseinov R. et al "Programming temporal morphing of
// self-actuated shells"
//
// Copyright (C) 2019 Ruslan Guseinov <guseynov.ruslan@gmail.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include "utils/umath.h"
#include "utils/smp_geogram.h"
#include "geogram/mesh/mesh_repair.h"

#include "igl/edge_lengths.h"
#include "igl/vertex_triangle_adjacency.h"
#include "igl/per_vertex_normals.h"
#include "igl/per_face_normals.h"
#include "igl/boundary_loop.h"
#include "igl/writeOBJ.h"
#include "igl/readOBJ.h"
#include "igl/triangle_triangle_adjacency.h"
#include "igl/remove_unreferenced.h"
#include "igl/slice.h"
#include "igl/point_mesh_squared_distance.h"
#include "igl/barycentric_to_global.h"

#include "cetm/Mesh.h"

namespace smpup
{

	// Create GEO mesh, 2D input points assume z == 0
	void create_geo_mesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, GEO::Mesh& mesh)
	{
		mesh.clear();
		mesh.vertices.create_vertices(V.rows());
		for (int i = 0; i < mesh.vertices.nb(); ++i)
		{
			mesh.vertices.point(i).x = V(i, 0);
			mesh.vertices.point(i).y = V(i, 1);
			if (V.cols() == 3) mesh.vertices.point(i).z = V(i, 2);
			else mesh.vertices.point(i).z = 0;
		}
		mesh.facets.create_facets(F.rows(), 3);
		for (int i = 0; i < mesh.facets.nb(); ++i)
		{
			mesh.facets.set_vertex(i, 0, F(i, 0));
			mesh.facets.set_vertex(i, 1, F(i, 1));
			mesh.facets.set_vertex(i, 2, F(i, 2));
		}
	}

	void extract_geo_mesh(const GEO::Mesh& mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
	{
		V.resize(mesh.vertices.nb(), 3);
		for (int i = 0; i < mesh.vertices.nb(); ++i)
		{
			V(i, 0) = mesh.vertices.point(i).x;
			V(i, 1) = mesh.vertices.point(i).y;
			V(i, 2) = mesh.vertices.point(i).z;
		}
		F.resize(mesh.facets.nb(), 3);
		for (int i = 0; i < mesh.facets.nb(); ++i)
		{
			F(i, 0) = mesh.facets.vertex(i, 0);
			F(i, 1) = mesh.facets.vertex(i, 1);
			F(i, 2) = mesh.facets.vertex(i, 2);
		}
	}

	bool GeoRemesher::geogram_initialized = false;

	bool GeoRemesher::load_mesh(std::string fname, bool recenter)
	{
		Eigen::MatrixXd TC, N;
		Eigen::MatrixXi FTC, FN;
		if (!igl::readOBJ(fname, V_imported, TC, N, F_imported, FTC, FN)) return false;
		V_imported.rowwise() -= V_imported.colwise().mean();

		if (TC.size() == 0) flatten_imported_mesh();
		else
		{
			if (FTC.size() != 0)
			{
				Eigen::VectorXi perm = Eigen::VectorXi::Constant(V_imported.rows(), -1);
				for (int i = 0; i < FTC.rows(); ++i)
					for (int j = 0; j < FTC.cols(); ++j)
					{
						if (perm(FTC(i, j)) == -1) perm(FTC(i, j)) = F_imported(i, j);
						else if (perm(FTC(i, j)) != F_imported(i, j)) return false; // wrong uv
					}
				if (perm.minCoeff() < 0) return false; // wrong uv
				TC = perm.asPermutation() * TC;
			}
			TC.rowwise() -= TC.colwise().minCoeff();
			TC.array().rowwise() /= TC.array().colwise().maxCoeff() / 2.0;
			TC.rowwise() -= Eigen::Vector2d::Ones().transpose();
			uv_imported = TC;
		}

		u_imported = compute_scaling_factors(V_imported, uv_imported, F_imported);
		auto scaling_factors = (2.0 * u_imported).array().exp();
		scaling_bounds(0) = scaling_factors.minCoeff();
		scaling_bounds(1) = scaling_factors.maxCoeff();
		max_stretch = scaling_bounds(1) / scaling_bounds(0);

		return true;
	}

	bool GeoRemesher::flatten_imported_mesh(Eigen::VectorXd* u_test)
	{
		if (V_imported.size() == 0) return false;

		cetm::Mesh mesh(V_imported, F_imported);
		mesh.parameterize();
		uv_imported.resize(V_imported.rows(), 2);
		for (int i = 0; i < uv_imported.rows(); ++i)
			uv_imported.row(i) = mesh.vertices[i].uv;

		return true;
	}

	void GeoRemesher::remesh(int vertices, int lloyd_iter, int newton_iter, int newton_m)
	{
		GEO::Mesh mesh_flattened;
		create_geo_mesh(uv_imported, F_imported, mesh_flattened);

		GEO::Attribute<double> weight(mesh_flattened.vertices.attributes(), "weight");
		for (int i = 0; i < mesh_flattened.vertices.nb(); ++i)
			weight[i] = exp(-2 * u_imported(i)); //TODO: maybe times 4 is better?

		GEO::CmdLine::import_arg_group("remesh");
		GEO::Mesh mesh_remeshed;
		try {
			GEO::Mesh mesh_flattened_copy;
			mesh_flattened_copy.copy(mesh_flattened); // a copy is needed since GEO::mesh_repair and GEO::remesh_smooth affect the input mesh
			GEO::mesh_repair(mesh_flattened_copy); //TODO: remesh fails if not repaired (weird...)
			GEO::remesh_smooth(mesh_flattened_copy, mesh_remeshed, vertices, 6, lloyd_iter, newton_iter, newton_m);
		}
		catch (const std::exception& e) {
			std::cerr << "Received an exception: " << e.what() << std::endl;
			return;
		}

		// Lift back; probably, not ideal implementation, use igl::AABB::find() for 2D triangle search
		extract_geo_mesh(mesh_remeshed, uv_remeshed, F_remeshed);
		Eigen::VectorXi I;
		Eigen::MatrixXd C;
		Eigen::VectorXd sqrD;
		Eigen::MatrixXd V_flat_imported(uv_imported.rows(), 3);
		V_flat_imported.leftCols(2) = uv_imported;
		V_flat_imported.rightCols(1).setZero();
		igl::point_mesh_squared_distance(uv_remeshed, V_flat_imported, F_imported, sqrD, I, C);
		Eigen::MatrixXd br;
		Eigen::MatrixXd Va, Vb, Vc;
		Eigen::MatrixXi F2;
		Eigen::Vector3i xyz{ 0, 1, 2 };
		igl::slice(F_imported, I, xyz, F2);
		igl::slice(V_flat_imported, F2.col(0), xyz, Va);
		igl::slice(V_flat_imported, F2.col(1), xyz, Vb);
		igl::slice(V_flat_imported, F2.col(2), xyz, Vc);
		igl::barycentric_coordinates(C, Va, Vb, Vc, br);
		uv_remeshed.conservativeResize(Eigen::NoChange, 2);

		br.col(0) = I.cast<double>();
		V_remeshed = igl::barycentric_to_global(V_imported, F_imported, br);
	}

	int GeoRemesher::vertices_by_target_edgelen(double edgelen)
	{
		if (!has_imported_mesh()) return 0;
		Eigen::MatrixXd V;
		Eigen::MatrixXi F;
		extract_imported_mesh(V, F);
		double area2 = 0;
		for (int i = 0; i < F.rows(); ++i)
		{
			Eigen::Vector3d a = (V.row(F(i, 1)) - V.row(F(i, 0))).head(3);
			Eigen::Vector3d b = (V.row(F(i, 2)) - V.row(F(i, 0))).head(3);
			area2 += a.cross(b).norm();
		}
		double face_vertex_ratio = 1.8;
		return int(4. / sqrt(3) * area2 / (edgelen * edgelen) / face_vertex_ratio / 2.25);
	}

	int GeoRemesher::get_remeshed_mean_edgelen() const
	{
		Eigen::MatrixXd L;
		igl::edge_lengths(V_remeshed, F_remeshed, L); // Counts shared edges twice
		return L.mean();
	}

	// TODO: reimplement without geogram
	void GeoRemesher::remove_isolated_trangles()
	{
		GEO::Mesh mesh_remeshed;
		create_geo_mesh(V_remeshed, F_remeshed, mesh_remeshed);
		GEO::vector<GEO::index_t> faces_to_delete(mesh_remeshed.facets.nb(), 0);
		for (int i = 0; i < mesh_remeshed.facets.nb(); ++i)
		{
			int brd = 0;
			for (int j = 0; j < 3; ++j)
				brd += (mesh_remeshed.facets.adjacent(i, j) == -1 ? 1 : 0);
			if (1 < brd) faces_to_delete[i] = 1;
		}
		mesh_remeshed.facets.delete_elements(faces_to_delete);
		extract_geo_mesh(mesh_remeshed, V_remeshed, F_remeshed);
	}

	bool GeoRemesher::write_remeshed(std::string fname)
	{
		Eigen::Matrix2d d;
		d.col(0) = uv_remeshed.row(F_remeshed(0, 1)) - uv_remeshed.row(F_remeshed(0, 0));
		d.col(1) = uv_remeshed.row(F_remeshed(0, 2)) - uv_remeshed.row(F_remeshed(0, 0));
		if (d.determinant() < 0)
			uv_remeshed.col(1) *= -1;

		//Eigen::MatrixXd CN, FN;
		//igl::writeOBJ(fname, V_remeshed, F_remeshed, CN, FN, uv_remeshed, F_remeshed);
		std::ofstream fout(fname);
		fout << V_remeshed.format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, " ", "\n", "v ", "", "", "\n"));
		fout << uv_remeshed.format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, " ", "\n", "vt ", "", "", "\n"));
		fout << (F_remeshed.array() + 1).format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, " ", "\n", "f ", "", "", "\n"));
		return true;
	}
}
