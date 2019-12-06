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

#include "smp_mesh.h"
#include "utils/umath.h"

#include "cetm/Mesh.h"

#include "igl/boundary_facets.h"
#include "igl/unique.h"
#include "igl/per_vertex_normals.h"
#include "igl/triangle_triangle_adjacency.h"

namespace smpup
{
	void SmpMesh::compute_parameters(bool is_post_load)
	{
		igl::triangle_triangle_adjacency(F, ff);
		ff *= Eigen::PermutationMatrix<3>(Eigen::Vector3i{ 1, 2, 0 }); // to agree with matlab convention

		FF.resize((ff.size() - ff.array().cwiseEqual(-1).count()) / 2, 2);
		fi.resizeLike(FF);
		fi.setConstant(-1);
		int ind = 0;
		for (int i = 0; i < F.rows(); ++i)
			for (int j = 0; j < 3; ++j)
			{
				if (ff(i, j) < 0 || ff(i, j) < i) continue;
				FF.row(ind) = Eigen::Vector2i{ i, ff(i, j) };
				fi(ind, 0) = j;
				int k = 0;
				while (ff(ff(i, j), k) != i) ++k;
				fi(ind, 1) = k;
				++ind;
			}
		Ff.resizeLike(F);
		Ff.setConstant(-1);
		for (int i = 0; i < FF.rows(); ++i)
		{
			Ff(FF(i, 0), fi(i, 0)) = i;
			Ff(FF(i, 1), fi(i, 1)) = i;
		}

		igl::per_vertex_normals(V, F, vnorm);
		igl::per_face_normals(V, F, fnorm);

		//TODO: test
		theta.resize(FF.rows());
		for (int i = 0; i < FF.rows(); ++i)
		{
			int f0 = FF(i, 0);
			int f1 = FF(i, 1);
			theta(i) = angle_vectors(fnorm.row(f0), fnorm.row(f1));
			if (theta(i) == 0.) continue;
			Eigen::Vector3d dir = V.row(F(f0, prev3(fi(i, 0)))) - V.row(F(f0, next3(fi(i, 0))));
			if (fnorm.row(f0).cross(fnorm.row(f1)).dot(dir) < 0)
				theta(i) *= -1;
		}

		Eigen::MatrixXi bf;
		igl::boundary_facets(F, bf);
		Eigen::VectorXi bv;
		igl::unique(bf, bv);
		Vb.resize(V.rows());
		Vb.setZero();
		for (int i = 0; i < bv.size(); ++i)
			Vb(bv(i)) = 1;
		FFb.resize(FF.rows());
		for (int i = 0; i < FFb.size(); ++i)
			FFb(i) = Vb(F(FF(i, 0), prev3(fi(i, 0)))) == 1 || Vb(F(FF(i, 0), next3(fi(i, 0)))) == 1;

		if (!is_post_load)
		{
			if (uv.size() == 0) compute_uv();
			uv.rowwise() -= uv.colwise().mean();

			eu = compute_scaling_factors(V, uv, F);
			eu = (2.0 * eu).array().exp();

			// Center 3D mesh
			V.rowwise() -= V.colwise().sum() / V.rows();

			update_parameters();

			// Scale 3D mesh
			scale3d.setZero(FF.rows());
			double rescale_V = 0;
			for (int i = 0; i < F.rows(); ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					if (ff(i, j) < 0) continue;
					double ang = angle_vectors(fnorm.row(i), fnorm.row(ff(i, j))) / 2;

					double eFc3dx = options.base_r + (options.base_th + options.base_p) * tan(ang) + options.min_bumper;
					double eFc3d = (Fc.row(i) - Fc.row(ff(i, j))).norm() / cos(ang) / 2;

					rescale_V = std::max(rescale_V, eFc3dx / eFc3d);
					scale3d(Ff(i, j)) = eFc3dx / eFc3d;
				}
			}
			V *= rescale_V;

			scale3d /= scale3d.minCoeff();

			update_parameters();

			// Scale 2D mesh
			double rescale_uv = 0;
			for (int i = 0; i < F.rows(); ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					if (ff(i, j) < 0) continue;
					double ang = angle_vectors(fnorm.row(i), fnorm.row(ff(i, j))) / 2;

					double eFc2dx = (Fc.row(i) - Fc.row(ff(i, j))).norm() / cos(ang) / 2 + (options.base_th + options.base_p) * tan(ang) + options.min_gap / 2;
					double eFc2d = (Fc2d.row(i) - Fc2d.row(ff(i, j))).norm() / 2;
					rescale_uv = std::max(rescale_uv, eFc2dx / eFc2d);
				}
			}
			uv *= rescale_uv;

			eu /= sqrt(eu(F(0, 0)) * eu(F(0, 1)));
			eu *= (uv.row(F(0, 0)) - uv.row(F(0, 1))).norm() / (V.row(F(0, 0)) - V.row(F(0, 1))).norm();
		}

		update_parameters();
		scale3d.maxCoeff(&act_max_scale);
	}

	void SmpMesh::compute_uv()
	{
		if (V.rows() < 0 || F.rows() < 0)
		{
			uv.resize(0, 2);
			return;
		}

		cetm::Mesh mesh(V, F);
		mesh.parameterize();
		uv.resize(mesh.vertices.size(), 2);
		int i = 0;
		for (cetm::VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); ++v)
		{
			uv.row(i) = v->uv;
			++i;
		}

		// make sure the flat stencil is not inverted
		//Eigen::Matrix2d hd;
		//hd.col(0) = (uv.row(F(0, 1)) - uv.row(F(0, 0))).transpose();
		//hd.col(1) = (uv.row(F(0, 2)) - uv.row(F(0, 0))).transpose();
		//if (hd.determinant() < 0) uv *= -1;
		// ^^^ Commented since it should never be the case
	}

	void SmpMesh::update_parameters()
	{
		int m = F.rows();

		Fc.resize(m, 3);
		Fco.resizeLike(Fc);
		Fcu.resizeLike(Fc);
		if (uv.size() > 0)
			Fc2d.resize(m, 2);

		for (int i = 0; i < F.rows(); ++i)
		{
			fnorm.row(i) = (V.row(F(i, 1)) - V.row(F(i, 0))).cross(V.row(F(i, 2)) - V.row(F(i, 0))).normalized();
			Fc.row(i) = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3.;
			Fco.row(i) = Fc.row(i) + fnorm.row(i) * (options.base_th / 2 + options.base_p);
			Fcu.row(i) = Fc.row(i) - fnorm.row(i) * (options.base_th / 2 + options.base_p);
			if (uv.size() > 0)
				Fc2d.row(i) = (uv.row(F(i, 0)) + uv.row(F(i, 1)) + uv.row(F(i, 2))) / 3.;
		}
	}

	void SmpMesh::init_single()
	{
		V.resize(4, 3);
		V.row(0) = Eigen::Vector3d{ 0,  1, 0 };
		V.row(1) = Eigen::Vector3d{ 0, -1, 0 };
		V.row(2) = Eigen::Vector3d{ 1,  0, 0 };
		V.row(3) = Eigen::Vector3d{ -1,  0, 0 };
		F.resize(2, 3);
		F.row(0) = Eigen::Vector3i{ 2, 0, 1 };
		F.row(1) = Eigen::Vector3i{ 3, 1, 0 };
		compute_parameters();
	}

	void SmpMesh::init_regular(int m, int n)
	{
		V.resize((3 + m * 2) * n + 1 + m, 3);
		V.rightCols(1).setZero();
		int ind = 0;
		for (int i = 0; i < 2 * n + 1; ++i)
			for (int j = (i + 1) % 2; j < m + 2; ++j)
			{
				V(ind, 0) = j + (i % 2) * 0.5;
				V(ind, 1) = i * sqrt(3) / 2;
				++ind;
			}
		F.resize((2 + m * 4) * n, 3);
		ind = 0;
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m + 1; ++j)
			{
				F(ind, 0) = 0;
				F(ind, 1) = m + 2;
				F(ind, 2) = m + 1;
				F.array().row(ind) += j + (3 + m * 2) * i;
				++ind;
				F(ind, 0) = m + 1;
				F(ind, 1) = m + 2;
				F(ind, 2) = 3 + m * 2;
				F.array().row(ind) += j + (3 + m * 2) * i;
				++ind;
				if (j < m)
				{
					F(ind, 0) = 0;
					F(ind, 1) = 1;
					F(ind, 2) = m + 2;
					F.array().row(ind) += j + (3 + m * 2) * i;
					++ind;
					F(ind, 0) = m + 2;
					F(ind, 1) = 4 + m * 2;
					F(ind, 2) = 3 + m * 2;
					F.array().row(ind) += j + (3 + m * 2) * i;
					++ind;
				}
			}
		}
		compute_parameters();

		Eigen::Vector2d dir = uv.row(1) - uv.row(0);
		double rot = atan2(dir(1), dir(0));
		for (int i = 0; i < uv.rows(); ++i)
			uv.row(i) *= Eigen::Rotation2Dd(rot).matrix();

		update_parameters();
	}

}
