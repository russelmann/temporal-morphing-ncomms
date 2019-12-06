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

#include "smp_addons.h"
#include "Eigen/Eigen"

#include "geogram/mesh/mesh.h"
#include "geogram/mesh/mesh_io.h"
#include "geogram/mesh/mesh_remesh.h"
#include "geogram/mesh/mesh_geometry.h"
#include "geogram/basic/command_line.h"
#include "geogram/basic/command_line_args.h"

namespace smpup
{
	void extract_geo_mesh(const GEO::Mesh& mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
	double get_mean_edgelen(const GEO::Mesh& mesh);
	void flip_geo_faces(GEO::Mesh& mesh);

	class GeoRemesher
	{
	private:
		static bool geogram_initialized;
		//GEO::Mesh mesh_imported;
		//GEO::Mesh mesh_flattened;
		//GEO::Mesh mesh_remeshed;

		Eigen::MatrixXd V_imported;
		Eigen::MatrixXd uv_imported;
		Eigen::VectorXd u_imported; // scaling exp factors
		Eigen::MatrixXi F_imported;

		Eigen::MatrixXd V_remeshed;
		Eigen::MatrixXd uv_remeshed;
		Eigen::MatrixXi F_remeshed;

		Eigen::Vector2f scaling_bounds;
		float max_stretch; // scaling_bounds(1) / scaling_bounds(0)

	public:
		GeoRemesher()
		{
			scaling_bounds.setZero();
			max_stretch = 0;
			if (!geogram_initialized)
			{
				GEO::initialize();
				//GEO::mesh_io_initialize();
				GEO::CmdLine::import_arg_group("algo");
				GEO::CmdLine::set_arg("algo:delaunay", "delaunay");
				geogram_initialized = true;
			}
		}

		inline bool has_imported_mesh() const { return 0 < V_imported.size(); }
		inline bool has_remeshed_mesh() const { return 0 < V_remeshed.size(); }

		bool load_mesh(std::string fname, bool recenter = true);
		bool flatten_imported_mesh(Eigen::VectorXd* u_test = nullptr);

		inline float* get_scaling_bounds() { return scaling_bounds.data(); }
		inline float* get_max_stretch() { return &max_stretch; }

		inline void extract_imported_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F)
		{
			V = V_imported;
			F = F_imported;
		}
		inline void extract_flattened_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F)
		{
			V.resize(V_imported.rows(), 3);
			V.leftCols(2) = uv_imported;
			V.rightCols(1).setZero();
			F = F_imported;
		}
		inline void extract_remeshed_mesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F)
		{
			V = V_remeshed;
			F = F_remeshed;
		}

		void remesh(int vertices, int lloyd_iter = 50, int newton_iter = 3000, int newton_m = 30);

		int vertices_by_target_edgelen(double edgelen);

		inline const Eigen::MatrixXd& get_V_imported() const { return V_imported; }
		inline const Eigen::MatrixXi& get_F_imported() const { return F_imported; }

		inline const Eigen::MatrixXd& get_V_remeshed() const { return V_remeshed; }
		inline const Eigen::MatrixXi& get_F_remeshed() const { return F_remeshed; }

		int get_remeshed_mean_edgelen() const;

		void flip_remeshed_faces()
		{
			F_remeshed *= Eigen::Vector3i{ 2, 1, 0 }.asPermutation();
		}

		void remove_isolated_trangles();

		bool write_remeshed(std::string fname);
	};

}
