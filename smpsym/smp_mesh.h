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

#ifndef SMP_MESH_H
#define SMP_MESH_H

#include "smp_addons.h"
#include "interface.h"
//#include "utils/self_serializer.h"

namespace smpup
{
	class SmpMesh
	{
	private:
		// Compute all parameters
		void compute_parameters(bool is_post_load = false);

		// Compute uv map
		void compute_uv();

		// Update some internal parameters
		void update_parameters();

	public:
		SmpOptions options;

		Eigen::MatrixX3d V;    // #V by 3: actuated stencil mesh vertices
		Eigen::MatrixX3i F;    // #F by 3: actuated stencil mesh faces

		// Indexing

		Eigen::MatrixXi  FF;   // #FF by 2: pairs of neighboring faces
		Eigen::MatrixXi  ff;   // #F by 3: for each face, indices of neighbor faces
		Eigen::MatrixXi  fi;   // #FF by 2: for each face pair, internal (0, 1, 2) indices of corresponding neighbor from ff
		Eigen::MatrixXd  Ff;   // #F by 3: for each face, indices of neighbor FF pair

		// Scale-independent 

		Eigen::MatrixX3d vnorm;    // #V by 3: vertex normals
		Eigen::MatrixX3d fnorm;    // #F by 3: face normals
		Eigen::VectorXd  theta;    // #FF : signed dihedral angles between neighboring faces
		Eigen::VectorXi  Vb;       // #V : 1 if boundary vertex, 0 if internal vertex
		Eigen::VectorXi  FFb;       // #FF : 1 if any of two shared vertices is boundary, 0 otherwise

		// Scale-dependent

		Eigen::MatrixX3d Fc;       // #F by 3: face centers for actuated stencil
		Eigen::MatrixX3d Fco;      // #F by 3: face centers for actuated stencil, outer side
		Eigen::MatrixX3d Fcu;      // #F by 3: face centers for actuated stencil, inner side

		// - Flat stentil

		Eigen::MatrixX2d uv;       // #V by 2: flat stencil mesh vertices
		Eigen::VectorXd  eu;       // #V conformal scaling factors
		Eigen::MatrixX2d Fc2d;     // #F by 2: face centers for flat stencil

								   // Analysis
		Eigen::VectorXd scale3d;   // #FF: scaling requited by actuated state
		int act_max_scale;         // actuator with maximal scaling

		void init(const Eigen::MatrixXd& V_in, const Eigen::MatrixXi& F_in, const Eigen::MatrixXd const* uv_in = nullptr)
		{
			V = V_in;
			F = F_in;
			if (uv_in) uv = *uv_in;
			compute_parameters();
		}

		void init_single();

		void init_regular(int m, int n);

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(cereal::make_nvp("options", options));
			ar(cereal::make_nvp("vertices", V));
			ar(cereal::make_nvp("faces", F));
			ar(cereal::make_nvp("uv", uv));
			ar(cereal::make_nvp("eu", eu));
			ar(cereal::make_nvp("scale3d", scale3d));

			if (Archive::is_loading::value) compute_parameters(true);
		}
	};

}

#endif
