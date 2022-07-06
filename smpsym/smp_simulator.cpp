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

#include "smp_simulator.h"

#include <igl/per_vertex_normals.h>
#include "igl/exterior_edges.h"
#include "igl/triangle/triangulate.h"
#include "igl/unique.h"
#include "igl/arap.h"

#include "cetm/Mesh.h"

#include "utils/umath.h"

namespace smpup
{

	// Segment spring representation
	struct Spring
	{
		int edge;
		double rest_length;
		double thickness;
		double bumper;
		Eigen::VectorXd poly;
		Spring(int edge, double rest_length, double thickness, double bumper)
			: edge(edge), rest_length(rest_length), thickness(thickness), bumper(bumper)
		{
			poly.setZero(9);
		}
	};


	// Prism base representation for node indexing
	struct BasePrism
	{
		std::vector<int> index; // structure (16 nodes):  0: origin; 1-12: base nodes; 13-15: mid-sideface nodes

		BasePrism() : index(get_num_core_nodes()) { }

		inline void set_orig_index(int ind) { index[0] = ind; }
		inline int get_orig_index() { return index[0]; }

		void set_base_node_index(int ind, int k, int j, int s)
		{
			int ix = 1 + k + j * 2 + s * 6;
			index[ix] = ind;
		}

		inline int get_base_node_index(int k, int j, int s)
		{
			return index[1 + k + j * 2 + s * 6];
		}

		// Set mid-sideface nodes for segment membrane
		inline void set_mid_node_index(int ind, int j) { index[13 + j] = ind; }

		// Get mid-sideface nodes for segment membrane
		inline int get_mid_node_index(int j) { return index[13 + j]; }

		inline int get_num_nodes() { return index.size(); }
		inline static int get_num_core_nodes() { return 13 + 3; }
	};


	// Prism visualization mesh data
	struct BodyPrism
	{
		Eigen::MatrixXd V;
		BodyPrism()
		{
			V.resize(get_num_vertices(), 3);
		}

		void complete_mid_vertices(double indent_ratio, bool heightfield = false)
		{
			double alpha = (indent_ratio + 1) / 2;
			V.block( 0 + 24, 0, 6, 3) = alpha * V.block(0, 0, 6, 3) + (1 - alpha) * V.block(6, 0, 6, 3);
			V.block( 6 + 24, 0, 6, 3) = (1 - alpha) * V.block(0, 0, 6, 3) + alpha * V.block(6, 0, 6, 3);
			V.block(12 + 24, 0, 6, 3) = alpha * V.block(12, 0, 6, 3) + (1 - alpha) * V.block(6 + 12, 0, 6, 3);
			V.block(18 + 24, 0, 6, 3) = (1 - alpha) * V.block(12, 0, 6, 3) + alpha * V.block(6 + 12, 0, 6, 3);

			for (int i = 0; i < 12; ++i)
			{
				Eigen::Vector3d inner = V.row(i + 36) - V.row(i + 24);
				double inner_len = inner.norm();
				double outer_len = (V.row(i + 12) - V.row(i +  0)).norm();
				if (outer_len < inner_len)
				{
					inner *= outer_len / inner_len;
					V.row(i + 36) = V.row(i + 24) + inner.transpose();
				}
			}
		}

		static int get_num_vertices() { return (12 + 12) * 2; }
		static int get_num_faces() { return 12 + 8 + 12 + 12 + 36; }

		static int build_faces(Eigen::ArrayXXi& F)
		{
			F.resize(get_num_faces(), 4);

			// body side faces
			int ind = 0;
			for (int i = 0; i < 6; ++i)
			{
				F.row(ind++) = Eigen::Vector4i{ i, i + 6, (i + 1) % 6 + 6, 0 } + (i % 2 == 0 ? Eigen::Vector4i{ 24, 24, 24, 0} : Eigen::Vector4i::Zero());
				F.row(ind++) = Eigen::Vector4i{ (i + 1) % 6, i, (i + 1) % 6 + 6, 0 } + (i % 2 == 0 ? Eigen::Vector4i{ 24, 24, 24, 0 } : Eigen::Vector4i::Zero());
			}

			// body up and down faces
			F.row(ind++) = Eigen::Vector4i{ 0, 1, 2, 0 };
			F.row(ind++) = Eigen::Vector4i{ 0, 2, 3, 0 };
			F.row(ind++) = Eigen::Vector4i{ 0, 3, 5, 0 };
			F.row(ind++) = Eigen::Vector4i{ 3, 4, 5, 0 };

			F.row(ind++) = Eigen::Vector4i{ 7,  6,  8, 0 };
			F.row(ind++) = Eigen::Vector4i{ 8,  6,  9, 0 };
			F.row(ind++) = Eigen::Vector4i{ 9,  6, 11, 0 };
			F.row(ind++) = Eigen::Vector4i{ 10, 9, 11, 0 };

			// bumpers faces
			F.row(ind++) = Eigen::Vector4i{ 1,  0, 12, 0 };
			F.row(ind++) = Eigen::Vector4i{ 1, 12, 13, 0 };
			F.row(ind++) = Eigen::Vector4i{ 3,  2, 14, 0 };
			F.row(ind++) = Eigen::Vector4i{ 3, 14, 15, 0 };
			F.row(ind++) = Eigen::Vector4i{ 5,  4, 16, 0 };
			F.row(ind++) = Eigen::Vector4i{ 5, 16, 17, 0 };

			F.row(ind++) = Eigen::Vector4i{  6,  7, 18, 0 };
			F.row(ind++) = Eigen::Vector4i{ 18,  7, 19, 0 };
			F.row(ind++) = Eigen::Vector4i{  8,  9, 20, 0 };
			F.row(ind++) = Eigen::Vector4i{ 20,  9, 21, 0 };
			F.row(ind++) = Eigen::Vector4i{ 10, 11, 22, 0 };
			F.row(ind++) = Eigen::Vector4i{ 22, 11, 23, 0 };

			F.row(ind++) = Eigen::Vector4i{  0, 1, 12, 0 } + Eigen::Vector4i{ 24, 24, 24, 0};
			F.row(ind++) = Eigen::Vector4i{ 12, 1, 13, 0 } + Eigen::Vector4i{ 24, 24, 24, 0};
			F.row(ind++) = Eigen::Vector4i{  2, 3, 14, 0 } + Eigen::Vector4i{ 24, 24, 24, 0};
			F.row(ind++) = Eigen::Vector4i{ 14, 3, 15, 0 } + Eigen::Vector4i{ 24, 24, 24, 0};
			F.row(ind++) = Eigen::Vector4i{  4, 5, 16, 0 } + Eigen::Vector4i{ 24, 24, 24, 0};
			F.row(ind++) = Eigen::Vector4i{ 16, 5, 17, 0 } + Eigen::Vector4i{ 24, 24, 24, 0};

			F.row(ind++) = Eigen::Vector4i{  7,  6, 18, 0 } + Eigen::Vector4i{ 24, 24, 24, 0 };
			F.row(ind++) = Eigen::Vector4i{  7, 18, 19, 0 } + Eigen::Vector4i{ 24, 24, 24, 0 };
			F.row(ind++) = Eigen::Vector4i{  9,  8, 20, 0 } + Eigen::Vector4i{ 24, 24, 24, 0 };
			F.row(ind++) = Eigen::Vector4i{  9, 20, 21, 0 } + Eigen::Vector4i{ 24, 24, 24, 0 };
			F.row(ind++) = Eigen::Vector4i{ 11, 10, 22, 0 } + Eigen::Vector4i{ 24, 24, 24, 0 };
			F.row(ind++) = Eigen::Vector4i{ 11, 22, 23, 0 } + Eigen::Vector4i{ 24, 24, 24, 0 };

			for (int j = 0; j < 3; ++j)
				for (int s = 0; s < 2; ++s)
				{
					for (int k = 0; k < 2; ++k)
					{
						int vt = k + j * 2 + s * 6;
						int inv = (s == k) ? 0 : 1;
						F.row(ind++) = Eigen::Vector4i{ vt + 0 + 24 * inv, vt + 0 + 24 * (1 - inv), vt + 12, 0 };
						F.row(ind++) = Eigen::Vector4i{ vt + 12 + 12 * inv, vt + 12 + 12 * (1 - inv), vt + 36, 0 };
					}
					int vt = j * 2 + s * 6;
					F.row(ind++) = Eigen::Vector4i{ vt + 12, vt + 36 + s, vt + 37 - s, 0 };
					F.row(ind++) = Eigen::Vector4i{ vt + 12 + s, vt + 37, vt + 13 - s, 0 };
				}

			return 0;
		}
	};

	int SmpSimulator::load_smp(const Smp& smpd)
	{
		StaticSimulatorStruct::MEM_MODEL mode = settings.mem_model;

		smp = std::make_unique<Smp>(smpd);
		//smp = new Smp(*smpd);

		int m = smp->F.rows();

		StaticSimulatorStruct sss;

		sss.set_smp_options(smp->options);
		sss.set_settings(settings.simulator_settings);
		sss.set_linear_mem_model(smp_physim->get_linear_mem_model());
		sss.set_bracket_model(smp_physim->get_bracket_model());

		Eigen::MatrixXd nodes_flt(m * BasePrism::get_num_core_nodes(), 3);
		Eigen::MatrixXd nodes_act(m * BasePrism::get_num_core_nodes(), 3);

		std::vector<BasePrism> bprisms(m);
		std::vector<BodyPrism> body_prisms(m);

		// Initialize origins

		nodes_flt.topLeftCorner(m, 2) = smp->Fc2d;
		nodes_flt.topRightCorner(m, 1).setZero();
		nodes_act.topRows(m) = smp->Fc;
		for (int i = 0; i < m; ++i)
			bprisms[i].set_orig_index(i);

		// Initialize prism nodes

		Eigen::Vector3d zvec;
		zvec << 0, 0, 1;
		Eigen::Vector3d sideB;
		sideB(2) = 0;

		int indx = m;
		for (int i = 0; i < m; ++i)
		{
			int dr = 0;
			while (smp->ff(i, dr) < 0) ++dr;
			sideB.head(2) = smp->Fc2d.row(smp->ff(i, dr)) - smp->Fc2d.row(i);
			Eigen::Matrix3d R = vectors2rotation(smp->fnorm.row(i), smp->Fc.row(smp->ff(i, dr)) - smp->Fc.row(i), zvec, sideB);

			for (int s = 0; s < 2; ++s)
			{
				for (int j = 0; j < 3; ++j)
				{
					for (int k = 0; k < 2; ++k)
					{
						Eigen::Vector3d dir;
						if (smp->ff(i, j) < 0)
							dir = ( smp->V.row(smp->F(i, prev3(j))) + smp->V.row(smp->F(i, next3(j))) ) / 2;
						else
							dir = smp->Fc.row(smp->ff(i, j));
						dir -= smp->Fc.row(i);
						Eigen::Vector3d prp = dir.cross(smp->fnorm.row(i)).normalized();
						dir = smp->fnorm.row(i).cross(prp).normalized();
						double alpha = asin(smp->options.base_th / smp->options.base_r / 2) + 7 * M_PI / 180;
						Eigen::Vector3d vert = (dir * cos(alpha) + prp * (1 - k * 2) * sin(alpha)) * smp->options.base_r;
						vert += (1 - 2 * s) * smp->fnorm.row(i) * (smp->options.base_th + smp->options.base_p);
						nodes_act.row(indx) = vert + smp->Fc.row(i).transpose();

						vert = R.transpose() * vert;
						vert.head(2) += smp->Fc2d.row(i);
						nodes_flt.row(indx) = vert;

						bprisms[i].set_base_node_index(indx++, k, j, s);

						body_prisms[i].V.row(k + j * 2 + s * 6) = vert;
						Eigen::Vector3d vbump = R.transpose() * dir;

						if (smp->ff(i, j) >= 0)
						{
							double bump = ( smp->Fc.row(smp->ff(i, j)) - smp->Fc.row(i) ).norm() / 2;
							double ang = acos(smp->fnorm.row(i).dot(smp->fnorm.row(smp->ff(i, j)))) / 2;
							bump /= cos(ang);
							double sig = 1 - 2 * s;
							if (smp->fnorm.row(i).cross(smp->fnorm.row(smp->ff(i, j))).dot(smp->V.row(smp->F(i, prev3(j))) - smp->V.row(smp->F(i, next3(j)))) < 0)
								sig *= -1;
							bump += tan(ang) * (smp->options.base_th + smp->options.base_p) * sig;
							bump -= smp->options.base_r * cos(alpha);
							vbump *= bump;
						}
						else
						{
							//vbump *= 1e-6;
							vbump *= smp->options.min_bumper;
						}
						body_prisms[i].V.row(k + j * 2 + s * 6 + 12) = vert + vbump * 0.99; //N.B. trick to avoid flickering at rendering, magic number
					}
				}
			}

			for (int j = 0; j < 3; ++j)
			{
				nodes_act.row(indx).setZero();
				for (int k = 0; k < 2; ++k)
					for (int s = 0; s < 2; ++s)
						nodes_act.row(indx) += nodes_act.row(bprisms[i].get_base_node_index(k, j, s));
				nodes_act.row(indx) /= 4;

				nodes_flt.row(indx).setZero();
				for (int k = 0; k < 2; ++k)
					for (int s = 0; s < 2; ++s)
						nodes_flt.row(indx) += nodes_flt.row(bprisms[i].get_base_node_index(k, j, s));
				nodes_flt.row(indx) /= 4;

				bprisms[i].set_mid_node_index(indx++, j);
			}

			body_prisms[i].complete_mid_vertices(smp->options.base_p / smp->options.base_th, true);
		}

		std::vector<Eigen::Vector2i> vec_edges;
		std::vector<Spring> vec_springs;
		std::vector<Eigen::Vector2i> vec_scissors;
		std::vector<Eigen::Vector4i> vec_flips;
		std::vector<double> vec_flips_d;

		actuators.resize(smp->FF.rows());
		for (int i = 0; i < smp->FF.rows(); ++i)
			actuators[i] = new Actuator();

		//TODO: more realistic evaluation of number of springs/scissors
		vec_edges.reserve(smp->FF.rows() * 17);
		vec_springs.reserve(smp->FF.rows() * 5);
		vec_scissors.reserve(smp->FF.rows() * 12);
		vec_flips.reserve(12 * m);
		vec_flips_d.reserve(vec_flips.capacity());

		// Generate membrane springs

		if (mode == StaticSimulatorStruct::MEM_CONSTANT || mode == StaticSimulatorStruct::MEM_LINEAR)
		{
			sss.mspring_edges.resize(smp->FF.rows());
			for (int i = 0; i < smp->FF.rows(); ++i)
			{
				int f0 = smp->FF(i, 0);
				int f1 = smp->FF(i, 1);

				//Eigen::Vector2i edg(bprisms[f0].get_orig_index(), bprisms[f1].get_orig_index());
				Eigen::Vector2i edg(bprisms[f0].get_mid_node_index(smp->fi(i, 0)), bprisms[f1].get_mid_node_index(smp->fi(i, 1)));

				double rest_len_memsp = 0;
				vec_springs.emplace_back(vec_edges.size(), rest_len_memsp, 0, -1);
				
				actuators[i]->membrane_edge = vec_edges.size();
				sss.mspring_edges(i) = vec_edges.size();
				vec_edges.emplace_back(edg);
				
				double strength = 1;
				if (settings.mimic_bounary_weakening && smp->FFb(i))
					strength = settings.boundary_weakening;

				double spring_len = (nodes_flt.row(vec_edges[actuators[i]->membrane_edge](1)) - nodes_flt.row(vec_edges[actuators[i]->membrane_edge](0))).norm();
				vec_springs.back().poly.setZero(9);
				if (mode == StaticSimulatorStruct::MEM_CONSTANT)
				{
					vec_springs.back().poly(1) = -settings.constant_force;
				}
				else
				{
					auto m_poly = sss.get_linear_mem_model()->make_energy_poly_coeff(spring_len);
					vec_springs.back().poly.head(m_poly.size()) = m_poly;
				}
				vec_springs.back().poly *= strength;
			}
		}

		// Generate springs

		for (int i = 0; i < smp->FF.rows(); ++i)
		{
			int f0 = smp->FF(i, 0);
			int f1 = smp->FF(i, 1);
			int j0 = smp->fi(i, 0);
			int j1 = smp->fi(i, 1);

			for (int s = 0; s < 2; ++s)
			{
				Eigen::Vector2i sch, scv, scs;
				for (int q = 0; q < 2; ++q)
				{
					Eigen::Vector2i edg;
		
					edg(0) = bprisms[f0].get_base_node_index(q, j0, s);

					scs(q) = vec_edges.size();
					edg(1) = bprisms[f1].get_base_node_index(1 - q, j1, s);
					actuators[i]->segment_springs(q + s * 2) = vec_springs.size();
					vec_springs.emplace_back(vec_edges.size(), (nodes_flt.row(edg(0)) - nodes_flt.row(edg(1))).norm(),
						smp->spth(i, s), (nodes_act.row(edg(0)) - nodes_act.row(edg(1))).norm());
					vec_edges.emplace_back(edg);

					sch(q) = vec_edges.size();
					edg(1) = bprisms[f1].get_base_node_index(q, j1, s);
					vec_edges.emplace_back(edg);

					scv(q) = vec_edges.size();
					edg(0) = bprisms[f0].get_base_node_index(s, j0, q);
					edg(1) = bprisms[f1].get_base_node_index(1-s, j1, 1-q);
					vec_edges.emplace_back(edg);
				}
				vec_scissors.emplace_back(sch);
				vec_scissors.emplace_back(scv);
				//vec_scissors.emplace_back(scs); // no parallel cross springs
			}

			for (int s = 0; s < 2; ++s)
			{
				for (int q = 0; q < 2; ++q)
				{
					vec_flips.emplace_back(
						bprisms[f0].get_base_node_index(1, j0, 0),
						bprisms[f0].get_base_node_index(0, j0, 0),
						bprisms[f0].get_base_node_index(0, j0, 1),
						bprisms[f1].get_base_node_index(q, j1, s)
					);
					Eigen::Vector3d v1 = (nodes_act.row(vec_flips.back()[1]) - nodes_act.row(vec_flips.back()[0])).normalized();
					Eigen::Vector3d v2 = (nodes_act.row(vec_flips.back()[2]) - nodes_act.row(vec_flips.back()[0])).normalized();
					Eigen::Vector3d v3 = nodes_act.row(vec_flips.back()[3]) - nodes_act.row(vec_flips.back()[0]);
					vec_flips_d.push_back(v1.cross(v2).dot(v3));
				}
			}
		}

		// Generate membrane

		if (mode == SmpPhySim::MEM_COARSE)
		{
			int nn = nodes_flt.rows();
			nodes_flt.conservativeResize(nodes_flt.rows() + smp->uv.rows(), 3);
			nodes_flt.bottomLeftCorner(smp->uv.rows(), 2) = smp->uv;
			nodes_flt.bottomRightCorner(smp->uv.rows(), 1).setZero();

			Eigen::MatrixXi Fm(smp->FF.rows() * 2, 3);
			int ind = 0;
			for (int i = 0; i < smp->FF.rows(); ++i)
			{
				Fm.row(ind++) = Eigen::Vector3i{ smp->FF(i, 0), smp->FF(i, 1), nn + smp->F(smp->FF(i, 0), prev3(smp->fi(i, 0))) }.transpose();
				Fm.row(ind++) = Eigen::Vector3i{ smp->FF(i, 1), smp->FF(i, 0), nn + smp->F(smp->FF(i, 1), prev3(smp->fi(i, 1))) }.transpose();
			}
			
			sss.memtri = Fm;
			sss.memrest.resize(sss.memtri.rows());
			for (int i = 0; i < sss.memrest.size(); ++i)
				sss.memrest[i] = Eigen::Matrix2d::Identity() * smp->options.tau;
		}

		if (mode == SmpPhySim::MEM_FINE)
		{
			Eigen::VectorXi base_corner_nodes(m * 6);
			int ind = 0;
			for (int i = 0; i < m; ++i)
				for (int j = 0; j < 6; ++j)
					base_corner_nodes(ind++) = bprisms[i].get_base_node_index(j % 2, j / 2, 0);

			Eigen::MatrixXi Em(smp->FF.rows() * 4, 2);
			Eigen::VectorXi EMark(Em.rows());
			ind = 0;
			for (int i = 0; i < smp->FF.rows(); ++i)
			{
				for (int k = 0; k < 2; ++k)
				{
					EMark(ind) = 100000 + smp->FF(i, k);
					Em.row(ind++) = Eigen::Vector2i(0 + smp->fi(i, k) * 2 + smp->FF(i, k) * 6, 1 + smp->fi(i, k) * 2 + smp->FF(i, k) * 6);
					if (smp->Vb(smp->F(smp->FF(i, k), prev3(smp->fi(i, k)))))
					{
						EMark(ind) = 0;
						Em.row(ind++) = Eigen::Vector2i(1 + smp->fi(i, k) * 2 + smp->FF(i, k) * 6, 0 + smp->fi(i, 1 - k) * 2 + smp->FF(i, 1 - k) * 6);
					}
					else
					{
						EMark(ind) = 100000 + smp->FF(i, k);
						Em.row(ind++) = Eigen::Vector2i(1 + smp->fi(i, k) * 2 + smp->FF(i, k) * 6, 0 + next3(smp->fi(i, k)) * 2 + smp->FF(i, k) * 6);
					}
				}
			}

			// Remove isolated vertices
			Eigen::VectorXi vuni;
			igl::unique(Em, vuni); // vuni must be ordered
			Eigen::VectorXi mem_boundary_nodes;
			mem_boundary_nodes.resizeLike(vuni);
			Eigen::VectorXi vuni_inv(m * 6);
			vuni_inv.setConstant(-1);
			for (int i = 0; i < vuni.size(); ++i)
			{
				mem_boundary_nodes(i) = base_corner_nodes(vuni(i));
				vuni_inv(vuni(i)) = i;
			}
			for (int i = 0; i < Em.rows(); ++i)
			{
				Em(i, 0) = vuni_inv(Em(i, 0));
				Em(i, 1) = vuni_inv(Em(i, 1));
			}

			Eigen::MatrixXd Vm(mem_boundary_nodes.size(), 2);
			for (int i = 0; i < Vm.rows(); ++i)
				Vm.row(i) = nodes_flt.row(mem_boundary_nodes(i)).head(2);

			// Add vertices to open membrane boundary
			if (settings.vert_per_edge > 0)
			{
				ind = Vm.rows();
				int open_edges_count = (EMark.array() == 0).count(); // count membrane open boundary edges
				int vert_per_edge = settings.vert_per_edge; // how many vetrices per open edge to add
				Vm.conservativeResize(Vm.rows() + open_edges_count * vert_per_edge, 2);
				Em.conservativeResize(Em.rows() + open_edges_count * vert_per_edge, 2);
				for (int i = 0; i < EMark.size(); ++i)
				{
					if (EMark(i) > 0) continue;
					int last_vert = Em(i, 1);
					Eigen::Vector2d edge_step = (Vm.row(Em(i, 1)) - Vm.row(Em(i, 0))) / (vert_per_edge + 1);
					Em(i, 1) = ind;
					for (int j = 0; j < vert_per_edge; ++j)
					{
						Vm.row(ind) = Vm.row(Em(i, 0)) + (j + 1) * edge_step.transpose();
						Em(ind, 0) = ind;
						Em(ind, 1) = ind + 1;
						++ind;
					}
					Em(ind - 1, 1) = last_vert;
				}
				EMark.conservativeResize(Em.rows());
				EMark.bottomRows(open_edges_count * vert_per_edge).setZero();
			}

			Eigen::MatrixXd Hf = smp->Fc2d; //TODO: no need for scaling in this implementation
			Eigen::MatrixXd TV;
			Eigen::VectorXi VMark(Vm.rows());
			VMark.setZero();
			Eigen::MatrixXi Fm;
			Eigen::VectorXi VMark2;
			Eigen::VectorXi EMark2;
			igl::triangle::triangulate(Vm, Em, Hf, VMark, EMark, settings.triangle.c_str(), TV, Fm, VMark2, EMark2);

			// Extend nodes
			for (int i = 0; i < TV.rows(); ++i)
			{
				if (VMark2(i) < 100000) continue;
				bprisms[VMark2(i) - 100000].index.push_back(nodes_flt.rows() + i);
			}

			sss.memtri = Fm.array() + nodes_flt.rows();
			sss.memrest.resize(sss.memtri.rows());
			for (int i = 0; i < sss.memrest.size(); ++i)
				sss.memrest[i] = Eigen::Matrix2d::Identity() * smp->options.tau;

			nodes_flt.conservativeResize(nodes_flt.rows() + TV.rows(), 3);
			nodes_flt.bottomLeftCorner(TV.rows(), 2) = TV;
			nodes_flt.bottomRightCorner(TV.rows(), 1).setZero();
		}

		// Generate actuator faces
		Fa.resize(smp->FF.rows() * 8, 4);
		int ind = 0;
		for (int i = 0; i < smp->FF.rows(); ++i)
		{
			int f0 = smp->FF(i, 0);
			int f1 = smp->FF(i, 1);
			for (int s = 0; s < 2; ++s)
			{
				Fa(ind, 0) = bprisms[f0].get_base_node_index(1 - s, smp->fi(i, 0), s);
				Fa(ind, 1) = bprisms[f0].get_base_node_index(    s, smp->fi(i, 0), s);
				Fa(ind, 2) = bprisms[f1].get_base_node_index(    s, smp->fi(i, 1), s);
				Fa(ind, 3) = i;
				++ind;
				Fa(ind, 0) = bprisms[f0].get_base_node_index(    s, smp->fi(i, 0), s);
				Fa(ind, 1) = bprisms[f1].get_base_node_index(1 - s, smp->fi(i, 1), s);
				Fa(ind, 2) = bprisms[f1].get_base_node_index(    s, smp->fi(i, 1), s);
				Fa(ind, 3) = i;
				++ind;
				Fa(ind, 0) = bprisms[f0].get_base_node_index(    s, smp->fi(i, 0),     s);
				Fa(ind, 1) = bprisms[f0].get_base_node_index(    s, smp->fi(i, 0), 1 - s);
				Fa(ind, 2) = bprisms[f1].get_base_node_index(1 - s, smp->fi(i, 1),     s);
				Fa(ind, 3) = i;
				++ind;
				Fa(ind, 0) = bprisms[f0].get_base_node_index(1 - s, smp->fi(i, 0),     s);
				Fa(ind, 1) = bprisms[f1].get_base_node_index(    s, smp->fi(i, 1),     s);
				Fa(ind, 2) = bprisms[f1].get_base_node_index(    s, smp->fi(i, 1), 1 - s);
				Fa(ind, 3) = i;
				++ind;
			}
		}

		sss.nodes = nodes_flt;
		sss.edges.resize(vec_edges.size(), 2);
		sss.linear.resize(vec_springs.size());
		sss.cross.resize(vec_scissors.size(), 2);
		sss.flip.resize(vec_flips.size(), 4);
		sss.flip_d.resize(vec_flips_d.size());
		sss.bodies = std::vector<std::vector<int>>(m);

		sss.bumpers.resizeLike(sss.linear);
		sss.thickness.resize(vec_springs.size());
		sss.segment_initial_lengths.resizeLike(sss.linear);

		sss.base_orig_nodes.resize(m);

		Eigen::VectorXd sim_rest_lengths(vec_springs.size());
		Eigen::VectorXd sim_thickness(vec_springs.size());
		std::vector<Eigen::VectorXd> sim_stiffness(vec_springs.size());

		for (int i = 0; i < vec_edges.size(); ++i)
		{
			sss.edges.row(i) = vec_edges[i];
		}

		for (int i = 0; i < vec_springs.size(); ++i)
		{
			sss.segment_initial_lengths(i) = sim_rest_lengths(i) = vec_springs[i].rest_length;
			sim_thickness(i) = vec_springs[i].thickness;
			sss.bumpers(i) = vec_springs[i].bumper;
			sim_stiffness[i] = vec_springs[i].poly;
			// Membrane mimicking springs
			if (vec_springs[i].bumper < 0)
				sss.segment_initial_lengths(i) = (sss.nodes.row(sss.edges(vec_springs[i].edge, 0)) - sss.nodes.row(sss.edges(vec_springs[i].edge, 1))).norm();
		}

		for (int i = 0; i < vec_springs.size(); ++i)
			sss.linear(i) = vec_springs[i].edge;

		//Eigen::MatrixXi cross(scissors.size(), 2);
		for (int i = 0; i < vec_scissors.size(); ++i)
			sss.cross.row(i) = vec_scissors[i];

		Eigen::VectorXd flip_d(vec_flips_d.size());
		for (int i = 0; i < vec_flips.size(); ++i)
		{
			sss.flip.row(i) = vec_flips[i];
			sss.flip_d(i) = vec_flips_d[i];
		}

		// Assemble body nodes in groups
		for (int i = 0; i < m; ++i)
		{
			sss.bodies[i] = bprisms[i].index;
			sss.base_orig_nodes(i) = bprisms[i].get_orig_index();
		}

		sss.mem_model = mode;

		// Initialize static simulator

		smp_physim->init(sss);
		smp_physim->set_rest_lengths(sim_rest_lengths);
		smp_physim->set_thickness(sim_thickness);
		smp_physim->set_rest_crossratios(Eigen::VectorXd::Zero(sss.cross.rows()));
		smp_physim->set_stiffness_polynomial(sim_stiffness);
		smp_physim->precompute();
		smp_physim->register_callback(&SmpSimulator::physim_callback, this);

		// Generate body elements

		body_vertices.resize(BodyPrism::get_num_vertices() * m, 3);
		body_index.resize(body_vertices.rows());
		body_faces.resize(BodyPrism::get_num_faces() * m, 4);

		Eigen::ArrayXXi bodyf;
		BodyPrism::build_faces(bodyf);

		ind = 0;
		for (int i = 0; i < m; ++i)
		{
			body_faces.block(i * BodyPrism::get_num_faces(), 0, BodyPrism::get_num_faces(), 4) = bodyf;
			body_faces.block(i * BodyPrism::get_num_faces(), 0, BodyPrism::get_num_faces(), 3).array() += ind;
			body_faces.block(i * BodyPrism::get_num_faces(), 3, BodyPrism::get_num_faces(), 1).array() = i;
			ind += BodyPrism::get_num_vertices();
		}

		Eigen::MatrixXd body_cnt = Eigen::MatrixXd::Zero(m, 3);
		for (int i = 0; i < m; ++i)
		{
			for (int k = 0; k < bprisms[i].get_num_nodes(); ++k)
				body_cnt.row(i) += nodes_flt.row(bprisms[i].index[k]);
			body_cnt.row(i) /= bprisms[i].get_num_nodes();
			for (int k = 0; k < BodyPrism::get_num_vertices(); ++k)
			{
				body_vertices.row(k + i * BodyPrism::get_num_vertices()) = body_prisms[i].V.row(k) - body_cnt.row(i);
				body_index(k + i * BodyPrism::get_num_vertices()) = i;
			}
		}

		// Add bracket meshes for visualization
		bracket_vertices.resize(6 * 2 * 4 * smp->FF.rows(), 3);
		bracket_body_index.resize(bracket_vertices.rows());
		bracket_faces.resize(4 * 4 * 4 * smp->FF.rows(), 3);
		Eigen::MatrixXi br_faces(16, 3);
		Eigen::MatrixXi vp(3, 5);
		vp << 4, 10, 11, 5, 4, 2, 8, 9, 3, 2, 0, 6, 7, 1, 0;
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 4; ++j)
			{
				br_faces.row(0 + j * 2 + i * 8) = Eigen::RowVector3i{ vp(i + 1, j + 1), vp(i, j), vp(i, j + 1) };
				br_faces.row(1 + j * 2 + i * 8) = Eigen::RowVector3i{ vp(i + 1, j + 1), vp(i + 1, j), vp(i, j) };
			}
		Eigen::PermutationMatrix<-1, -1> perm_tri(3);
		perm_tri.indices() = Eigen::Vector3i{ 0, 2, 1 };
		for (int i = 0; i < smp->FF.rows(); ++i)
			for (int s = 0; s < 2; ++s)
				for (int k = 0; k < 2; ++k)
				{
					int offs = ((i * 2 + s) * 2 + k) * 12; // offset of current bracket vertices
					for (int t = 0; t < 2; ++t)
					{
						int f0 = smp->FF(i, k);
						int f1 = smp->FF(i, 1 - k);
						int j0 = smp->fi(i, k);
						int j1 = smp->fi(i, 1 - k);

						Eigen::MatrixXd V(6, 3);

						V.row(0) = body_prisms[f0].V.row(0 + j0 * 2 + s * 6);
						V.row(1) = body_prisms[f0].V.row(1 + prev3(j0) * 2 + s * 6) * 0.33 + V.row(0) * 0.67;
						V.row(4) = body_prisms[f1].V.row(1 + j1 * 2 + s * 6);
						V.row(5) = body_prisms[f1].V.row(0 + next3(j1) * 2 + s * 6) * 0.33 + V.row(4) * 0.67;

						V.row(2) = (V.row(0) + V.row(4)) / 2;
						Eigen::Vector3d dir = Eigen::Vector3d::Zero();
						dir.head(2) = Eigen::Rotation2Dd(-M_PI_2) * (V.row(4) - V.row(0)).head(2).normalized().transpose();
						V.row(3) = V.row(2) + dir.transpose() * 0.25;
						double leng = tan(37 * M_PI / 180) * (V.row(0) - V.row(4)).norm() / 2;
						V.row(2) += dir * leng;
						V.row(3) += dir * (leng + 0.4);

						V.rowwise() += Eigen::RowVector3d{ 0, 0, 1 } *t * (2 * s - 1);

						// transform to body coodrinate frame, the ones attached
						V.row(0) -= body_cnt.row(f0);
						V.row(1) -= body_cnt.row(f0);
						V.row(4) -= body_cnt.row(f1);
						V.row(5) -= body_cnt.row(f1);

						bracket_vertices.block(offs + t * 6, 0, 6, 3) = V;
						bracket_body_index.segment(offs + t * 6, 6) << f0, f0, -1, -1, f1, f1;
					}
					bracket_faces.block(k * 16 + s * 32 + i * 64, 0, 16, 3) = (s ? br_faces * perm_tri : br_faces).array() + offs;
				}

		store_frame_info();

		return 0;
	}

	void SmpSimulator::single_async_task(void* user_params)
	{
		std::clog << "Simulation thread started\n";
		for (int frame = get_simulated_frames(); frame < get_max_frames(); ++frame)
		{
			if (task_stop_requested()) break;
			this->simulate_timestep(frame == 1); // frame == 0 is for rest config
		}
		std::clog << "Simulation thread shut down\n";
	}

	void SmpSimulator::store_frame_info()
	{
		auto finfo = std::make_shared<FrameInfo>();
		smp_physim->get_frame_info(*finfo);

		compute_frame_info_derived(finfo.get());

		for (int i = 0; i < settings.timestep_mul; ++i)
			frame_info.emplace_back(finfo);
	}

	void SmpSimulator::compute_frame_info_derived(FrameInfo* finfo)
	{
		// Compute rigid body vertices
		auto vertices = body_vertices;
		for (int i = 0; i < vertices.rows(); ++i)
			vertices.row(i) = vertices.row(i) * finfo->body_rotation[body_index(i)].transpose() + finfo->body_centroids.row(body_index(i));
		finfo->body_vertices = vertices.cast<float>();
		igl::per_vertex_normals(finfo->body_vertices, body_faces, finfo->body_normals);

		static igl::ARAPData arap_data;
		static Eigen::MatrixXd V;
		static Eigen::VectorXi b;

		if (frame_info.empty())
		{
			V = bracket_vertices;
			for (int i = 0; i < V.rows(); ++i)
			{
				if (bracket_body_index(i) < 0) continue;
				V.row(i) += finfo->body_centroids.row(bracket_body_index(i));
			}
			finfo->bracket_vertices = V.cast<float>();

			// Precompute arap data
			b.resize((0 <= bracket_body_index.array()).count());
			int ind = 0;
			for (int i = 0; i < bracket_body_index.size(); ++i)
			{
				if (bracket_body_index(i) < 0) continue;
				b(ind) = i;
				++ind;
			}
			igl::arap_precomputation(V, bracket_faces, 3, b, arap_data);
		}
		else
		{
			// Compute deformation of brackets for visualization
			Eigen::MatrixXd bc(b.size(), 3);
			for (int i = 0; i < b.size(); ++i)
				bc.row(i) = bracket_vertices.row(b(i)) * finfo->body_rotation[bracket_body_index(b(i))].transpose() + finfo->body_centroids.row(bracket_body_index(b(i)));
			// Find any prior deformation of visualized brackets.
			Eigen::MatrixXd U;
			for (int i = frame_info.size() - 1; i >= 0; --i) {
				if (0 < frame_info[i]->bracket_vertices.size()) {
					U = frame_info[i]->bracket_vertices.cast<double>();
					break;
				}
			}
			if (U.size() == 0)
				U = bracket_vertices.cast<double>();
			igl::arap_solve(bc, arap_data, U);
			finfo->bracket_vertices = U.cast<float>();
		}
		igl::per_vertex_normals(finfo->bracket_vertices, bracket_faces, finfo->bracket_normals);

		if (smp_physim->memtri.size() != 0)
			igl::per_vertex_normals(finfo->nodes.cast<float>(), smp_physim->memtri, finfo->mem_normals);
	}

	int SmpSimulator::simulate_timestep(bool dry)
	{		
		double t = t_simulated + t_step * settings.timestep_mul;
		double plastic_fraction = (dry ? 0.0 : settings.plastic_fraction);
		smp_physim->update_brackets(t, plastic_fraction); // no plasticity for the dry simulation step

		if (smp_physim->step(t_step * settings.timestep_mul, dry))
		{
			store_frame_info();
			frame_info.back()->time = t;
			frame_info.back()->plastic_fraction = plastic_fraction;
			t_simulated = t;
			return 0;
		}

		return 1;
	}

}
