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

#include "smpup.h"

#include "igl/readOBJ.h"
#include "igl/boundary_loop.h"

#include <algorithm>

namespace smpup {

	void Smp::merge_meshes(const std::vector<Eigen::MatrixXd>& Vm, const std::vector<Eigen::MatrixXi>& Fm, Eigen::MatrixXd& VV, Eigen::MatrixXi& FF)
	{
		int Vcount = 0;
		int Fcount = 0;

		for (int i = 0; i < Vm.size(); ++i)
		{
			Vcount += Vm[i].rows();
			Fcount += Fm[i].rows();
		}

		VV.resize(Vcount, 3);
		FF.resize(Fcount, 3);
		int Vco = 0;
		int Fco = 0;
		for (int i = 0; i < Vm.size(); ++i)
		{
			if (Vm[i].rows() == 0) continue;
			VV.block(Vco, 0, Vm[i].rows(), 3) = Vm[i];
			FF.block(Fco, 0, Fm[i].rows(), 3) = Fm[i].array() + Vco;
			Vco += Vm[i].rows();
			Fco += Fm[i].rows();
		}
	}

	void Smp::single_async_task(void* user_params) const
	{
		SmpBuildParameters smp_build_parameters = *(SmpBuildParameters*)user_params;
		for (int s = 0; s < 2; ++s)
		{
			if (task_stop_requested()) break;

			std::vector<Eigen::MatrixXd> Vm(F.rows());
			std::vector<Eigen::MatrixXi> Fm(F.rows());

			int nthreads = 6;
			std::vector<std::vector<int>> vth(nthreads);
			for (int i = 0; i < Vm.size(); ++i)
				vth[i % nthreads].push_back(i);

			std::vector<std::thread> thpool;
			thpool.reserve(nthreads);
			for (int j = 0; j < nthreads; ++j)
			{
				thpool.emplace_back([this, s, smp_build_parameters](std::vector<int>& vs, std::vector<Eigen::MatrixXd>* Vm, std::vector<Eigen::MatrixXi>* Fm) {
					for (int i = 0; i < vs.size(); ++i)
					{
						if (task_stop_requested()) break;
			
						std::string tmpscadname = working_dir + "/tmp_s" + std::to_string(s) + "d_" + std::to_string(vs[i]) + ".scad";
						std::FILE* fid = std::fopen(tmpscadname.c_str(), "wt");
						write_node(fid, vs[i], s, smp_build_parameters.pos == SMP_ACT, smp_build_parameters.print_options);
						fclose(fid);
						std::string tmpoffname = working_dir + "/tmp_s" + std::to_string(s) + "d_" + std::to_string(vs[i]) + ".off";
			
						run_process('"' + core_settings->openscad_path + '"' + ' ' + tmpscadname + " -o " + tmpoffname);
			
						igl::readOFF(tmpoffname, (*Vm)[vs[i]], (*Fm)[vs[i]]);
						remove(tmpscadname.c_str());
						remove(tmpoffname.c_str());
			
						if (smp_build_parameters.base_done) smp_build_parameters.base_done[vs[i]] += 1;
					}
				}, std::ref(vth[j]), &Vm, &Fm);
			}
			
			for (auto &th : thpool)
				th.join();

			if (task_stop_requested()) break;

			Eigen::MatrixXd VV;
			Eigen::MatrixXi FF;
			merge_meshes(Vm, Fm, VV, FF);

			if (s == 0 && smp_build_parameters.pos == SMP_FLT) VV.rightCols(2) *= -1;
			igl::writeSTL(working_dir + "/" + smp_build_parameters.local_name + "_" + std::to_string(smp_build_parameters.pos) + "d" + schar[s] + ".stl", VV, FF);

			if (smp_build_parameters.pos == SMP_FLT && smp_build_parameters.print_options->markers)
			{
				// Markers
				std::string tmpscadname = working_dir + "/tmp_marker.scad";
				std::string tmpoffname = working_dir + "/tmp_marker.off";
				std::FILE* fid = std::fopen(tmpscadname.c_str(), "wt");
				fprintf(fid, "cylinder(%f,%f,%f,$fs=0.1);\n", 0.5, 0.3, 0.3);
				std::fclose(fid);

				run_process('"' + core_settings->openscad_path + '"' + ' ' + tmpscadname + " -o " + tmpoffname);
				Eigen::MatrixXd Vmark;
				Eigen::MatrixXi Fmark;
				igl::readOFF(tmpoffname, Vmark, Fmark);
				remove(tmpscadname.c_str());
				remove(tmpoffname.c_str());

				std::vector<Eigen::MatrixXd> Vm(F.rows());
				std::vector<Eigen::MatrixXi> Fm(F.rows());
				for (int i = 0; i < Vm.size(); ++i)
				{
					if (smp_build_parameters.print_options->alignment && base_align.find(i) != base_align.end()) continue;
					Vm[i] = Vmark;
					Fm[i] = Fmark;
					Vm[i].array().col(0) += Fc2d(i, 0);
					Vm[i].array().col(1) += Fc2d(i, 1);
					Vm[i].array().col(2) += -0.5 / 2 + (options.base_th / 2 - 0.5 / 2 + 0.002) * (1 - 2 * s);
					if (s == 0) Vm[i].rightCols(2) *= -1;
				}

				Eigen::MatrixXd VV;
				Eigen::MatrixXi FF;
				merge_meshes(Vm, Fm, VV, FF);
				if (VV.size() != 0) igl::writeSTL(working_dir + "/" + smp_build_parameters.local_name + "_" + schar[s] + "_markers.stl", VV, FF);
			}
		}
	}

	bool Smp::read_mesh(const std::string fname)
	{
		filesystem::path path(fname);
		if (path.extension().compare(".obj"))
		{
			std::clog << "Error: not an OBJ file\n";
			return false;
		}
		Eigen::MatrixXd V, TC, N;
		Eigen::MatrixXi F, FTC, FN;
		if (!igl::readOBJ(fname, V, TC, N, F, FTC, FN))
		{
			std::clog << "Error reading file <" << fname << ">\n";
			return false;
		}
		if (FTC.size() != 0)
		{
			Eigen::VectorXi perm = Eigen::VectorXi::Constant(V.rows(), -1);
			for (int i = 0; i < FTC.rows(); ++i)
				for (int j = 0; j < FTC.cols(); ++j)
				{
					if (perm(FTC(i, j)) == -1) perm(FTC(i, j)) = F(i, j);
					else if (perm(FTC(i, j)) != F(i, j)) return false; // wrong uv
				}
			if (perm.minCoeff() < 0) return false; // wrong uv
			TC = perm.asPermutation() * TC;
		}
		init(V, F, &TC);
		name = path.stem().string();
		working_dir = path.parent_path().string();
		select_align();
		//spth = Eigen::MatrixXd::Ones(FF.rows(), 2) * 0.4;
		spth = Eigen::MatrixXd::Zero(FF.rows(), 2);
		time_landscape = Eigen::VectorXd::Zero(V.rows());
		return true;
	}

	void Smp::write_obj(const std::string fname) const
	{
		// Output OBJ with flat configuration
		std::ofstream s(fname);
		if (!s.is_open())
		{
			fprintf(stderr, "IOError: writeOBJ() could not open %s\n", fname.c_str());
			return;
		}
		s << V.format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, " ", "\n", "v ", "", "", "\n"));
		s << uv.format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, " ", "\n", "vt ", "", "", "\n"));
		s << (F.array() + 1).format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, " ", "\n", "f ", "", "", "\n"));
		s.close();
	}

	void Smp::select_align()
	{
		std::vector<int> qalign(4, -1);

		int pn = (ff.array() < 0).rowwise().any().count();
		std::vector<int> pind;
		pind.reserve(pn);
		for (int i = 0; i < ff.rows(); ++i)
			if ((ff.row(i).array() < 0).any()) pind.push_back(i);
		Eigen::MatrixXd pdist(pn, pn);
		pdist.setZero();
		Eigen::Vector2i maxij(-1, -1);
		double dmax = 0;
		for (int i = 0; i < pind.size(); ++i)
		{
			for (int j = i + 1; j < pind.size(); ++j)
			{
				pdist(i, j) = pdist(j, i) = (Fc2d.row(pind[i]) - Fc2d.row(pind[j])).norm();
				if (dmax < pdist(i, j))
				{
					maxij = Eigen::Vector2i{ i, j };
					qalign[0] = pind[i];
					qalign[2] = pind[j];
					dmax = pdist(i, j);
				}
			}
		}

		for (int i = 0; i < pind.size(); ++i)
			if (pdist(i, maxij(0)) < 10 || pdist(i, maxij(1)) < 10)
				pind[i] = -1;

		dmax = 0;
		for (int i = 0; i < pind.size(); ++i)
		{
			if (pind[i] < 0) continue;
			for (int j = i + 1; j < pind.size(); ++j)
			{
				if (pind[j] < 0) continue;
				if (dmax < pdist(i, j))
				{
					qalign[1] = pind[i];
					qalign[3] = pind[j];
					dmax = pdist(i, j);
				}
			}
		}

		base_align.clear();
		for (auto i : qalign)
			if (-1 < i) base_align.insert(i);
	}

	int Smp::write_node(std::FILE* fid, int i, int s, bool act, const SmpPrintOptions* print_options) const
	{
		// N.B.: rotation in openscad is the opposite to Eigen convention, use rotate(a = %f*180/PI, v = -[%f, %f, %f]); with negative axis vector
		// N.B.: side in matlab is here is 1 - 2 * s !!!

		double eps = 1e-3;

		const Eigen::MatrixX3d& Fcs = (s == 0) ? Fco : Fcu;

		//TODO: is dr, a randomly chosen direction a good way to find rotation of the base in flat layout
		int dr = 0;
		while (ff(i, dr) < 0) ++dr;

		Eigen::Vector3d fc2dir;
		fc2dir.head(2) = Fc2d.row(ff(i, dr)) - Fc2d.row(i);
		fc2dir(2) = 0;
		Eigen::Matrix3d R = vectors2rotation(fnorm.row(i), Fcs.row(ff(i, dr)) - Fcs.row(i), Eigen::Vector3d{ 0, 0, 1 }, fc2dir);
		Eigen::AngleAxisd rot(R);

		fprintf(fid, "\n\n// Node %04d\n\n", i);

		if (act)
		{
			fprintf(fid, "translate([%f, %f, %f]) {\n", Fcs(i, 0), Fcs(i, 1), Fcs(i, 2));
		}
		else
		{
			fprintf(fid, "translate([%f, %f, 0]) {\n", Fc2d(i, 0), Fc2d(i, 1));
			fprintf(fid, "rotate(a = %f*180/PI, v = -[%f, %f, %f]) {\n", rot.angle(), rot.axis()(0), rot.axis()(1), rot.axis()(2));
		}

		fprintf(fid, "difference() {\n");
		fprintf(fid, "union() {\n");


		// Base
		Eigen::AngleAxisd rot_cnt(vectors2rotation(Eigen::Vector3d{ 0, 0, 1 }, Eigen::Vector3d{ 1, 0, 0 }, fnorm.row(i), Fcs.row(ff(i, dr)) - Fcs.row(i)));
		fprintf(fid, "rotate(a = %f*180/PI, v = -[%f, %f, %f]) {\n", rot_cnt.angle(), rot_cnt.axis()(0), rot_cnt.axis()(1), rot_cnt.axis()(2));
		fprintf(fid, "// Base\n");
		fprintf(fid, "translate([0, 0, %f])\n", -options.base_th / 2 - options.base_p * (1 - s) - eps * s);
		fprintf(fid, "cylinder(%f,%f,%f,$fs=0.5);\n", options.base_th + options.base_p + eps, options.base_r, options.base_r);

		// [film]

		fprintf(fid, "}\n");

		// Bumpers
		for (int j = 0; j < 3; ++j)
		{
			if (ff(i, j) < 0) continue;

			Eigen::Vector3d vrf = Fcs.row(ff(i, j)) - Fcs.row(i);
			Eigen::Vector3d vr = vrf - vrf.dot(fnorm.row(i)) * fnorm.row(i).transpose();
			double alph = M_PI / 2 - angle_vectors(vrf, fnorm.row(i));
			double lng = vrf.norm() / cos(alph) * 0.5;

			Eigen::AngleAxisd rot_cyl(vectors2rotation(Eigen::Vector3d{ 0, 0, 1 }, Eigen::Vector3d{ 1, 0, 0 }, vr, fnorm.row(i)));
			fprintf(fid, "// Bumper %d\n", j);
			fprintf(fid, "rotate(a = %f*180/PI, v = -[%f, %f, %f]) {\n", rot_cyl.angle(), rot_cyl.axis()(0), rot_cyl.axis()(1), rot_cyl.axis()(2));
			fprintf(fid, "difference() {\n");
	//		fprintf(fid, "difference() {\n"); // anti-glossy
			//fprintf(fid, 'intersection() {\n');
			//fprintf(fid, 'scale([1, 1.5, 1])\n');
			//fprintf(fid, 'cylinder(%f,%f,%f,$fs=0.2);\n', 2 * lng, base_h / 2, base_h / 2); % too long, but cut away further
			fprintf(fid, "rotate([0,0,45])\n");
			fprintf(fid, "cylinder(%f,%f,%f,$fn=4);\n", 2 * lng, sqrt(2) * options.base_th / 2, sqrt(2) * options.base_th / 2); // too long, but cut away further
																																//fprintf(fid, '}\n');

			double ang = M_PI_2 - angle_vectors(fnorm.row(ff(i, j)), Fcs.row(i) - Fcs.row(ff(i, j)));

            bool overhang = !((ang < 0) ^ s); // if bumper has an overhang, used for turning flat structure into a height field
			fprintf(fid, "translate([0,0,%f])\n", lng - (overhang ? options.base_th / 2 * tan(abs(ang)) : 0) - 0.05); // -0.1 printing gap
			fprintf(fid, "rotate([0,%f*180/PI,0])\n", overhang ? 0 : ang);
			fprintf(fid, "translate([-5,-5,0])\n");
			fprintf(fid, "cube([10,10,10]);\n");
			fprintf(fid, "}\n");

			// Anti-glossy artifact
	//		fprintf(fid, "// Anti-glossy artifact\n");
	//		fprintf(fid, "translate([0,0,%f])\n", -0.1);
	//		fprintf(fid, "translate([0,0,%f])\n", lng + options.base_th / 2 * tan(abs(ang)) * (overhang ? -1 : 1));
	//		fprintf(fid, "translate([%f,0,0])\n", options.base_th / 2 * (1 - 2 * s));
	//		fprintf(fid, "rotate([0,%f*180/PI,0])\n", M_PI_4 * (1 - 2 * s));
	//		fprintf(fid, "translate([%f,%f,0])\n", -options.base_th / 2, -options.base_th / 2);
	//		fprintf(fid, "cube([%f,%f,1]);\n", options.base_th, options.base_th);
	//		fprintf(fid, "}\n");

			// Small plates enforcing colliding parts as matte
			//         if pos == 2
			// fprintf(fid, 'translate([%f,%f,%f])\n', -base_p / 4 + (-base_h / 2 - base_p / 2 - base_p * 3 / 8) * side, -base_h / 2, base_r + base_p);
			//             fprintf(fid, 'cube([%f,%f,%f]);\n', base_p / 4, base_h, norm(Fc2d(ff(i, j), :) - Fc2d(i, :)) / 2 - base_r);
			//         end

			fprintf(fid, "}\n");
		}

		fprintf(fid, "}\n");


		// Voids to decrease mean mass density of the structure
		//if ((qalign.array() != i).all())
		//{
		//	fprintf(fid, "// Mass density void\n");
		//	fprintf(fid, "rotate(a = %f*180/PI, v = -[%f, %f, %f])\n", rot_cnt.angle(), rot_cnt.axis()(0), rot_cnt.axis()(1), rot_cnt.axis()(2));
		//	double hh = 1.4; //TODO: this needs to be comuted via total structure volume and mass density
		//	fprintf(fid, "translate([0, 0, %f])\n", -hh / 2 - (options.base_th / 2 + options.base_p - hh / 2 + eps) * (1 - 2 * s));
		//	fprintf(fid, "cylinder(%f,%f,%f,$fs=0.5);\n", hh, options.base_r / (2 * s + 2), options.base_r / (4 - 2 * s));
		//}

		// Voids for markers
		if (print_options->markers)
		{
			fprintf(fid, "// Marker void\n");
			fprintf(fid, "rotate(a = %f*180/PI, v = -[%f, %f, %f])\n", rot_cnt.angle(), rot_cnt.axis()(0), rot_cnt.axis()(1), rot_cnt.axis()(2));
			fprintf(fid, "translate([0, 0, %f])\n", -0.5 / 2 + (options.base_th / 2 - 0.5 / 2 + eps * 2) * (1 - 2 * s));
			fprintf(fid, "cylinder(%f,%f,%f,$fs=0.1);\n", 0.5, 0.3, 0.3);
		}

		// Holes for alignment rods
		if (print_options->alignment && base_align.find(i) != base_align.end())
		{
			fprintf(fid, "// Alignment hole\n");
			fprintf(fid, "rotate(a = %f*180/PI, v = -[%f, %f, %f])\n", rot_cnt.angle(), rot_cnt.axis()(0), rot_cnt.axis()(1), rot_cnt.axis()(2));
			fprintf(fid, "translate([0, 0, %f])\n", -options.base_th / 2 - eps * 10 - options.base_p * (1 - s));
			fprintf(fid, "cylinder(%f,%f,%f,$fs=0.1);\n", options.base_th + options.base_p + eps * 20, print_options->align_diam, print_options->align_diam);
		}

		fprintf(fid, "}\n");
		if (act)
			fprintf(fid, "rotate(a = %f*180/PI, v = -[%f, %f, %f]) {\n", -rot.angle(), rot.axis()(0), rot.axis()(1), rot.axis()(2));
		else
			fprintf(fid, "}\n");

		// Brackets
		for (int j = 0; j < 3; ++j)
		{
			if (ff(i, j) < 0) continue;
			double acr = 0;

			Eigen::Vector2d vr = (Fc2d.row(ff(i, j)) - Fc2d.row(i)).transpose();
			double ang = atan2(vr(1), vr(0));
			double alpha = asin(options.base_th / options.base_r / 2) + 7 * M_PI / 180;
			double lng = (vr.norm() - options.base_r * 2 * cos(alpha)) / 2 / cos(alpha + acr);

			if (act) acr = acos(((Fcs.row(ff(i, j)) - Fcs.row(i)).norm() - options.base_r * 2 * cos(alpha)) / 2 / lng) - alpha;

			double bracket_th = spth(Ff(i, j), s);

			fprintf(fid, "// Bracket %d\n", j);
			for (int k = -1; k < 2; k += 2)
			{
				fprintf(fid, "rotate(180/PI*[0, 0, %f])\n", ang + alpha * k);
				fprintf(fid, "translate([%f, 0, 0])\n", options.base_r);
				fprintf(fid, "rotate(180/PI*[0, 0, %f])\n", acr * k);
				double sadj = (1 - 2 * s) * options.base_th / 2;
				if (s == 0) sadj -= options.bracket_w;
				fprintf(fid, "translate([%f, %f, %f])\n", -bracket_th, -bracket_th / 2, sadj);
				fprintf(fid, "cube([%f, %f, %f]);\n", bracket_th + lng + bracket_th / 2 * tan(alpha), bracket_th, options.bracket_w);
			}
		}

		if (act) fprintf(fid, "}\n");

		fprintf(fid, "}\n");

		// Handle
		if (~act && i == print_options->handle)
		{
			double handle_length = 9;
			double slot_indent = 5;
			double slot_width = 2;

			Eigen::Vector3d v1 = fnorm.row(i);
			int fxi = 0;
			if (ff(i, fxi) >= 0) fxi = 1;
			if (ff(i, fxi) >= 0) fxi = 2;
			int vi0 = F(i, next3(fxi));
			int vi1 = F(i, prev3(fxi));
			Eigen::Vector3d v2;
			v2.head(2) = (uv.row(vi0) - uv.row(vi1)).normalized();
			v2(2) = 0;
			Eigen::Vector3d dir = v1.cross(v2);
			double ang = atan2(dir(1), dir(0));
			fprintf(fid, "// Handle for node %04d\n", i);
			fprintf(fid, "translate([%f, %f, 0])\n", Fc2d(i, 0), Fc2d(i, 1));
			fprintf(fid, "rotate(a = %f*180/PI, v = -[0, 0, 1])\n", -ang);
			fprintf(fid, "difference() {\n");
			fprintf(fid, "translate([0, %f, %f])\n", -options.base_r, -options.base_th / 2 - options.base_p * (1 - s));
			fprintf(fid, "cube([%f, %f, %f]);\n", handle_length, options.base_r * 2, options.base_th + options.base_p);
			fprintf(fid, "{\n");
			fprintf(fid, "translate([0, 0, %f])\n", -options.base_th / 2 - options.base_p);
			fprintf(fid, "cylinder(%f, %f, %f);\n", (options.base_th + options.base_p) * 2, options.base_r - eps, options.base_r - eps);
			fprintf(fid, "translate([%f, %f, %f])\n", handle_length - slot_indent, -options.base_r - eps, (-options.base_th / 2 - eps) * s);
			fprintf(fid, "cube([%f, %f, %f]);\n", slot_width, options.base_r * 2 + eps * 2, options.base_th / 2 + eps);
			fprintf(fid, "}\n");
			fprintf(fid, "}\n");
		}

		return 0;
	}

	void Smp::export_scad(std::string name, const SmpPrintOptions* print_options) const
	{
		if (!print_options) print_options = &options.default_print_options;
		std::FILE* fid3d = std::fopen((working_dir + "/" + name + "_3d.scad").c_str(), "w");
		for (int s = 0; s < 2; ++s)
		{
			std::FILE* fid2d = std::fopen((working_dir + "/" + name + "_2d" + schar[s] + ".scad").c_str(), "w");
			for (int i = 0; i < F.rows(); ++i)
			{
				write_node(fid3d, i, s, true, print_options);
				write_node(fid2d, i, s, false, print_options);
			}
			std::fclose(fid2d);
		}
		std::fclose(fid3d);
	}

	void Smp::build_support_frame(int side) const
	{
		double min_connector_dist = 30;
		double outer_radius = 9;
		double inner_radius = 6;
		double extra_h = 0; // extra thickness for transportation robustness

		Eigen::VectorXi loop;
		igl::boundary_loop(F, loop);

		int pn = (ff.array() < 0).rowwise().any().count();
		std::vector<int> pind;
		pind.reserve(pn);
		for (int i = 0; i < ff.rows(); ++i)
			if ((ff.row(i).array() < 0).any()) pind.push_back(i);

		std::string scad_name = working_dir + "/support_frame_" + schar[side] + ".scad";
		std::FILE* fid = std::fopen(scad_name.c_str(), "w");

		fprintf(fid, "module connector(a, b)\n{\nc = b - a;\n");
	    fprintf(fid, "translate([0, 0, %f])\n", (1 - 2 * side) * (options.base_th/2 - 0.25));
		fprintf(fid, "translate((a + b) / 2)\nrotate(atan2(c[1], c[0]))\n"
			"cube([norm([c[0], c[1]]), 0.5, 0.5], center = true);\n}\n");

		fprintf(fid, "module support_frame()\n{\n");
		//for (int i = 0; i < Fc2d.rows(); ++i)
		//	fprintf(fid, "translate([%f, %f, 0]) circle(r=0.01);\n", Fc2d(i, 0), Fc2d(i, 1));
		fprintf(fid, "polygon([\n");
		for (int i = 0; i < loop.size(); ++i)
			fprintf(fid, "\t[%f, %f]%s\n", uv(loop(i), 0), uv(loop(i), 1), i < loop.size() - 1 ? "," : "");
		fprintf(fid, "]);\n");


		fprintf(fid, "}\n\n");

		if (side == 0) fprintf(fid, "rotate([180, 0, 0]) {\n");

		std::vector<bool> base_connected(pind.size(), true);
		for (int i = 0; i < pind.size(); ++i)
		{
			for (int j = 0; j < i; ++j)
				if (base_connected[j] && (Fc2d.row(pind[i]) - Fc2d.row(pind[j])).norm() < min_connector_dist)
				{
					base_connected[i] = false;
					break;
				}
			if (!base_connected[i]) continue;

			Eigen::Vector2d vcon = Fc2d.row(pind[i]);
			int ind = 0;
			for (; ff(pind[i], ind) >= 0; ++ind);
			Eigen::Vector2d vdir = (uv.row(F(pind[i], next3(ind))) - uv.row(F(pind[i], prev3(ind)))).normalized();
			vdir = Eigen::Rotation2D<double>(M_PI_2) * vdir;
			vcon = vcon + vdir * (options.base_r - 1);
			vdir *= 12 - (options.base_r - 1);
			fprintf(fid, "connector([%f, %f], [%f, %f]);\n", vcon(0), vcon(1), vcon(0) + vdir(0), vcon(1) + vdir(1));
		}

		fprintf(fid, "translate([0, 0, %f])\n", -options.base_th / 2 - (1 - side) * (options.base_p + extra_h));
		fprintf(fid, "linear_extrude(%f)\n", options.base_th + options.base_p + extra_h);
		fprintf(fid, "difference()\n{\n");
		fprintf(fid, "offset(%f)\n", outer_radius);
		fprintf(fid, "support_frame();\n");
		fprintf(fid, "offset(%f)\n", inner_radius);
		fprintf(fid, "support_frame();\n");
		fprintf(fid, "}\n");

		if (side == 0) fprintf(fid, "}\n");

		fclose(fid);

		std::string stl_name = working_dir + "/support_frame_" + schar[side] + ".stl";
		run_process('"' + core_settings->openscad_path + '"' + ' ' + scad_name + " -o " + stl_name);
	}

	void Smp::init_regular(int rows, int cols)
	{
		name = "regular_" + std::to_string(rows) + 'x' + std::to_string(cols);

		if (core_settings)
		{
			working_dir = core_settings->get_smpup_path().append("reg").string();
			if (!filesystem::exists(working_dir)) filesystem::create_directories(working_dir);
		}
		SmpMesh::init_regular(rows, cols);
		options.default_print_options.handle = 2 + 4 * (rows / 2);
		options.default_print_options.alignment = false;
		select_align();
		spth = Eigen::MatrixXd::Zero(FF.rows(), 2);
		time_landscape = Eigen::VectorXd::Zero(V.rows());
	}

	void Smp::init_single()
	{
		name = "single";
		if (core_settings)
		{
			working_dir = core_settings->get_smpup_path().append("reg").string();
			if (!filesystem::exists(working_dir)) filesystem::create_directories(working_dir);
		}

		SmpMesh::init_single();
		select_align();
		spth = Eigen::MatrixXd::Zero(FF.rows(), 2);
		time_landscape = Eigen::VectorXd::Zero(V.rows());
	}

}
