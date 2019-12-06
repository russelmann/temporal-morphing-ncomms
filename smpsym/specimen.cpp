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

#include "specimen.h"
#include "smp_addons.h"
#include "utils/utils.h"

#include "igl/readSTL.h"
#include "igl/writeSTL.h"

#include <fstream>
#include <vector>

namespace smpup {

	void specimen_block(std::string openscad_path, std::string dir)
	{
		//dir = "C:/Research/smpup/source/smpsym/";

		int number_copies = 5;
		std::vector<double> lengths{ 6, 7, 8 };
		std::vector<double> thickness{ 0.4, 0.5, 0.6 }; // we should do a smarter selection of thicknesses vs lengths

		std::string scadname = dir + "scad/bracket.scad";

		std::ofstream ofs;

		for (int i = 0; i < lengths.size(); ++i)
		{
			ofs.open(dir + "scad/inputs.scad");
			ofs << "spring_th = 0.4;\n";
			ofs << "d = " << lengths[i] << ";\n";
			ofs.close();

			std::string stlname = dir + "specimen/bracket_" + std::to_string(lengths[i]) + ".stl";

			run_process('"' + openscad_path + '"' + ' ' + scadname + " -o " + stlname);

			Eigen::MatrixXd V;
			Eigen::MatrixXi F;
			Eigen::MatrixXd N;
			igl::readSTL(stlname, V, F, N);
			std::remove(stlname.c_str());

			V.col(0).array() += i * 50;
			for (int j = 0; j < number_copies; ++j)
			{
				V.col(1).array() += 20;
				igl::writeSTL(stlname.substr(0, stlname.size()-4) + "_" + std::to_string(j) + ".stl", V, F, N);
			}
		}

		std::remove((dir + "scad/inputs.scad").c_str());
	}

}