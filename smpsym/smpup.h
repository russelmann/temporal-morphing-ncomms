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

#ifndef SMPUP_H
#define SMPUP_H

#include "utils/umath.h"
#include "interface.h"
#include "smp_mesh.h"
#include "utils/umath.h"
#include "utils/utils.h"

#include "igl/readOFF.h"
#include "igl/writeSTL.h"

#include "utils/single_task_thread.h"

namespace smpup {

	enum SMP_POS { SMP_0, SMP_1, SMP_FLT, SMP_ACT };

	struct SmpBuildParameters
	{
		std::string local_name;
		SMP_POS pos;
		const SmpPrintOptions* print_options;
		float* base_done;
		SmpBuildParameters(std::string local_name, SMP_POS pos, const SmpPrintOptions* print_options, float* base_done)
			: local_name(local_name), pos(pos), print_options(print_options), base_done(base_done)
		{ }
	};

	class SmpData : public SmpMesh
	{
	public:
		std::string name;
		std::string working_dir;

		Eigen::VectorXd time_landscape;
		Eigen::MatrixXd spth;      // #FF by 2: spring thicknesses per actuator sides
		std::set<int> base_align;  // #align: bases used for alignment rods

		template<typename Archive>
		void serialize(Archive& ar)
		{
			ar(cereal::make_nvp("name", name));
			ar(cereal::make_nvp("time_landscape", time_landscape));
			ar(cereal::make_nvp("spth", spth));
			ar(cereal::make_nvp("base_align", base_align));
			SmpMesh::serialize(ar);
		}
	};

	// Main Smp class
	class Smp : public SmpData, SingleAsyncTask
	{
	private:
		const char schar[2] = { 'o', 'u' };

		// Write a single node to file
		int write_node(std::FILE* fid, int i, int s, bool act, const SmpPrintOptions* print_options) const;

		static void merge_meshes(const std::vector<Eigen::MatrixXd>& Vm, const std::vector<Eigen::MatrixXi>& Fm, Eigen::MatrixXd& VV, Eigen::MatrixXi& FF);

		void single_async_task(void* user_params = nullptr) const;

	public:
		std::shared_ptr<const CoreSettings> core_settings;

		Smp() { }

		// Copy constructor, only SmpData is copied
		Smp(const Smp& smp) : SmpData(smp)
		{
			core_settings = smp.core_settings;
		}

		Smp& operator = (const Smp& smp)
		{
			core_settings = smp.core_settings;
			*(SmpData*)this = smp; return *this;
		}

		// Read stencil from file and compute all internal parameters
		bool read_mesh(const std::string full_name);

		// Write stencil to OBJ file
		void write_obj(const std::string full_name) const;

		// Select alignment bases (for fabrication)
		void select_align();

		// Export SCAD file in working directory
		void export_scad(std::string local_name, const SmpPrintOptions* print_options = nullptr) const;

		inline bool is_build_in_progress() const { return is_running_task(); }

		// Build mesh for printing by running OpenSCAD in working directory
		void build_scad(std::string local_name, SMP_POS pos, const SmpPrintOptions* print_options = nullptr, float* base_done = nullptr) const
		{
			SmpBuildParameters* smp_build_parameters = new SmpBuildParameters(local_name, pos, print_options, base_done);
			run_task(smp_build_parameters);
		}

		// Stop building OpenSCAD
		inline void stop_build() const { shut_down_task(); }

		// Build support frame in working directory
		void build_support_frame(int side) const;

		// Create regular structure
		void init_regular(int rows, int cols);

		// Create regular structure
		void init_single();

		~Smp()
		{
			shut_down_task();
		}

		template<typename Archive>
		void serialize(Archive& ar)
		{
			SmpData::serialize(ar);
		}
	};
	
}

#endif
