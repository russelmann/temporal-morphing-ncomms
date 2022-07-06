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

#ifndef SMPDATA_H
#define SMPDATA_H

#include "smp_addons.h"
#include "bracket_model.h"
#include "utils/self_serializer.h"
#include <set>

namespace smpup
{

	struct CoreSettings : public cereal::SelfSerializer<CoreSettings>
	{
		int num_cores = 4; // Number of cores to use for thread pool
#ifdef _WIN32
        std::string openscad_path = "C:/Program Files/OpenSCAD/openscad.exe";
		std::string smpup_path = filesystem::canonical("../").string();
#else
        std::string openscad_path = "/Applications/OpenSCAD.app/Contents/MacOS/./OpenSCAD";
        std::string smpup_path = "../../..";
#endif

		float timestep = 1; // sec

		CoreSettings(filesystem::path default_dir = filesystem::path()) : SelfSerializer("core_settings", default_dir) { }

		filesystem::path get_openscad_path() const { return openscad_path; }
		filesystem::path get_smpup_path() const { return smpup_path; }

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(cereal::make_nvp("num_cores", num_cores));
			ar(cereal::make_nvp("openscad_path", openscad_path));
			ar(cereal::make_nvp("smpup_path", smpup_path));
			ar(cereal::make_nvp("time step", timestep));
		}
	};

	struct SmpPrintOptions
	{
		bool markers = true;       // output markers
		int handle = -1;           // handle index, -1 if none
		bool alignment = true;     // output alignment holes
		float align_diam = 0.75;   // diameter of alignment holes

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(cereal::make_nvp("markers", markers));
			ar(cereal::make_nvp("alignment diameter", align_diam));
			ar(cereal::make_nvp("handle", handle));
			ar(cereal::make_nvp("alignment", alignment));
		}
	};

	struct SmpOptions : public cereal::SelfSerializer<SmpOptions>
	{
		float base_r      = 2;     // base radius
		float base_th     = 2;     // base thickness
		float base_p      = 0.3;   // base indent
		float bracket_w   = 1;     // bracket width
		float min_gap     = 0.1;   // minimal gap between bumpers in the flat state
		float min_bumper  = 0.1;   // minimal bumper size

		float tau         = 3;     // membrane stretch factor
		float lame2       = 0.72;  // Shear modulus (second Lame parameter), MPa
		float mem_h       = 0.5;   // membrane thickness, mm
		float scissor     = 0.1;   // scissor stiffness

		SmpPrintOptions default_print_options;

		SmpOptions(filesystem::path default_dir = "") : SelfSerializer("smp_options", default_dir)
		{ }

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(cereal::make_nvp("base radius", base_r));
			ar(cereal::make_nvp("base thickness", base_th));
			ar(cereal::make_nvp("base indent", base_p));
			ar(cereal::make_nvp("bracket width", bracket_w));
			ar(cereal::make_nvp("minimal gap", min_gap));
			ar(cereal::make_nvp("minimal bumper", min_bumper));
			ar(cereal::make_nvp("membrane stretch factor", tau));
			ar(cereal::make_nvp("membrane shear modulus", lame2));
			ar(cereal::make_nvp("membrane thickness", mem_h));
			ar(cereal::make_nvp("scissor", scissor));
			ar(cereal::make_nvp("default print options", default_print_options));
		}
	};

	struct SimulatorSettings
	{
		// Solver settings
		float  tol = 1e-3;
		int    max_iter_dry = 1000;
		int    max_iter_wet = 500;
		float  max_step_size = 10;
		bool   freeze_centroid = true;
		bool   compute_condition = false; // compute condition number after each time step (time costly)

		template<typename Archive>
		void serialize(Archive& ar)
		{
			ar(cereal::make_nvp("tol", tol));
			ar(cereal::make_nvp("max_iter_dry", max_iter_dry));
			ar(cereal::make_nvp("max_iter_wet", max_iter_wet));
			ar(cereal::make_nvp("max_step_size", max_step_size));
			ar(cereal::make_nvp("freeze_centroid", freeze_centroid));
			ar(cereal::make_nvp("compute_condition", compute_condition));
		}
	};

	class StaticSimulatorStruct
	{
	private:
		SmpOptions        options;
		SimulatorSettings settings;
		std::shared_ptr<LinearMemModel> linear_mem_model;
		std::shared_ptr<BracketModel> bracket_model;

	public:
		enum MEM_MODEL { MEM_CONSTANT, MEM_LINEAR, MEM_COARSE, MEM_FINE, MEM_MODELS };
		MEM_MODEL mem_model = MEM_FINE;

		static const char* MEM_MODEL_COMBO()
		{
			static const char list[] = "Constant\0Linear\0Coarse\0Fine\0\0";
			return list;
		}

		// Caution: this function is not returning const value since imgui requires non-const input
		static std::string& MEM_MODEL_NAME(int i)
		{
			static std::vector<std::string> vstr;
			if (vstr.size() == 0)
			{
				const char* mem_model_combo = MEM_MODEL_COMBO();
				vstr.reserve(MEM_MODELS);
				vstr.emplace_back(mem_model_combo);
				for (int i = 1; vstr.size() < MEM_MODELS; ++i)
					if (*mem_model_combo++ == '\0') vstr.emplace_back(mem_model_combo);
			}
			return vstr[i];
		}

		// Initialization data
		Eigen::MatrixXd nodes;     // #nodes by 3 : simulation nodes
		Eigen::MatrixXi edges;     // #edges by 2 : edges for linear and scissor springs
		Eigen::VectorXi linear;    // #linear by 1 : edge indices for segment springs
		Eigen::MatrixXi cross;     // #cross by 2 : edge indices for cross springs
		Eigen::MatrixXi flip;      // #flip by 4 : node indices for flip elements
		Eigen::VectorXd flip_d;    // #flip by 1 : distances for flip elements
		std::vector<std::vector<int>> bodies;  // #bodies by #nodes_per_body : rigid bodies node indices
		Eigen::MatrixXi memtri;                // #membrane triangles by 3 : membrane triangulation node indices
		std::vector<Eigen::MatrixXd> memrest;  // #membrane triangles : rest strains for triangles

        // Operation data
		Eigen::VectorXd bumpers;                   // #linear by 1 : edge bumper distances (minimal allowed length)
		Eigen::VectorXd thickness;                 // #linear by 1 : bracket thicknesses / membrane mimics strength, e.g. for boundary weakening
		Eigen::VectorXd segment_initial_lengths;   // #linear by 1 : bracket original rest length

		// Segment spring membrane special parameters
		Eigen::VectorXi mspring_edges;     // Edges of segment springs mimicing the membrane
		Eigen::VectorXi base_orig_nodes;   // #F : nodes of base origings

		StaticSimulatorStruct(std::shared_ptr<LinearMemModel> linear_mem_model = nullptr, std::shared_ptr<BracketModel> bracket_model = nullptr)
			: linear_mem_model(linear_mem_model), bracket_model(bracket_model) { }

		StaticSimulatorStruct(SmpOptions options, SimulatorSettings settings, std::shared_ptr<LinearMemModel> linear_mem_model, std::shared_ptr<BracketModel> bracket_model)
			: options(options), settings(settings), linear_mem_model(linear_mem_model), bracket_model(bracket_model) { }

		const SmpOptions& get_smp_options() const { return options; }
		void set_smp_options(const SmpOptions& smp_options) { options = smp_options; }
		const SimulatorSettings& get_settings() const { return settings; }
		void set_settings(const SimulatorSettings& simulator_settings) { settings = simulator_settings; }

		std::shared_ptr<LinearMemModel> get_linear_mem_model() const { return linear_mem_model; }
		void set_linear_mem_model(std::shared_ptr<LinearMemModel> linear_mem_model) { this->linear_mem_model = linear_mem_model; }

		std::shared_ptr<BracketModel> get_bracket_model() const { return bracket_model; }
		void set_bracket_model(std::shared_ptr<BracketModel> bracket_model) { this->bracket_model = bracket_model; }
	};

	struct SmpSettings
	{
		SimulatorSettings simulator_settings;

		Eigen::Vector3f color;
		StaticSimulatorStruct::MEM_MODEL mem_model = StaticSimulatorStruct::MEM_FINE;
		std::string triangle = "a5qQ";
		int vert_per_edge = 0;
		float plastic_fraction = 0.2;
		bool mimic_bounary_weakening = true;
		float boundary_weakening = 0.6;
		float constant_force = 0;
		int timestep_mul = 1;

		SmpSettings()
		{
			color = Eigen::Vector3f::Random().array().abs();
			if (color.norm() < 0.5) color /= color.norm() * 2.0f;
		}

		template<typename Archive>
		void serialize(Archive& ar)
		{
			ar(cereal::make_nvp("simulator_settings", simulator_settings));
			ar(cereal::make_nvp("color", color));
			ar(cereal::make_nvp("mem_model", mem_model));
			ar(cereal::make_nvp("triangle", triangle));
			ar(cereal::make_nvp("vert_per_edge", vert_per_edge));
			ar(cereal::make_nvp("plastic_fraction", plastic_fraction));
			ar(cereal::make_nvp("mimic_bounary_weakening", mimic_bounary_weakening));
			ar(cereal::make_nvp("boundary_weakening", boundary_weakening));
			ar(cereal::make_nvp("constant_force", constant_force));
			ar(cereal::make_nvp("timestep_mul", timestep_mul));
		}
	};

}

#endif
