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

#ifndef SMPSYM_H
#define SMPSYM_H

#include "smpup.h"
#include "smp_physim.h"
#include "interface.h"
#include "bracket_configurer.h"

#include "PhySim/Geometry/Node.h"
#include "PhySim/Geometry/Edge.h"
#include "PhySim/Geometry/Cell_Poly.h"
#include "PhySim/Energies/EnergyElement_SpringLinear.h"
#include "PhySim/Energies/EnergyElement_SpringCross.h"
#include "PhySim/Models/Model_MassSpring_SMP.h"
#include "PhySim/Models/Model_Reduction_RBCloud.h"
#include "PhySim/Solvers/OptimSolver_USQP_LS.h"
#include "PhySim/Solvers/OptimProblem_BasicStatic.h"

#include "igl/readOBJ.h"

#include "utils/single_task_thread.h"
#define _USE_MATH_DEFINES
#include <cmath> 

namespace smpup
{

	// Actuator data
	struct Actuator
	{
		int membrane_edge;
		Eigen::VectorXi segment_springs;
		Eigen::VectorXi scissors;
		Actuator()
		{
			segment_springs.resize(4);
			scissors.resize(6);
		}
	};


	// Smp simulator
	class SmpSimulator : public SingleAsyncTask
	{
	private:
		std::string name;

		double t_start;      // starting time
		double t_end;        // end time
		double t_simulated;  // simulated time
		double t_step;       // time step

		std::vector<std::shared_ptr<FrameInfo>> frame_info; // Simulation frame tracking

		std::unique_ptr<Smp> smp;
		std::unique_ptr<SmpPhySim> smp_physim;
		//std::shared_ptr<BracketModel> bracket_model;
		SmpSettings settings;

		// Prism visualization data
		Eigen::VectorXi body_index;
		Eigen::MatrixXd body_vertices;
		Eigen::MatrixXi body_faces;
		Eigen::MatrixXd bracket_vertices;
		Eigen::VectorXi bracket_body_index;
		Eigen::MatrixXi bracket_faces;

		// Info about elements
		std::vector<Actuator*> actuators;    // #FF : indexes of elements for each actuator
		Eigen::MatrixXi Fa;                  // #FF by 4 : actuator faces for drawing and corresponding indexes of actuators back to FF

		void store_frame_info();
		void compute_frame_info_derived(FrameInfo* finfo);
		void single_async_task(void* user_params = nullptr);
		int simulate_timestep(bool dry);

	public:

		SmpSimulator(std::shared_ptr<LinearMemModel> linear_mem_model, std::shared_ptr<BracketModel> bracket_model)// = nullptr)
		{
			smp = std::make_unique<Smp>();
			smp_physim = std::make_unique<SmpPhySim>(linear_mem_model, bracket_model);
		}

		SmpSimulator(double t_start, double t_end, double t_step, SmpSettings settings, std::shared_ptr<LinearMemModel> linear_mem_model, std::shared_ptr<BracketModel> bracket_model)
			: t_start(t_start), t_end(t_end), t_step(t_step), settings(settings)
		{
			smp_physim = std::make_unique<SmpPhySim>();
			smp_physim->set_linear_mem_model(linear_mem_model);
			smp_physim->set_bracket_model(bracket_model);
			t_simulated = -t_step;
			frame_info.reserve(int(t_end / t_step) + 1);
		}

		std::string& get_name() { return name; }
		void set_name(std::string name) { this->name = name; }

		// Load smp object and prepare everything for simulation
		int load_smp(const Smp& smpd);

		// Get pointer to smp interface
		inline const Smp* get_smp() { return smp.get(); }

		//inline std::shared_ptr<const CoreSettings> get_smp_core_settings() { return smp->core_settings; }
		inline void set_smp_core_settings(std::shared_ptr<const CoreSettings> core_settings) { smp->core_settings = core_settings; }
		inline void set_smp_working_dir(std::string working_dir) { smp->working_dir = working_dir; }

		inline const BracketModel& get_bracket_model() { return *smp_physim->get_bracket_model(); }
		inline void set_bracket_model(std::shared_ptr<BracketModel> bracket_model) { smp_physim->set_bracket_model(bracket_model); }

		inline std::set<int>& get_align() { return smp->base_align; }
		
		// Get pointer to physim interface
		inline const SmpPhySim* psim() const { return smp_physim.get(); }

		// Configure bracket thicknesses
		std::unique_ptr<SmpSimulator> configure_brackets(const Eigen::VectorXd& time_landscape)
		{
			Eigen::VectorXd actuator_times(smp->FF.rows());
			for (int i = 0; i < smp->FF.rows(); ++i)
				actuator_times(i) = (time_landscape(smp->F(smp->FF(i, 0), next3(smp->fi(i, 0)))) + time_landscape(smp->F(smp->FF(i, 0), prev3(smp->fi(i, 0))))) * 0.5;

			Smp new_smp(*smp);
			new_smp.time_landscape = time_landscape;
			for (int i = 0; i < new_smp.FF.rows(); ++i)
			{
				//double actuator_length = get_actuator_length_at_frame(i, 0);
				Eigen::Vector2i act_edg = smp_physim->edges.row(smp_physim->linear(actuators[i]->segment_springs(0))).transpose(); // since all actuator brackets are equal length in flat config
				double init_mem_length = (frame_info[0]->nodes.row(act_edg(0)) - frame_info[0]->nodes.row(act_edg(1))).norm();
				Eigen::VectorXd mem_poly = smp_physim->get_linear_mem_model()->make_force_poly_coeff(init_mem_length);

				Eigen::Vector2d target_length;
				for (int j = 0; j < 2; ++j)
				{
					Eigen::Vector2i seg_id = get_actuators()[i]->segment_springs.segment(j * 2, 2);
					target_length(j) = (psim()->bumpers(seg_id(0)) + psim()->bumpers(seg_id(1))) * 0.5;
				}
				for (int j = 0; j < 2; ++j)
				{
					double target_mem_length = target_length.mean();
					Eigen::VectorXd target_mem_length_pow = Eigen::ArrayXd::Constant(2, -target_mem_length).pow(natural_seriesd(2));
					double target_mem_force = -mem_poly.dot(target_mem_length_pow); // negative since the force is along negative direction

					bool weak = new_smp.FFb(i);
					if (settings.mem_model == StaticSimulatorStruct::MEM_LINEAR && !settings.mimic_bounary_weakening) weak = false;
					if (weak) target_mem_force *= settings.boundary_weakening;

					new_smp.spth(i, j) = smp_physim->get_bracket_model()->configure_thickness(init_mem_length, target_length(j), actuator_times(i), target_mem_force / 4.0);
				}
			}
			auto new_simulator = std::make_unique<SmpSimulator>(t_start, t_end, t_step, settings, smp_physim->get_linear_mem_model(), smp_physim->get_bracket_model());
			new_simulator->name = name;
			new_simulator->load_smp(new_smp);
			return new_simulator;
		}

		// Create a new thread and run simulation
		inline void simulate() { run_task(); }
		
		// Interrupt simulation thread if it is running
		inline void stop_simulation() { request_task_stop(); }

		static inline void physim_callback(PhySim::IOptimSolver* solver, void* receiver)
		{
			auto simulator = (SmpSimulator*)receiver;
			if (simulator->task_stop_requested())
				solver->ForceSolveStop() = true;
		}

		// Get settings
		inline const SmpSettings& get_settings() const { return settings; }


		// Get end time
		inline double get_time_end() const { return t_end; }

		// Get simulated time
		inline double get_time_simulated() const { return t_simulated; }

		// Get the timestep
		inline double get_timestep() const { return t_step; }

		// Get number of simulated frames
		inline int get_simulated_frames() const { return frame_info.size(); }

		inline int get_max_frames() const { return int(t_end / t_step) + 1; }

		// Get total planned simulation frames
		inline int get_total_frames() const { return frame_info.capacity(); }

		inline const std::vector<Actuator*>& get_actuators() const
		{
			return actuators;
		}
		
		// Get actuator faces
		inline Eigen::MatrixXi get_actuator_faces() const { return Fa; }


		// Get membrane simulation model
		inline StaticSimulatorStruct::MEM_MODEL get_membrane_model()
		{
			return smp_physim->mem_model;
		}

		// Get membrane simulation model name
		inline std::string& get_membrane_model_name()
		{
			return StaticSimulatorStruct::MEM_MODEL_NAME(smp_physim->mem_model);
		}

		// Get actuator length for the given frame
		inline double get_actuator_length_at_frame(int actuator_id, int frame) const
		{
			return (frame_info[frame]->body_origins.row(smp->FF(actuator_id, 0)) - frame_info[frame]->body_origins.row(smp->FF(actuator_id, 1))).norm();
		}

		// Get actuator lengths for simulated frames
		Eigen::VectorXd get_actuator_lengths(int actuator_id) const
		{
			Eigen::VectorXd actuator_lengths(get_simulated_frames());
			for (int i = 0; i < actuator_lengths.size(); ++i)
				actuator_lengths(i) = get_actuator_length_at_frame(actuator_id, i);
			return actuator_lengths; // compiler RVO is assumed
		}

		// Get actuator lengths for the given frame
		Eigen::VectorXd get_actuator_lengths_at_frame(int frame) const
		{
			Eigen::VectorXd actuator_lengths(smp->FF.rows());
			for (int i = 0; i < actuator_lengths.size(); ++i)
				actuator_lengths(i) = get_actuator_length_at_frame(i, frame);
			return actuator_lengths; // compiler RVO is assumed
		}

		// Get recorded node positions for the given frame, interpolate with next frame if frame_alpha > 0
		Eigen::MatrixXd get_positions_at_frame(int frame, double frame_alpha = 0) const
		{
			Eigen::MatrixXd pos;

			frame = std::min(frame, get_simulated_frames() - 1);
			pos = frame_info[frame]->nodes;

			if (frame_alpha != 0)
			{
				int frame_b = std::min(frame + 1, get_simulated_frames() - 1);
				if (frame_b != frame) pos = (1 - frame_alpha) * pos + frame_alpha * frame_info[frame_b]->nodes;
			}

			return pos;
		}

		// Get recorded total strain of scissor springs for the given frame
		double get_cross_strain_at_frame(int frame) const
		{
			frame = std::min(frame, get_simulated_frames() - 1);
			return frame_info[frame]->total_cross_strain;
		}

		// Get linear springs from node data recorded for the given frame: rowwise from a to b and collision flags
		// interpolate with next frame if frame_alpha > 0
		void get_springs_at_frame(Eigen::MatrixXd& a, Eigen::MatrixXd& b, Eigen::VectorXi& collided, int frame, double frame_alpha = 0) const
		{
			frame = std::min(frame, get_simulated_frames() - 1);
			a.resize(smp_physim->linear.rows(), 3);
			b.resizeLike(a);
			for (int i = 0; i < a.rows(); ++i)
			{
				a.row(i) = frame_info[frame]->nodes.row(smp_physim->edges(smp_physim->linear(i), 0));
				b.row(i) = frame_info[frame]->nodes.row(smp_physim->edges(smp_physim->linear(i), 1));
			}
			collided = frame_info[frame]->collided;

			if (frame_alpha == 0) return;

			int frame_b = std::min(frame + 1, get_simulated_frames() - 1);
			if (frame_b == frame) return;

			Eigen::MatrixXd a_b;
			Eigen::MatrixXd b_b;
			a_b.resizeLike(a);
			b_b.resizeLike(b);
			for (int i = 0; i < a.rows(); ++i)
			{
				a_b.row(i) = frame_info[frame_b]->nodes.row(smp_physim->edges(smp_physim->linear(i), 0));
				b_b.row(i) = frame_info[frame_b]->nodes.row(smp_physim->edges(smp_physim->linear(i), 1));
			}

			a = (1 - frame_alpha) * a + frame_alpha * a_b;
			b = (1 - frame_alpha) * b + frame_alpha * b_b;
		}

		// Get scissor springs from node data recorded for the given frame: rowwise from a to b has stiffness k
		void get_scissors_at_frame(Eigen::MatrixXd& a, Eigen::MatrixXd& b, int frame, double frame_alpha = 0) const
		{
			frame = std::min(frame, get_simulated_frames() - 1);
			Eigen::MatrixXd positions = get_positions_at_frame(frame, frame_alpha);
			a.resize(smp_physim->cross.rows() * 2, 3);
			b.resize(smp_physim->cross.rows() * 2, 3);
			for (int i = 0; i < smp_physim->cross.rows(); ++i)
			{
				a.row(i) = positions.row(smp_physim->edges(smp_physim->cross(i, 0), 0));
				b.row(i) = positions.row(smp_physim->edges(smp_physim->cross(i, 0), 1));
				a.row(i + smp_physim->cross.rows()) = positions.row(smp_physim->edges(smp_physim->cross(i, 1), 0));
				b.row(i + smp_physim->cross.rows()) = positions.row(smp_physim->edges(smp_physim->cross(i, 1), 1));
			}
		}

		// Get body vertices from data recorded for the given frame, interpolate with next frame if frame_alpha > 0
		auto get_body_vertices_at_frame(int frame, double frame_alpha = 0) const
		{
			frame = std::min(frame, get_simulated_frames() - 1);
			auto vertices = frame_info[frame]->body_vertices;
			if (0 < frame_alpha)
			{
				int frame_b = std::min(frame + 1, get_simulated_frames() - 1);
				if (frame_b != frame)
				{
					auto vertices_b = frame_info[frame_b]->body_vertices;
					vertices = (1 - frame_alpha) * vertices + frame_alpha * vertices_b;
				}
			}
			return vertices;
		}

		// Get body vertex normals from data recorded for the given frame, interpolate with next frame if frame_alpha > 0
		auto get_body_normals_at_frame(int frame, double frame_alpha = 0) const
		{
			frame = std::min(frame, get_simulated_frames() - 1);
			auto normals = frame_info[frame]->body_normals;
			if (0 < frame_alpha)
			{
				int frame_b = std::min(frame + 1, get_simulated_frames() - 1);
				if (frame_b != frame)
				{
					auto normals_b = frame_info[frame_b]->body_normals;
					normals = (1 - frame_alpha) * normals + frame_alpha * normals_b;
				}
			}
			return normals;
		}

		const Eigen::MatrixXi& get_body_faces_indexed() { return body_faces; }

		// Get rigid body mesh faces
		Eigen::MatrixXi get_body_faces(int flag = -1) const
		{
			Eigen::MatrixXi faces;
			if (flag == -1) faces = body_faces.leftCols(3);
			else
			{
				int nrows = (body_faces.rightCols(1).array() == flag).count();
				faces.resize(nrows, 3);
				int ind = 0;
				for (int i = 0; i < body_faces.rows(); ++i)
					if (body_faces(i, 3) == flag) faces.row(ind++) = body_faces.block(i, 0, 1, 3);
			}
			return faces;
		}

		auto get_bracket_vertices_at_frame(int frame, double frame_alpha = 0) const
		{
			frame = std::min(frame, get_simulated_frames() - 1);
			auto vertices = frame_info[frame]->bracket_vertices;
			if (0 < frame_alpha)
			{
				int frame_b = std::min(frame + 1, get_simulated_frames() - 1);
				if (frame_b != frame)
				{
					auto vertices_b = frame_info[frame_b]->bracket_vertices;
					vertices = (1 - frame_alpha) * vertices + frame_alpha * vertices_b;
				}
			}
			return vertices;
		}

		auto get_bracket_normals_at_frame(int frame, double frame_alpha = 0) const
		{
			frame = std::min(frame, get_simulated_frames() - 1);
			auto normals = frame_info[frame]->bracket_normals;
			if (0 < frame_alpha)
			{
				int frame_b = std::min(frame + 1, get_simulated_frames() - 1);
				if (frame_b != frame)
				{
					auto normals_b = frame_info[frame_b]->bracket_normals;
					normals = (1 - frame_alpha) * normals + frame_alpha * normals_b;
				}
			}
			return normals;
		}

		auto get_mem_normals_at_frame(int frame, double frame_alpha = 0) const
		{
			frame = std::min(frame, get_simulated_frames() - 1);
			auto normals = frame_info[frame]->mem_normals;
			if (0 < frame_alpha)
			{
				int frame_b = std::min(frame + 1, get_simulated_frames() - 1);
				if (frame_b != frame)
				{
					auto normals_b = frame_info[frame_b]->mem_normals;
					normals = (1 - frame_alpha) * normals + frame_alpha * normals_b;
				}
			}
			return normals;
		}

		inline Eigen::MatrixXi get_bracket_faces() const { return bracket_faces; }

		// Get bracket lengths for the given frame
		void get_bracket_length_at_frame(int actuator_id, int frame, Eigen::Vector4d& bracket_length) const
		{
			frame = std::min(frame, get_simulated_frames() - 1);
			for (int i = 0; i < 4; ++i)
			{
				Eigen::Vector2i edge = smp_physim->edges.row(smp_physim->linear(actuators[actuator_id]->segment_springs(i))).transpose();
				bracket_length(i) = (frame_info[frame]->nodes.row(edge(0)) - frame_info[frame]->nodes.row(edge(1))).norm();
			}
		}

		// Get bracket rest lengths for the given frame
		void get_bracket_restlen_at_frame(int actuator_id, int frame, Eigen::Vector4d& restlen) const
		{
			frame = std::min(frame, get_simulated_frames() - 1);
			for (int i = 0; i < 4; ++i)
				restlen(i) = frame_info[frame]->restlen(actuators[actuator_id]->segment_springs(i));
		}

		// Get bracket force magnitudes for the given frame
		void get_bracket_forces_at_frame(int actuator_id, int frame, Eigen::Vector4d& forces) const
		{
			frame = std::min(frame, get_simulated_frames() - 1);
			for (int i = 0; i < 4; ++i)
				forces(i) = frame_info[frame]->segment_forces(actuators[actuator_id]->segment_springs(i));
		}

		// Get actuator membrane energy recorded for the given frame
		double get_actuator_membrane_energy(int actuator_id, int frame) const
		{
			frame = std::min(frame, get_simulated_frames() - 1);
			return frame_info[frame]->linear_spring_energy(actuators[actuator_id]->membrane_edge);
		}

		// Get forces per base origin generated by linear springs mimicking membrane recorded for the given frame
		void get_mspring_forces_at_frame(Eigen::MatrixXd& vertices, Eigen::MatrixXd& forces, int frame)
		{
			frame = std::min(frame, get_simulated_frames() - 1);
			vertices.resize(smp_physim->base_orig_nodes.rows(), 3);
			for (int i = 0; i < smp_physim->base_orig_nodes.size(); ++i)
				vertices.row(i) = frame_info[frame]->nodes.row(smp_physim->base_orig_nodes(i));
			forces = frame_info[frame]->segment_forces;
		}

		// Get stored frame information by reference
		const std::vector<std::shared_ptr<FrameInfo>>& finfo() const { return frame_info; }

		~SmpSimulator()
		{
			shut_down_task();
		}

		bool save(filesystem::path path)
		{
			if (path.extension().string().compare(".bin") == 0)
			{
				std::ofstream ofs(path.string(), std::ios::binary);
				cereal::PortableBinaryOutputArchive ar(ofs);
				ar(cereal::make_nvp("simulation", *this));
			}
			else
			{
				std::ofstream ofs(path.string());
				cereal::JSONOutputArchive ar(ofs);
				ar(cereal::make_nvp("simulation", *this));
			}
			return true;
		}

		bool load(filesystem::path path)
		{
			if (path.extension().string().compare(".bin") == 0)
			{
				std::ifstream ifs(path.string(), std::ios::binary);
				cereal::PortableBinaryInputArchive ar(ifs);
				ar(cereal::make_nvp("simulation", *this));
			}
			else
			{
				std::ifstream ifs(path.string());
				cereal::JSONInputArchive ar(ifs);
				ar(cereal::make_nvp("simulation", *this));
			}
			return true;
		}

		template<typename Archive>
		void serialize(Archive& ar)
		{
			ar(cereal::make_nvp("name", name));
			ar(cereal::make_nvp("t_start", t_start));
			ar(cereal::make_nvp("t_end", t_end));
			ar(cereal::make_nvp("t_simulated", t_simulated));
			ar(cereal::make_nvp("t_step", t_step));
			ar(cereal::make_nvp("settings", settings));
			ar(cereal::make_nvp("smp", *smp));
			if (Archive::is_loading::value)
			{
				load_smp(*smp);
				frame_info.reserve(get_max_frames());
			}
			ar(cereal::make_nvp("frames", frame_info));
			if (Archive::is_loading::value)
			{
				for (auto finfo : frame_info)
				{
					smp_physim->get_frame_info(*finfo);
					compute_frame_info_derived(finfo.get());
				}
			}
		}

	};

}

#endif
