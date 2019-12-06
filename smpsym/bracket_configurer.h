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

#ifndef BRACKET_CONFIGURER_H
#define BRACKET_CONFIGURER_H

#include "interface.h"
#include "utils/self_serializer.h"
#include "utils/single_task_thread.h"
#include "Eigen/Core"

namespace smpup {

	class BracketConfigurer : public cereal::SelfSerializer<BracketConfigurer>, public SingleAsyncTask
	{
	private:
		std::shared_ptr<LinearMemModel> linear_mem_model;
		std::shared_ptr<BracketModel> bracket_model;

	public:
		double timestep;
		Eigen::VectorXd lengths;
		std::vector<Eigen::VectorXd> thicknesses;
		std::vector<Eigen::MatrixXd> displacements;
		std::vector<Eigen::MatrixXd> displacements_weak;

	private:

		int progress_length;
		int progress_thickness;
		void single_async_task(void* user_params = nullptr);

	public:
		
		BracketConfigurer(std::shared_ptr<LinearMemModel> linear_mem_model, std::shared_ptr<BracketModel> bracket_model, double timestep = 0.) :
			linear_mem_model(linear_mem_model), bracket_model(bracket_model), SelfSerializer("bracket_configurer"), timestep(timestep)
		{
			progress_length = 0;
			progress_thickness = 0;
		}

		// Perform sampling of single actuators deformation (assuming coherence of actuation)
		void single_actuator_sampling(const SmpOptions& options)
		{
			auto opt = new SmpOptions(options);
			run_task(opt);
		}

		inline float get_progress_length() { return float(progress_length) / lengths.size(); }
		inline float get_progress_thickness() { return (progress_length == lengths.size() ? (progress_length == 0 ? 0.0f : 1.0f) : float(progress_thickness) / thicknesses[progress_length].size()); }
		inline bool valid_computation_done() { return 0 < progress_length && progress_length == lengths.size(); }

		// TODO: this function rounds up to nearest length
		Eigen::Vector2d compute_time_range(double initial_length, double target_length, bool weak = false) const;

		// N.B.: initial_lengths is for actuator; target_length is for brackets
		double configure_thickness(double initial_length, double target_length, double time, bool weak = false) const;

	};

	template<class Archive>
	void serialize(Archive& ar, smpup::BracketConfigurer& bracket_configurer)
	{
		ar(cereal::make_nvp("timestep", bracket_configurer.timestep));
		ar(cereal::make_nvp("lengths", bracket_configurer.lengths));
		ar(cereal::make_nvp("thicknesses", bracket_configurer.thicknesses));
		ar(cereal::make_nvp("displacements", bracket_configurer.displacements));
		ar(cereal::make_nvp("displacements_weak", bracket_configurer.displacements_weak));
	}

}

#endif
