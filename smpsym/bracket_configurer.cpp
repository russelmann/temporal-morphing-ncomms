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

#include "bracket_configurer.h"

#include "smp_simulator.h"

namespace smpup {

	void BracketConfigurer::single_async_task(void* user_params)
	{
		SmpOptions& options = *(SmpOptions*)user_params;

		int timesteps = 120;

		displacements.resize(lengths.size());
		displacements_weak.resize(lengths.size());

		SmpSettings settings;
		settings.mem_model = StaticSimulatorStruct::MEM_LINEAR;

		progress_length = 0;
		for (int length_id = 0; length_id < lengths.size(); ++length_id)
		{
			if (task_stop_requested()) break;

			displacements[length_id].resize(timesteps + 1, thicknesses[length_id].size());
			displacements_weak[length_id].resize(timesteps + 1, thicknesses[length_id].size());

			progress_thickness = 0;
			for (int thickness_id = 0; thickness_id < thicknesses[length_id].size(); ++thickness_id)
			{
				if (task_stop_requested()) break;

				Smp* smp = new Smp();
				smp->options = options;
				smp->options.min_bumper = 0;
				smp->options.min_gap = lengths(length_id) - smp->options.base_r * 2;
				smp->init_single();
				smp->spth.setConstant(thicknesses[length_id](thickness_id));

				for (int k = 0; k < 2; ++k)
				{
					settings.mimic_bounary_weakening = (k == 1);
					settings.timestep_mul = 1;
					std::vector<Eigen::MatrixXd>& disp = (settings.mimic_bounary_weakening) ? displacements_weak : displacements;
					SmpSimulator* simulator = new SmpSimulator(0, timesteps, timestep, settings, linear_mem_model, bracket_model);
					simulator->load_smp(*smp);
					simulator->simulate();
					while (simulator->get_simulated_frames() <= timesteps) process_sleep(10);

					for (int frame_id = 0; frame_id < simulator->get_simulated_frames(); ++frame_id)
					{
						Eigen::Vector4d bracket_length;
						simulator->get_bracket_length_at_frame(0, frame_id, bracket_length);
						disp[length_id](frame_id, thickness_id) = bracket_length.mean();
					}
				}
				++progress_thickness;
			}
			++progress_length;
		}
	}

	// TODO: this function rounds up to nearest length
	Eigen::Vector2d BracketConfigurer::compute_time_range(double initial_length, double target_length, bool weak) const
	{
		const std::vector<Eigen::MatrixXd>& disp = weak ? displacements_weak : displacements;
		Eigen::Vector2d time_range;

		int length_id = 0;
		while (length_id < lengths.size() - 1 && lengths(length_id + 1) <= initial_length) ++length_id;

		int time_id = 0;
		while (time_id < disp[length_id].rows() - 1 && target_length < disp[length_id](time_id + 1, 0))
			++time_id;

		time_range(0) = time_id * timestep;

		time_id = 0;
		while (time_id < disp[length_id].rows() - 1 && target_length < disp[length_id](time_id + 1, thicknesses[length_id].size() - 1))
			++time_id;

		time_range(1) = time_id * timestep;

		return time_range;
	}

	// N.B.: initial_lengths is for actuator; target_length is for brackets
	double BracketConfigurer::configure_thickness(double initial_length, double target_length, double time, bool weak) const
	{
		const std::vector<Eigen::MatrixXd>& disp = weak ? displacements_weak : displacements;

		int length_id = -1;
		while (length_id < lengths.size() - 1 && lengths(length_id + 1) <= initial_length) ++length_id;

		double thickness[2];
		for (int i = 0; i < 2; ++i)
		{
			Eigen::DenseIndex time_id = std::floor(time / timestep);
			time_id = std::min(time_id, disp[length_id + i].rows() - 1);

			if (length_id + i < 0 || lengths.size() <= length_id + i) continue;
			int thickness_id = -1;
			while (thickness_id < thicknesses[length_id + i].size() - 1 && disp[length_id + i](time_id, thickness_id + 1) <= target_length)
				++thickness_id;
			if (thickness_id < 0) {
				thickness[i] = thicknesses[length_id + i](0);
			}
			else if (thicknesses[length_id + i].size() <= thickness_id + 1) {
				thickness[i] = thicknesses[length_id + i].tail(1)(0);
			}
			else {
				double alpha = (disp[length_id + i](time_id, thickness_id + 1) - target_length) /
					(disp[length_id + i](time_id, thickness_id + 1) - disp[length_id + i](time_id, thickness_id));
				thickness[i] = alpha * thicknesses[length_id + i](thickness_id) + (1 - alpha) * thicknesses[length_id + i](thickness_id + 1);
			}
		}

		if (length_id < 0) return thickness[1];
		if (lengths.size() <= length_id + 1) return thickness[0];
		double beta = (lengths(length_id + 1) - initial_length) / (lengths(length_id + 1) - lengths(length_id));
		return beta * thickness[0] + (1 - beta) * thickness[1];
	}

}
