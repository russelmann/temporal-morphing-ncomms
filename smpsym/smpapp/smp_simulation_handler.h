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

#ifndef SMP_SIMULATION_HANDLER
#define SMP_SIMULATION_HANDLER

#include "smpup_visual.h"
#include "smp_simulator.h"

using namespace smpup;
using namespace std;

class SimulationHandler : public SmpupVisual
{
public:
	std::shared_ptr<SmpSimulator> simulator;

	SimulationHandler(const SmpAppVisual* const smpapp_visual)
		: SmpupVisual(smpapp_visual)
	{ }

	void createStencils();
	void update_RigidBodyBatch(int current_frame, double current_frame_alpha);
	void update_MembraneBatch(int current_frame, double current_frame_alpha);
	void update_TimeLandscapeBatch();

	void draw_solid()
	{
		SmpupVisual::draw_solid();
	}
};

#endif