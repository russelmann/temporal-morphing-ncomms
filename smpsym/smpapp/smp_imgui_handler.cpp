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

#include "smp_imgui_handler.h"
#include "smp_imgui.h"

namespace ui = ImGui;

SmpGuiHandler::SmpGuiHandler()
{
	modal_clear_all_request = false;
	modal_core_settings_request = false;
	modal_new_single_request = false;
	modal_new_regular_request = false;
	model_open_stencil_request = false;
	modal_build_mesh_request = false;
	modal_build_video_request = false;
	modal_membrane_linearity_request = false;
	modal_build_specimen_block_request = false;
	modal_bracket_sampling_request = false;

	mShowSurfaceImport = false;
	mShowTimeLandscapeEditor = false;
	mShowCurrentSolverInfo = false;
	mShowConsole = false;
	mShowBracketExplorer = false;
	mShowPlottingWindow = false;
	mShowDataExportWindow = false;
	mShowDebugWindow = false;
}

void SmpGuiHandler::showGuiOptions(SmpOptions& options, int flags)
{
	bool read_only = flags & ImGuiInputTextFlags_ReadOnly;
	ui::PushItemWidth(60);
	if (ui::CollapsingHeader("Structure parameters", ImGuiTreeNodeFlags_DefaultOpen))
	{
		if (ui::Button("Save as default")) options.save();
		if (!read_only)
		{
			ui::SameLine();
			if (ui::Button("Restore default")) options.load();
		}
		ui::Text("Printed structure parameters");
		ui::InputFloat("Base radius", &options.base_r, 0, 0, -1, flags);
		ui::InputFloat("Base thickness", &options.base_th, 0, 0, -1, flags);
		ui::InputFloat("Base indent", &options.base_p, 0, 0, -1, flags);
		ui::InputFloat("Bracket width", &options.bracket_w, 0, 0, -1, flags);
		ui::InputFloat("Minimal gap", &options.min_gap, 0, 0, -1, flags);
		ui::InputFloat("Minimal bumper", &options.min_bumper, 0, 0, -1, flags);
		ui::Spacing();
		ui::Text("Membrane parameters");
		ui::InputFloat("Strech factor", &options.tau, 0, 0, -1, flags);
		ui::InputFloat("Shear modulus", &options.lame2, 0, 0, -1, flags);
		ui::InputFloat("Thickness", &options.mem_h, 0, 0, -1, flags);
		ui::Spacing();
		ui::Text("Other parameters");
		ui::InputFloat("Scissor stiffness", &options.scissor, 0, 0, -1, flags);
		if (ui::TreeNode("Default print options"))
		{
			CheckboxEx("Print markers", &options.default_print_options.markers, read_only);
			ui::InputFloat("Align diameter", &options.default_print_options.align_diam, 0, 0, -1, flags); tooltip("Diameter of alignment holes");
			ui::InputInt("Handle", &options.default_print_options.handle, 0, 0, flags); tooltip("Index of the base where handle must be attached");
			CheckboxEx("Print alignment holes", &options.default_print_options.alignment, read_only);
			ui::TreePop();
		}
	}
}

void SmpGuiHandler::showGuiSettings(SmpSettings& settings, int flags)
{
	bool read_only = flags & ImGuiInputTextFlags_ReadOnly;
	ui::PushItemWidth(60);
	if (ui::CollapsingHeader("Simulator parameters", ImGuiTreeNodeFlags_DefaultOpen))
	{
		if (read_only)
			ui::InputText("Membrane model", &StaticSimulatorStruct::MEM_MODEL_NAME(settings.mem_model), flags);
		else
			ui::Combo("Membrane model", (int*)&settings.mem_model, StaticSimulatorStruct::MEM_MODEL_COMBO());

		if (settings.mem_model == StaticSimulatorStruct::MEM_FINE)
		{
			ui::InputText("Triangle", &settings.triangle, flags);
			ui::InputInt("Open edge vertices", &settings.vert_per_edge, 0, 0, flags); tooltip("Numnber of vertices to split membanre open boundary edges");
		}
		else if (settings.mem_model == StaticSimulatorStruct::MEM_LINEAR)
		{
			CheckboxEx("Minic boundary weakening", &settings.mimic_bounary_weakening, read_only); tooltip("Reduce strength of segment springs approximating membrane for boundary actuators");
		}
		else if (settings.mem_model == StaticSimulatorStruct::MEM_CONSTANT)
		{
			ui::InputFloat("Force magnitude", &settings.constant_force, 0, 0, -1, flags);
			if (!read_only)
			{
				for (int i = 1; i < 6; ++i)
				{
					ui::SameLine();
					if (ui::Button(std::to_string(i).c_str())) settings.constant_force = i;
				}
			}
		}
		ui::InputFloat("Boundary weakening", &settings.boundary_weakening, 0, 0, -1, flags);
		ui::InputFloat("Plastic fraction", &settings.plastic_fraction, 0, 0, -1, flags);
		ui::InputInt("Time step multiplier", &settings.timestep_mul, 0, 0, flags);
		if (!read_only) settings.timestep_mul = std::max(settings.timestep_mul, 1);

		if (ui::TreeNodeEx("Solver settings", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ui::InputFloat("Tolerance", &settings.simulator_settings.tol, 0, 0, -1, flags);
			ui::InputInt("Max iterations dry", &settings.simulator_settings.max_iter_dry, 0, 0, flags);
			ui::InputInt("Max iterations wet", &settings.simulator_settings.max_iter_wet, 0, 0, flags);
			ui::InputFloat("Max step size", &settings.simulator_settings.max_step_size, 0, 0, -1, flags);
			CheckboxEx("Freeze centroid", &settings.simulator_settings.freeze_centroid, read_only);
			CheckboxEx("Compute condition number", &settings.simulator_settings.compute_condition, read_only);
			ui::TreePop();
		}
	}
}
