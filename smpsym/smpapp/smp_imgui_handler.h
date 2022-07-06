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

#ifndef SMP_IMGUI_HANDLER_H
#define SMP_IMGUI_HANDLER_H

#include "interface.h"
#include "cinder/CinderImGui.h"

using namespace smpup;

// Structure for internal imgui handling
struct SmpGuiHandler
{
public:
	// Popups to run
	bool modal_clear_all_request;
	bool modal_core_settings_request;
	bool modal_new_single_request;
	bool modal_new_regular_request;
	bool model_open_stencil_request;
	bool modal_build_mesh_request;
	bool modal_build_video_request;
	bool modal_membrane_linearity_request;
	bool modal_build_specimen_block_request;
	bool modal_bracket_sampling_request;

	// Window visibility
	bool mShowSurfaceImport;
	bool mShowTimeLandscapeEditor;
	bool mShowCurrentSolverInfo;
	bool mShowConsole;
	bool mShowBracketExplorer;
	bool mShowPlottingWindow;
	bool mShowDataExportWindow;
	bool mShowDebugWindow;

	// Controls for options
	static void showGuiOptions(SmpOptions& options, int flags = 0);

	// Controls for options read-only
	template<class TypeOptions>
	static void showGuiOptions(const TypeOptions& options)
	{
		showGuiOptions(const_cast<TypeOptions&>(options), ImGuiInputTextFlags_ReadOnly);
	}

	// Controls for settings
	static void showGuiSettings(SmpSettings& settings, int flags = 0);

	// Controls for settings read-only
	template<class TypeSettings>
	static void showGuiSettings(const TypeSettings& settings)
	{
		showGuiSettings(const_cast<TypeSettings&>(settings), ImGuiInputTextFlags_ReadOnly);
	}

public:
	SmpGuiHandler();
};

#endif