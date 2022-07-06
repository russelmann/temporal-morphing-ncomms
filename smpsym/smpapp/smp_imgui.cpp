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

#include "smpapp.h"

#include "igl/writeOBJ.h"

#include "imgui/imgui_internal.h"

using namespace ci;
using namespace ci::app;
using namespace std;
using namespace smpup;

namespace ui = ImGui;

void tooltip(string tip)
{
	if (ImGui::IsItemHovered())
	{
		ImGui::BeginTooltip();
		ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
		ImGui::TextUnformatted(tip.c_str());
		ImGui::PopTextWrapPos();
		ImGui::EndTooltip();
	}
}

void CheckboxEx(const char* label, bool* v, bool read_only)
{
	if (read_only) ui::PushItemFlag(ImGuiItemFlags_Disabled, true);
	ui::Checkbox(label, v);
	if (read_only) ui::PopItemFlag();
}

void SmpApp::GuiPopups()
{
	PopupClearAll();
	PopupNewSingle();
	PopupNewRegular();
	PopupOpenStencil();
	PopupCloseApp();
	PopupBuildMesh();
	PopupBuildVideo();
	PopupCoreSettings();
	PopupMembraneLinearity();
	PopupBuildSpecimenBlock();
	PopupBracketSampling();
}

void SmpApp::showGui()
{
	GuiMainMenu();
	GuiPopups();
	GuiMainWindow();
	if (mShowSurfaceImport) GuiSurfaceImport();
	if (mShowTimeLandscapeEditor) GuiLandscapeEditor();
	if (mShowCurrentSolverInfo) GuiSolverInfo();
	if (mShowConsole) GuiConsole();
	if (mShowBracketExplorer) GuiBracketExplorer();
	if (mShowPlottingWindow) GuiPlottingWindow();
	if (mShowDataExportWindow) GuiDataExportWindow();
	if (mShowDebugWindow) GuiDebugWindow();
#ifdef USE_GEOGRAM
	if (mShowSurfaceFlattened) mSmpRemeshHandler.GuiWindowFlattened();
#endif
	
	//ui::ShowDemoWindow();
}
