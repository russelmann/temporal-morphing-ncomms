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

#ifndef SMP_IMGUI_H
#define SMP_IMGUI_H

#include "cinder/CinderImGui.h"
#include <string>

namespace ui = ImGui;

void tooltip(std::string tip);
void CheckboxEx(const char* label, bool* v, bool read_only);

#endif