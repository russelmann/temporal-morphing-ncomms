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

#ifndef SMP_REMESH_HANDLER
#define SMP_REMESH_HANDLER

#include "smpapp_visual.h"
#include "utils/smp_geogram.h"

using namespace smpup;
using namespace std;

// This class combines Handler and Visual since it is small
class SmpRemeshHandler
{
public:
	SmpAppVisual* smpapp_visual;

	shared_ptr<GeoRemesher> geo_remesher;

	gl::BatchRef        mImportedSurfaceBatch;
	gl::BatchRef        mFlattenedSurfaceBatch;
	gl::BatchRef        mRemeshedSurfaceBatch;

	vec2                mWinPos;
	float               mWinWidth;

	SmpRemeshHandler(SmpAppVisual* smpapp_visual)
		: smpapp_visual(smpapp_visual)
	{
		geo_remesher = make_shared<GeoRemesher>();
	}

	void GuiWindowFlattened();

	void draw();

	void draw_2d();

	void update_ImportSurfacesBatch();
};

#endif