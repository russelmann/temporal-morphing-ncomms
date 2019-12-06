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

#ifndef SMPUP_VISUAL_H
#define SMPUP_VISUAL_H

#include "smpapp_visual.h"

#include "cinder/gl/gl.h"
#include "cinder/gl/Batch.h"
#include "cinder/gl/Context.h"
#include "cinder/gl/GlslProg.h"
#include "cinder/gl/Texture.h"
#include "cinder/gl/VboMesh.h"

#include "cinder/app/App.h"
#include "cinder/app/RendererGl.h"
#include "cinder/params/Params.h"
#include "cinder/Log.h"
#include "cinder/app/Platform.h"

#include "smp_addons.h"
#include "Eigen/Eigen"

class SmpupVisual
{
public:
	const SmpAppVisual* const smpapp_visual;

	std::string id;
	bool visible;                              // geomtery is visible
	//std::string name;                          // simulation name
	Eigen::MatrixXd time_ranges;               // actuator time ranges
	Color color;                               // visualization color
	gl::BatchRef stencil_act_batch;
	gl::BatchRef stencil_flt_batch;
	gl::BatchRef rigid_bodies_batch;
	gl::BatchRef rigid_bodies_phong_batch;
	gl::BatchRef bracket_batch;
	gl::BatchRef bracket_phong_batch;
	gl::BatchRef membrane_batch;
	gl::BatchRef membrane_phong_batch;

	Eigen::VectorXd time_landscape;
	gl::BatchRef time_landscape_batch;

	SmpupVisual(const SmpAppVisual* const smpapp_visual)
		: smpapp_visual(smpapp_visual)
	{
		static int id_count = 0;
		id = std::to_string(id_count);
		++id_count;
	}

	void draw_solid();
	void draw_stencil();
	void draw_membrane();
	void draw_time_landscape();
};

#endif