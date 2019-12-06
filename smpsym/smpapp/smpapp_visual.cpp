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

#include "smpapp_visual.h"

using namespace ci;
using namespace ci::app;

void SmpAppVisual::createPhongShader()
{
	try {
#if defined( CINDER_GL_ES )
		mPhongShader = gl::GlslProg::create(loadAsset("phong_es2.vert"), loadAsset("phong_es2.frag"));
#else
		mPhongShader = gl::GlslProg::create(loadAsset("phong.vert"), loadAsset("phong.frag"));
#endif
	}
	catch (Exception &exc) {
		CI_LOG_E("error loading phong shader: " << exc.what());
	}
}

void SmpAppVisual::createWireShader()
{
	try {
		mWireShader = gl::context()->getStockShader(gl::ShaderDef().color());
	}
	catch (Exception &exc) {
		CI_LOG_E("error loading wire shader: " << exc.what());
	}
}

void SmpAppVisual::createWireframeShader()
{
#if ! defined( CINDER_GL_ES )
	try {
		auto format = gl::GlslProg::Format()
			.vertex(loadAsset("wireframe.vert"))
			.geometry(loadAsset("wireframe.geom"))
			.fragment(loadAsset("wireframe.frag"));

		mWireframeShader = gl::GlslProg::create(format);
	}
	catch (Exception &exc) {
		CI_LOG_E("error loading wireframe shader: " << exc.what());
	}
#endif // ! defined( CINDER_GL_ES )
}

void SmpAppVisual::createSelectionShader()
{
	try {
		mSelectionShader = gl::GlslProg::create(loadAsset("select.vert"), loadAsset("select.frag"));
	}
	catch (Exception& exc) {
		CI_LOG_E("error loading selection shader: " << exc.what());
	}
}

void SmpAppVisual::update_grid()
{
	mGrid = gl::VertBatch::create(GL_LINES);
	mGrid->begin(GL_LINES);
	for (int i = -150; i <= 150; i += 10)
	{
		for (int j = 0; j < 4; ++j)
			mGrid->color(Color(0.25, 0.25, 0.25) * (i == 0 ? 2 : (i % 50 == 0 ? 1.5 : 1)));
		mGrid->vertex(i, -150, 0);
		mGrid->vertex(i, 150, 0);
		mGrid->vertex(-150, i, 0);
		mGrid->vertex(150, i, 0);
	}
	mGrid->end();
}

void SmpAppVisual::draw_grid()
{
	if (mShowGrid)
	{
		mGrid->draw();
		gl::color(0.3, 0.3, 0.3);
		gl::drawStrokedCircle(vec2(0, 0), 150);
	}
}
