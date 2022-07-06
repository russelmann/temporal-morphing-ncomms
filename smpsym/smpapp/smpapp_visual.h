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

#ifndef SMPAPP_VISUAL_H
#define SMPAPP_VISUAL_H

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

#include "cinder/CinderImGui.h"

using namespace ci;
using namespace ci::app;

class SmpAppVisual
{
public:
	enum ViewMode { SHADED, WIREFRAME };
	enum TexturingMode { NONE, PROCEDURAL, SAMPLER };

	ViewMode			mViewMode;
	TexturingMode		mTexturingMode;

	bool                mShowCoordinateFrame;
	bool                mShowGrid;
	bool                mShowStencilAct;
	bool                mShowStencilFlt;
	bool				mShowLinearSprings;
	bool				mShowScissorSprings;
	bool				mShowSolidPrimitive;
	bool                mShowAlign;
	bool                mShowMembrane;

	bool                mShowSurfaceImported;
	bool                mShowSurfaceFlattened;
	bool                mShowSurfaceRemeshed;

	bool                mShowTimeLandscape;
	bool                mShowTimeRanges;

	gl::GlslProgRef		mPhongShader;
	gl::GlslProgRef		mWireShader;
	gl::GlslProgRef		mWireframeShader;
	gl::GlslProgRef     mSelectionShader;

	gl::TextureRef		mTexture;

	gl::VertBatchRef    mGrid;

	SmpAppVisual()
	{
		mViewMode = WIREFRAME;
		mTexturingMode = PROCEDURAL;

		mShowCoordinateFrame = false;
		mShowGrid = false;
		mShowStencilAct = false;
		mShowStencilFlt = false;
		mShowLinearSprings = false;
		mShowScissorSprings = false;
		mShowSolidPrimitive = true;
		mShowAlign = false;
		mShowMembrane = true;

		mShowSurfaceImported = true;
		mShowSurfaceFlattened = false;
		mShowSurfaceRemeshed = true;

		mShowTimeLandscape = true;
		mShowTimeRanges = true;

		// Load the textures.
		gl::Texture::Format fmt;
		fmt.setAutoInternalFormat();
		fmt.setWrap(GL_REPEAT, GL_REPEAT);
		mTexture = gl::Texture::create(loadImage(loadAsset("texture.jpg")), fmt);

		// Load and compile the shaders.
		createPhongShader();
		createWireShader();
		createWireframeShader();
		createSelectionShader();
	}

	void createPhongShader();
	void createWireShader();
	void createWireframeShader();
	void createSelectionShader();

	void update_grid();
	void draw_grid();

};

#endif