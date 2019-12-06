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

#include "smpup_visual.h"

void SmpupVisual::draw_solid()
{
	if (!visible) return;

	smpapp_visual->mWireframeShader->uniform("uOpaqueness", 0.65f);
	smpapp_visual->mPhongShader->uniform("uTexturingMode", SmpAppVisual::NONE);

	// Draw the primitive.
	if (smpapp_visual->mShowSolidPrimitive)
	{
		if (smpapp_visual->mViewMode == SmpAppVisual::WIREFRAME)
		{
			// We're using alpha blending, so render the back side first.
			gl::ScopedBlendAlpha blendScope;
			gl::ScopedFaceCulling cullScope(true, GL_FRONT);

			smpapp_visual->mWireframeShader->uniform("uBrightness", 0.5f);
			rigid_bodies_batch->draw();

			// Now render the front side.
			gl::cullFace(GL_BACK);

			smpapp_visual->mWireframeShader->uniform("uBrightness", 1.0f);
			rigid_bodies_batch->draw();
			bracket_batch->draw();
		}
		else
		{
			gl::ScopedBlendAlpha blendScope;
			gl::ScopedFaceCulling cullScope(true, GL_BACK);

			rigid_bodies_phong_batch->draw();
			bracket_phong_batch->draw();
		}
	}
}

void SmpupVisual::draw_stencil()
{
	if (!visible) return;

	if (smpapp_visual->mShowStencilAct)
	{
		gl::cullFace(GL_FRONT);
		smpapp_visual->mWireframeShader->uniform("uBrightness", 0.5f);
		stencil_act_batch->draw();
		gl::cullFace(GL_BACK);
		smpapp_visual->mWireframeShader->uniform("uBrightness", 1.0f);
		stencil_act_batch->draw();
	}
	if (smpapp_visual->mShowStencilFlt)
	{
		gl::cullFace(GL_FRONT);
		smpapp_visual->mWireframeShader->uniform("uBrightness", 0.5f);
		stencil_flt_batch->draw();
		gl::cullFace(GL_BACK);
		smpapp_visual->mWireframeShader->uniform("uBrightness", 1.0f);
		stencil_flt_batch->draw();
	}
}

void SmpupVisual::draw_membrane()
{
	if (!visible) return;

	if (smpapp_visual->mShowMembrane)
	{
		//learn how to use gl::ScopedFaceCulling();
		if (membrane_batch)
		{
			if (smpapp_visual->mViewMode == smpapp_visual->WIREFRAME)
			{
				smpapp_visual->mWireframeShader->uniform("uOpaqueness", 0.65f);

				gl::cullFace(GL_BACK);
				membrane_batch->draw();
				gl::cullFace(GL_FRONT);
				membrane_batch->draw();
			}
			else
			{
				gl::ScopedTextureBind scopedTextureBind(smpapp_visual->mTexture);
				smpapp_visual->mPhongShader->uniform("uTexturingMode", smpapp_visual->mTexturingMode);
				if (smpapp_visual->mTexturingMode == SmpAppVisual::PROCEDURAL) smpapp_visual->mPhongShader->uniform("uFreq", ivec2(1, 1));
				else if (smpapp_visual->mTexturingMode == SmpAppVisual::SAMPLER) smpapp_visual->mPhongShader->uniform("uFreq", ivec2(5, 5));

				gl::cullFace(GL_BACK);
				membrane_phong_batch->draw();
				gl::cullFace(GL_FRONT);
				membrane_phong_batch->draw();
			}
		}
	}
}

void SmpupVisual::draw_time_landscape()
{
	if (!visible) return;

	if (smpapp_visual->mShowTimeLandscape)
	{
		smpapp_visual->mWireframeShader->uniform("uBrightness", 1.0f);
		gl::cullFace(GL_BACK);
		time_landscape_batch->draw();
		gl::cullFace(GL_FRONT);
		time_landscape_batch->draw();
	}
}