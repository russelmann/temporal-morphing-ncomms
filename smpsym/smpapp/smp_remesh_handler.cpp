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

#include "smp_remesh_handler.h"
#include "imgui_internal.h"

void SmpRemeshHandler::GuiWindowFlattened()
{
	static bool popen = true;
	if (!popen)
	{
		smpapp_visual->mShowSurfaceFlattened = false;
		popen = true;
	}
	ui::SetNextWindowSizeConstraints(vec2(10, 35), vec2(2000, 35));

	if (!ui::Begin("Flattened surface", &popen, ImGuiWindowFlags_NoCollapse)) { ui::End(); return; }

	mWinPos = ui::GetWindowPos();
	mWinWidth = ui::GetWindowWidth();

	ui::End();
}

void SmpRemeshHandler::update_ImportSurfacesBatch()
{
	TriMesh::Format fmt = TriMesh::Format().positions().colors();

	if (geo_remesher->has_imported_mesh())
	{
		Eigen::MatrixXd V;
		Eigen::MatrixXi F;
		geo_remesher->extract_imported_mesh(V, F);

		TriMeshRef tm = TriMesh::create(fmt);

		Eigen::Matrix<float, -1, -1, Eigen::RowMajor> Vr = V.cast<float>();
		tm->appendPositions((vec3*)Vr.data(), Vr.rows());
		Eigen::Matrix<float, -1, -1, Eigen::RowMajor> Cr(Vr.rows(), 3);
		Cr.col(0).array() = 0.3;
		Cr.col(1).array() = 0.6;
		Cr.col(2).array() = 0.9;
		tm->appendColors((Color*)Cr.data(), Cr.rows());
		for (int i = 0; i < F.rows(); ++i)
			tm->appendTriangle(F(i, 0), F(i, 1), F(i, 2));
		tm->recalculateNormals();
		//			tm->recalculateTangents();
		mImportedSurfaceBatch = gl::Batch::create(*tm, smpapp_visual->mWireframeShader);
	}
	if (geo_remesher->has_imported_mesh())
	{
		Eigen::MatrixXd V;
		Eigen::MatrixXi F;
		geo_remesher->extract_flattened_mesh(V, F);

		TriMeshRef tm = TriMesh::create(fmt);

		Eigen::Matrix<float, -1, -1, Eigen::RowMajor> Vr = V.cast<float>();
		tm->appendPositions((vec3*)Vr.data(), Vr.rows());
		Eigen::Matrix<float, -1, -1, Eigen::RowMajor> Cr(Vr.rows(), 3);
		Cr.col(0).array() = 0.3;
		Cr.col(1).array() = 0.6;
		Cr.col(2).array() = 0.9;
		tm->appendColors((Color*)Cr.data(), Cr.rows());
		for (int i = 0; i < F.rows(); ++i)
			tm->appendTriangle(F(i, 0), F(i, 1), F(i, 2));
		tm->recalculateNormals();
		//			tm->recalculateTangents();
		mFlattenedSurfaceBatch = gl::Batch::create(*tm, smpapp_visual->mWireframeShader);
	}
	if (geo_remesher->has_remeshed_mesh())
	{
		Eigen::MatrixXd V;
		Eigen::MatrixXi F;
		geo_remesher->extract_remeshed_mesh(V, F);

		TriMeshRef tm = TriMesh::create(fmt);

		Eigen::Matrix<float, -1, -1, Eigen::RowMajor> Vr = V.cast<float>();
		tm->appendPositions((vec3*)Vr.data(), Vr.rows());
		Eigen::Matrix<float, -1, -1, Eigen::RowMajor> Cr = Vr;
		Cr.col(0).array() = 0.3;
		Cr.col(1).array() = 0.8;
		Cr.col(2).array() = 0.9;
		tm->appendColors((Color*)Cr.data(), Cr.rows());
		for (int i = 0; i < F.rows(); ++i)
			tm->appendTriangle(F(i, 0), F(i, 1), F(i, 2));
		tm->recalculateNormals();
		//			tm->recalculateTangents();
		mRemeshedSurfaceBatch = gl::Batch::create(*tm, smpapp_visual->mWireframeShader);
	}
}

void SmpRemeshHandler::draw()
{
	smpapp_visual->mWireframeShader->uniform("uOpaqueness", 1.0f);
	gl::ScopedBlendAlpha blendScope;
	if (smpapp_visual->mShowSurfaceImported && mImportedSurfaceBatch)
	{
		gl::cullFace(GL_FRONT);
		smpapp_visual->mWireframeShader->uniform("uBrightness", 0.5f);
		mImportedSurfaceBatch->draw();
		gl::cullFace(GL_BACK);
		smpapp_visual->mWireframeShader->uniform("uBrightness", 0.8f);
		mImportedSurfaceBatch->draw();
	}
	//if (smpapp_visual->mShowSurfaceFlattened && mFlattenedSurfaceBatch)
	//{
	//	gl::cullFace(GL_FRONT);
	//	smpapp_visual->mWireframeShader->uniform("uBrightness", 0.5f);
	//	mFlattenedSurfaceBatch->draw();
	//	gl::cullFace(GL_BACK);
	//	smpapp_visual->mWireframeShader->uniform("uBrightness", 0.8f);
	//	mFlattenedSurfaceBatch->draw();
	//}
	if (smpapp_visual->mShowSurfaceRemeshed && mRemeshedSurfaceBatch)
	{
		gl::cullFace(GL_FRONT);
		smpapp_visual->mWireframeShader->uniform("uBrightness", 0.5f);
		mRemeshedSurfaceBatch->draw();
		gl::cullFace(GL_BACK);
		smpapp_visual->mWireframeShader->uniform("uBrightness", 0.8f);
		mRemeshedSurfaceBatch->draw();
	}
}

void SmpRemeshHandler::draw_2d()
{
	if (!smpapp_visual->mShowSurfaceFlattened) return;

	gl::setMatricesWindow(getWindowSize());
	
	gl::translate(mWinPos + vec2(0, 35));
	gl::scale(mWinWidth, mWinWidth);
	
	if (mFlattenedSurfaceBatch)
	{
		auto x_min = -1;
		auto x_max =  1;
		auto y_min = -1;
		auto y_max =  1;
		float x_scl = 1. / (x_max - x_min);
		float y_scl = 1. / (y_max - y_min);
		float scl = min(x_scl, y_scl);
	
		gl::color(0.3, 0.3, 0.3, 0.5);
		gl::cullFace(GL_BACK);
		gl::drawSolidRect(Rectf(0, 0, 1, x_scl / y_scl));
	
		gl::scale(1, -1);
		gl::translate(0, -1);
		gl::translate(vec2(-x_min, -y_min) * scl);
		gl::scale(scl, scl);
		smpapp_visual->mWireframeShader->uniform("uOpaqueness", 0.8f);
		smpapp_visual->mWireframeShader->uniform("uBrightness", 1.0f);
		mFlattenedSurfaceBatch->draw();
	}
}