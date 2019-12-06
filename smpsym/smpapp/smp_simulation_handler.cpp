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

#include "smp_simulation_handler.h"
#include "igl/per_vertex_normals.h"

void update_vbo(gl::VboMeshRef vbo, Eigen::Matrix<float, -1, -1, Eigen::RowMajor>& V, Eigen::Matrix<float, -1, -1, Eigen::RowMajor>& N, Eigen::MatrixXi& F, Color color)
{
	auto vbo_map_pos = vbo->mapAttrib3f(geom::POSITION, false);
	for (int i = 0; i < V.rows(); ++i)
		vbo_map_pos[i] = vec3(V(i, 0), V(i, 1), V(i, 2));
	vbo_map_pos.unmap();

	auto vbo_map_norm = vbo->mapAttrib3f(geom::NORMAL, false);
	for (int i = 0; i < V.rows(); ++i)
		vbo_map_norm[i] = vec3(N(i, 0), N(i, 1), N(i, 2));
	vbo_map_norm.unmap();

	auto vbo_map_col = vbo->mapAttrib3f(geom::COLOR, false);
	for (int i = 0; i < V.rows(); ++i)
		vbo_map_col[i] = color;
	
	vbo_map_col.unmap();
}

void SimulationHandler::createStencils()
{
	const SmpData* smp = simulator->get_smp();
	Eigen::VectorXd sc = (3 - smp->eu.array()) / 6;

	TriMesh::Format fmt = TriMesh::Format().positions().colors();
	TriMeshRef tma = TriMesh::create(fmt);
	if (smp->uv.size() > 0)
	{
		for (int i = 0; i < smp->uv.rows(); ++i)
		{
			tma->appendPosition(vec3(smp->V(i, 0), smp->V(i, 1), smp->V(i, 2)));
			tma->appendColorRgb(Color(ColorModel::CM_HSV, max(0., sc(i)), 1, 1));
		}
		for (int i = 0; i < smp->F.rows(); ++i)
		{
			tma->appendTriangle(smp->F(i, 0), smp->F(i, 1), smp->F(i, 2));
		}
	}
	stencil_act_batch = gl::Batch::create(*tma.get(), smpapp_visual->mWireframeShader);

	TriMeshRef tmf = TriMesh::create(fmt);
	if (smp->uv.size() > 0)
	{
		for (int i = 0; i < smp->uv.rows(); ++i)
		{
			tmf->appendPosition(vec3(smp->uv(i, 0), smp->uv(i, 1), 1e-1)); // small offset from membrane
			tmf->appendColorRgb(Color(ColorModel::CM_HSV, max(0., sc(i)), 1, 1));
		}
		for (int i = 0; i < smp->F.rows(); ++i)
		{
			tmf->appendTriangle(smp->F(i, 0), smp->F(i, 1), smp->F(i, 2));
		}
	}
	stencil_flt_batch = gl::Batch::create(*tmf.get(), smpapp_visual->mWireframeShader);
}


void SimulationHandler::update_RigidBodyBatch(int current_frame, double current_frame_alpha)
{
	TriMesh::Format fmt = TriMesh::Format().positions().normals().colors();
	auto V = simulator->get_body_vertices_at_frame(current_frame, current_frame_alpha);
	auto N = simulator->get_body_normals_at_frame(current_frame, current_frame_alpha);
	auto F = simulator->get_body_faces();

	if (!rigid_bodies_batch)
	{
		TriMeshRef tm = TriMesh::create(fmt);
		tm->appendPositions(vector<vec3>(V.rows()).data(), V.rows());
		tm->appendNormals(vector<vec3>(N.rows()).data(), N.rows());
		tm->appendColors(vector<Color>(V.rows()).data(), V.rows());
		for (int i = 0; i < F.rows(); ++i)
			tm->appendTriangle(F(i, 0), F(i, 1), F(i, 2));
		gl::VboMeshRef vbo_ref = gl::VboMesh::create(*tm);
		rigid_bodies_batch = gl::Batch::create(vbo_ref, smpapp_visual->mWireframeShader);
		rigid_bodies_phong_batch = gl::Batch::create(vbo_ref, smpapp_visual->mPhongShader);
	}
	update_vbo(rigid_bodies_batch->getVboMesh(), V, N, F, color);

	{
		auto V = simulator->get_bracket_vertices_at_frame(current_frame, current_frame_alpha);
		auto N = simulator->get_bracket_normals_at_frame(current_frame, current_frame_alpha);
		auto F = simulator->get_bracket_faces();

		if (!bracket_batch)
		{
			TriMeshRef tm = TriMesh::create(fmt);
			tm->appendPositions(vector<vec3>(V.rows()).data(), V.rows());
			tm->appendNormals(vector<vec3>(N.rows()).data(), N.rows());
			tm->appendColors(vector<Color>(V.rows()).data(), V.rows());
			for (int i = 0; i < F.rows(); ++i)
				tm->appendTriangle(F(i, 0), F(i, 1), F(i, 2));
			gl::VboMeshRef vbo_ref = gl::VboMesh::create(*tm);
			bracket_batch = gl::Batch::create(vbo_ref, smpapp_visual->mWireframeShader);
			bracket_phong_batch = gl::Batch::create(vbo_ref, smpapp_visual->mPhongShader);

		}
		update_vbo(bracket_batch->getVboMesh(), V, N, F, color);
	}
}

void SimulationHandler::update_MembraneBatch(int current_frame, double current_frame_alpha)
{
	Eigen::MatrixXi Fm = simulator->psim()->memtri;
	if (Fm.rows() == 0)
	{
		membrane_batch.reset();
		return;
	}

	Eigen::MatrixXd Vm = simulator->get_positions_at_frame(current_frame, current_frame_alpha);
	auto N = simulator->get_mem_normals_at_frame(current_frame, current_frame_alpha);

	Eigen::MatrixXd Vm_tex = simulator->get_positions_at_frame(0);
	Vm_tex /= simulator->get_smp()->options.tau;

	TriMesh::Format fmt = TriMesh::Format().positions().normals().colors(4).texCoords();
	TriMeshRef mMembraneTriMesh = TriMesh::create(fmt);

	for (int i = 0; i < Vm.rows(); ++i)
	{
		glm::vec3 pos(Vm(i, 0), Vm(i, 1), Vm(i, 2));
		mMembraneTriMesh->appendPosition(pos);
		mMembraneTriMesh->appendNormal(vec3(N(i, 0), N(i, 1), N(i, 2)));
		mMembraneTriMesh->appendColorRgba(ColorA(0.9, 0.9, 0, 0.5));
		mMembraneTriMesh->appendTexCoord(vec2(Vm_tex(i, 0), Vm_tex(i, 1)));
	}
	//mMembraneTriMesh->appendNormals(vector<vec3>(N.rows()).data(), N.rows());

	for (int i = 0; i < Fm.rows(); ++i)
		mMembraneTriMesh->appendTriangle(Fm(i, 0), Fm(i, 1), Fm(i, 2));

	if (smpapp_visual->mWireframeShader) membrane_batch = gl::Batch::create(*mMembraneTriMesh.get(), smpapp_visual->mWireframeShader);
	if (smpapp_visual->mPhongShader) membrane_phong_batch = gl::Batch::create(*mMembraneTriMesh.get(), smpapp_visual->mPhongShader);
}

void SimulationHandler::update_TimeLandscapeBatch()
{
	const SmpData* smp = simulator->get_smp();
	TriMesh::Format fmt = TriMesh::Format().positions().colors();
	TriMeshRef tmt = TriMesh::create(fmt);
	if (smp->uv.size() > 0)
	{
		for (int i = 0; i < smp->uv.rows(); ++i)
		{
			tmt->appendPosition(vec3(smp->uv(i, 0), smp->uv(i, 1), time_landscape(i)));
			tmt->appendColorRgb(Color(ColorModel::CM_HSV, time_landscape(i) / 120, 1, 1)); //TODO: max time 120 is hard coded
		}
		for (int i = 0; i < smp->F.rows(); ++i)
			tmt->appendTriangle(smp->F(i, 0), smp->F(i, 1), smp->F(i, 2));
	}
	time_landscape_batch = gl::Batch::create(*tmt.get(), smpapp_visual->mWireframeShader);
}