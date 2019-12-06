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
#include "smp_imgui.h"
#include "igl/writeOBJ.h"

void SmpApp::GuiMainWindow()
{
	if (!ui::Begin("SmpApp")) { ui::End(); return; }

	ui::Text("FPS %.f", getAverageFps());

	if (ui::CollapsingHeader("Visualization", ImGuiTreeNodeFlags_DefaultOpen))
	{
		ui::PushItemWidth(100);

		static int u_camera_mode = mCameraActive == &mCameraPersp ? 0 : 1;
		if (ui::Combo("Camera mode", &u_camera_mode, "Perspective\0Orthotropic\0"))
		{
			if (u_camera_mode == 0) mCameraActive = &mCameraPersp;
			else mCameraActive = &mCameraOrtho;
		}
		ui::Combo("View mode", (int*)&mViewMode, "Shaded\0Wireframe\0");
		ui::Combo("Texturing mode", (int*)&mTexturingMode, "None\0Procedural\0Sampler\0");

		ui::Spacing();

		ui::Checkbox("Show actuated stencil", &mShowStencilAct);
		ui::Checkbox("Show flat stencil", &mShowStencilFlt);

		ui::Spacing();

		ui::Checkbox("Show segment springs", &mShowLinearSprings);
		ui::Checkbox("Show scissor springs", &mShowScissorSprings);
		ui::Checkbox("Show solids", &mShowSolidPrimitive);
		ui::Checkbox("Show alignment", &mShowAlign);
		ui::Checkbox("Show membrane", &mShowMembrane);
	}

	if (ui::CollapsingHeader("Navigation"))
	{
		if (ui::Button("Look at origin")) mCameraPersp.lookAt(vec3(0));

		if (!mSimulationHandlers.empty())
		{
			if (ui::Button("Zoom selected actuator") && mPickedActuator >= 0)
			{
				//TODO: clean up dirty implementation
				Eigen::MatrixXd V = mSelectedSimHandler->simulator->get_positions_at_frame(mCurrentFrame); // not accounting smooth playback
				Eigen::MatrixXi Fa = mSelectedSimHandler->simulator->get_actuator_faces();
				vec3 pos;
				for (int i = 0; i < Fa.rows(); ++i)
					if (Fa(i, 3) == mPickedActuator) pos = vec3(V(Fa(i, 0), 0), V(Fa(i, 0), 1), V(Fa(i, 0), 2));
				mCameraPersp.lookAt(pos);
				vec3 dir = pos - mCameraPersp.getEyePoint();
				if (length(dir) > 100) mCameraPersp.setEyePoint(mCameraPersp.getEyePoint() + normalize(dir) * (length(dir) - 100));
				mCameraPersp.setPivotDistance(min(float(100), length(dir)));
			}
			if (ui::Button("Select max scale 3d")) mPickedActuator = mSelectedSimHandler->simulator->get_smp()->act_max_scale;
			if (ui::Button("Select max length"))
				mSelectedSimHandler->simulator->get_actuator_lengths_at_frame(mCurrentFrame).maxCoeff(&mPickedActuator);
			if (ui::Button("Select min length"))
				mSelectedSimHandler->simulator->get_actuator_lengths_at_frame(mCurrentFrame).minCoeff(&mPickedActuator);
		}

		ui::Checkbox("Show coordinate frame", &mShowCoordinateFrame);
		ui::Checkbox("Show grid", &mShowGrid);
		ui::Checkbox("Rotation", &mRotation);
	}

	if (!mSimulationHandlers.empty())
	{
		//if (ui::CollapsingHeader("Info", ImGuiTreeNodeFlags_DefaultOpen))
		if (ui::CollapsingHeader("Info"))
		{
			ui::PushItemWidth(40);

			static int u_nbases;
			u_nbases = mSelectedSimHandler->simulator->get_smp()->F.rows();
			ui::InputInt("bases", &u_nbases, 0, 0, ImGuiInputTextFlags_ReadOnly);

			static int u_nactuators;
			u_nactuators = mSelectedSimHandler->simulator->get_smp()->FF.rows();
			ui::InputInt("actuators", &u_nactuators, 0, 0, ImGuiInputTextFlags_ReadOnly);

			ui::PopItemWidth();
		}

		if (ui::CollapsingHeader("Simulation", ImGuiTreeNodeFlags_DefaultOpen))
		{
			static shared_ptr<SimulationHandler> shandler_erase;
			for (auto shandler_iter = mSimulationHandlers.begin(); shandler_iter != mSimulationHandlers.end(); ++shandler_iter)
			{
				auto shandler = *shandler_iter;

				//static bool u_show;
				ui::Checkbox(("##show_sim_" + shandler->simulator->get_name()).c_str(), &shandler->visible);

				ui::SameLine(0, 2);
				if (ui::Button(("X##del_sim_" + shandler->simulator->get_name()).c_str())) shandler_erase = shandler;

				ui::SameLine(0, 2);
				ui::ColorEdit3(("MyColor##" + shandler->simulator->get_name()).c_str(),
					(float*)&shandler->color,
					ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel);

				ui::SameLine(0, 2);
				if (ui::Button(("*##set_sim_" + shandler->simulator->get_name()).c_str()))
					ui::OpenPopup(("Simulation properties##sim" + shandler->id).c_str());
				if (ui::BeginPopupModal(("Simulation properties##sim" + shandler->id).c_str(), NULL, ImGuiWindowFlags_AlwaysAutoResize))
				{
					ui::InputText("Name", &shandler->simulator->get_name());
					showGuiOptions(shandler->simulator->get_smp()->options);
					showGuiSettings(shandler->simulator->get_settings());
					if (ui::Button("Ok")) ui::CloseCurrentPopup();
					ui::EndPopup();
				}

				ui::SameLine(0, 7);
				if (ui::Selectable(shandler->simulator->get_name().c_str(), mSelectedSimHandler == shandler, 0, ImVec2(ui::GetContentRegionAvailWidth() - 30, 0)))
					mSelectedSimHandler = shandler;
				if (ui::BeginPopupContextItem(("move context menu##" + shandler->simulator->get_name()).c_str()))
				{
					if (ui::Selectable("move top") && shandler_iter != mSimulationHandlers.begin()) mSimulationHandlers.splice(mSimulationHandlers.begin(), mSimulationHandlers, next(--shandler_iter));
					if (ui::Selectable("move up") && shandler_iter != mSimulationHandlers.begin()) swap(*(prev(shandler_iter)), *shandler_iter);
					if (ui::Selectable("move down") && next(shandler_iter) != mSimulationHandlers.end()) swap(*(next(shandler_iter)), *shandler_iter);
					if (ui::Selectable("move bottom") && next(shandler_iter) != mSimulationHandlers.end()) mSimulationHandlers.splice(mSimulationHandlers.end(), mSimulationHandlers, prev(++shandler_iter));
					ui::Spacing();
					if (ui::Selectable("remove...")) shandler_erase = shandler;
					ImGui::EndPopup();
				}

				ui::SameLine(ui::GetWindowContentRegionWidth() - 15);
				if (ui::TreeNode(("##sim_detail_" + shandler->simulator->get_name()).c_str()))
				{
					static string u_smp_name;
					u_smp_name = shandler->simulator->get_smp()->name;
					ui::InputText("Smp name", &u_smp_name, ImGuiInputTextFlags_ReadOnly);
					static string u_smp_dir;
					u_smp_dir = shandler->simulator->get_smp()->working_dir;
					ui::InputText("##smp_dir", &u_smp_dir, ImGuiInputTextFlags_ReadOnly);
					ui::SameLine();
					if (ui::Button("Working dir")) os_reveal_folder(u_smp_dir);
					ui::PushItemWidth(50);
					ui::InputText(("Membrane model##" + shandler->simulator->get_name()).c_str(), &shandler->simulator->get_membrane_model_name(),
						ImGuiInputTextFlags_ReadOnly);
					if (shandler->simulator->get_membrane_model() == StaticSimulatorStruct::MEM_FINE)
					{
						string u_triang = shandler->simulator->get_settings().triangle;
						ui::InputText(("Triangle##" + shandler->simulator->get_name()).c_str(), &u_triang, ImGuiInputTextFlags_ReadOnly);

						static int u_nmemtri;
						u_nmemtri = mSelectedSimHandler->simulator->psim()->memtri.rows();
						ui::InputInt("Membrane triangles", &u_nmemtri, 0, 0, ImGuiInputTextFlags_ReadOnly);
					}
					ui::PopItemWidth();
					ui::TreePop();
				}
			}

			ui::PushItemWidth(-1);

			static string u_sim_name;
			static SmpSettings u_settings;
			if (ui::Button("add simulation...", ImVec2(-1, 0)))
			{
				u_settings = SmpSettings();
				u_sim_name = "Membrane new";
				ui::OpenPopup("Add simulation");
			}

			if (ui::BeginPopupModal("Add simulation", nullptr, ImGuiWindowFlags_AlwaysAutoResize))
			{
				ui::ColorEdit3("Color", u_settings.color.data(), ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel);
				ui::SameLine();
				ui::InputText("Simulation name", &u_sim_name);

				showGuiOptions(mSmpOptions);
				showGuiSettings(u_settings);

				if (ui::Button("Ok"))
				{
					auto shandler = mSelectedSimHandler;
					Smp smp = *mSelectedSimHandler->simulator->get_smp();
					smp.options = mSmpOptions; //TODO: this is very bad; make function set_options to update variables of SmpMesh and else...
					addSimulationHandler(u_sim_name, u_settings, smp);
					mSelectedSimHandler->time_landscape = shandler->time_landscape;
					mSelectedSimHandler->update_TimeLandscapeBatch();
					ui::CloseCurrentPopup();
				}
				ui::SameLine();
				if (ui::Button("Cancel")) ui::CloseCurrentPopup();
				ui::EndPopup();
			}

			ui::Separator();

			if (mPlayback) {
				if (ui::Button("PAUSE", vec2(50, 0))) mPlayback = false;
			}
			else {
				if (ui::Button("PLAY", vec2(50, 0))) mPlayback = true;
			}

			ui::PushItemWidth(60);
			ui::SameLine();
			static int u_play_speed = 3;
			if (ui::Combo("##combo_play_speed", &u_play_speed, "x8\0x4\0x2\0real\0x0.5"))
			{
				vector<double> playspeeds{ 8, 4, 2, 1, 0.5 };
				mPlaySpeed = playspeeds[u_play_speed];
			}
			ui::PopItemWidth();
			ui::PushItemWidth(40);
			float nTime = (mCurrentFrame - 1 + mCurrentFrameAlpha) * mSelectedSimHandler->simulator->get_timestep();
			ui::SameLine();
			ui::Checkbox("Loop", &mPlayLoop);
			ui::SameLine();
			ui::Checkbox("Smooth", &mPlaySmooth);

			ui::Checkbox("Lock init", &mLockRestConfiguration); tooltip("Exclude initial frame from output");
			ui::SameLine();
			ui::InputFloat("now", &nTime, 0, 0, 1, ImGuiInputTextFlags_ReadOnly);
			ui::PushItemWidth(-1);
			ui::SliderInt("", &mCurrentFrame, 0, mSelectedSimHandler->simulator->get_total_frames());
			if (mLockRestConfiguration) mCurrentFrame = max(mCurrentFrame, 1);
			mCurrentFrame = min(mCurrentFrame, mSelectedSimHandler->simulator->get_simulated_frames() - 1);

			static Eigen::VectorXf u_cstrain;
			static double u_cstrain_max = 1;
			u_cstrain = -Eigen::VectorXf::Ones(mSelectedSimHandler->simulator->get_total_frames());
			for (int i = 0; i < mSelectedSimHandler->simulator->get_simulated_frames(); ++i)
				u_cstrain[i] = fabs(mSelectedSimHandler->simulator->get_cross_strain_at_frame(i));
			u_cstrain_max = u_cstrain.maxCoeff();
			if (u_cstrain_max > 0) u_cstrain /= u_cstrain_max;
			ui::PlotLines("", u_cstrain.data(), u_cstrain.size(), 0, to_string(u_cstrain_max).c_str(), 0, 1, ImVec2(0, 50));

			ui::Text("Computed time");
			char buf[32];
			sprintf(buf, "%.1f/%g", mSelectedSimHandler->simulator->get_time_simulated(), mSelectedSimHandler->simulator->get_time_end());
			ui::ProgressBar(mSelectedSimHandler->simulator->get_time_simulated() / mSelectedSimHandler->simulator->get_time_end(), ImVec2(-1, 0), buf);
			
			if (ui::Button("Time editor")) mShowTimeLandscapeEditor = true;
			ui::SameLine();
			if (mSelectedSimHandler->simulator->is_running_task())
			{
				if (ui::Button("Stop simulation")) mSelectedSimHandler->simulator->stop_simulation();
			}
			else
			{
				if (ui::Button("Simulate")) mSelectedSimHandler->simulator->simulate();
			}
			if (mSelectedSimHandler->simulator->get_smp()->time_landscape(0) == 0.)
			{
				ui::SameLine();
				if (ui::Button("Default config")) mSelectedSimHandler->simulator.reset(mSelectedSimHandler->simulator->configure_brackets(mSelectedSimHandler->time_landscape).release());
			}

			if (ui::TreeNodeEx("Summary"))
			{
				double total_time = 0;
				for (int i = 1; i < mSelectedSimHandler->simulator->get_simulated_frames(); ++i)
					total_time += mSelectedSimHandler->simulator->finfo()[i]->solver_info->time;
				ui::Text("Total simulation time: %f", total_time);

				ui::Text("Membrane triangles: %d", mSelectedSimHandler->simulator->psim()->memtri.rows());

				ui::TreePop();
			}

			if (ui::TreeNodeEx("Solver info"))
			{
				shared_ptr<SolverInfo> solver_info = mSelectedSimHandler->simulator->finfo()[mCurrentFrame]->solver_info;
				ui::Text("Total iterations: %d", solver_info ? solver_info->iter.size() : 0);
				ui::Text("Final status: %s", solver_info ? solver_info->status.c_str() : "-");
				ui::Text("Energy decrement: %f", solver_info ? solver_info->iter.front().energy - solver_info->iter.back().energy : 0);

				ui::BeginChild("scroll_solver_iters", ImVec2(0, 200), true, ImGuiWindowFlags_HorizontalScrollbar);
				ui::Columns(3, "solver_info_columns");
				ui::Text("Iter#"); ui::NextColumn();
				ui::Text("Energy"); ui::NextColumn();
				ui::Text("Residual"); ui::NextColumn();
				ui::Separator();
				for (int i = 0; solver_info && i < solver_info->iter.size(); ++i)
				{
					ui::Text("%d", i); ui::NextColumn();
					ui::Text("%f", solver_info->iter[i].energy); ui::NextColumn();
					ui::Text("%f", solver_info->iter[i].residual); ui::NextColumn();
				}
				ui::EndChild();

				ui::TreePop();
			}

			if (ui::TreeNodeEx("Global info", ImGuiTreeNodeFlags_DefaultOpen))
			{
				ui::PushItemWidth(80);
				ui::InputFloat("Total segment energy", &mSelectedSimHandler->simulator->finfo()[mCurrentFrame]->total_segment_energy, -1, ImGuiInputTextFlags_ReadOnly);
				ui::InputFloat("Total cross energy", &mSelectedSimHandler->simulator->finfo()[mCurrentFrame]->total_cross_energy, -1, ImGuiInputTextFlags_ReadOnly);
				ui::InputFloat("Total membrane energy", &mSelectedSimHandler->simulator->finfo()[mCurrentFrame]->total_membrane_energy, -1, ImGuiInputTextFlags_ReadOnly);
				ui::TreePop();
			}

			if (ui::TreeNodeEx("Actuator info", ImGuiTreeNodeFlags_DefaultOpen))
			{
				ui::PushItemWidth(80);
				if (mPickedActuator < 0)
				{
					static string u_selected_actuator = "none";
					ui::InputText("Actuator ID", &u_selected_actuator, ImGuiInputTextFlags_ReadOnly);
				}
				else
				{
					ui::InputInt("Actuator ID", &mPickedActuator);
					mPickedActuator = min(mPickedActuator, int(mSelectedSimHandler->simulator->get_smp()->FF.rows() - 1));
					mPickedActuator = max(mPickedActuator, 0);
				}

				static float u_length = 0;
				u_length = (mPickedActuator < 0) ? 0 : mSelectedSimHandler->simulator->get_actuator_length_at_frame(mPickedActuator, mCurrentFrame);
				ui::InputFloat("Current length", &u_length);

				static float u_theta = 0;
				u_theta = (mPickedActuator < 0) ? 0 : mSelectedSimHandler->simulator->get_smp()->theta(mPickedActuator);
				ui::InputFloat("Target angle", &u_theta);

				static string u_boundary = "-";
				u_boundary = (mPickedActuator < 0) ? "-" : (mSelectedSimHandler->simulator->get_smp()->FFb(mPickedActuator) == 1) ? "YES" : "NO";
				ui::InputText("At boundary", &u_boundary);

				static float u_scale3d = 0;
				u_scale3d = (mPickedActuator < 0) ? 0 : mSelectedSimHandler->simulator->get_smp()->scale3d(mPickedActuator);
				ui::InputFloat("Scale 3d", &u_scale3d);

				ui::PushItemWidth(80);

				static int u_selected_actuator_faces[2];
				u_selected_actuator_faces[0] = (mPickedActuator == -1) ? -1 : mSelectedSimHandler->simulator->get_smp()->FF(mPickedActuator, 0);
				u_selected_actuator_faces[1] = (mPickedActuator == -1) ? -1 : mSelectedSimHandler->simulator->get_smp()->FF(mPickedActuator, 1);
				ui::InputInt2("Base IDs", u_selected_actuator_faces, ImGuiInputTextFlags_ReadOnly);

				static float spth[2];
				if (mPickedActuator < 0) spth[0] = spth[1] = 0;
				else
				{
					spth[0] = mSelectedSimHandler->simulator->get_smp()->spth(mPickedActuator, 0);
					spth[1] = mSelectedSimHandler->simulator->get_smp()->spth(mPickedActuator, 1);
				}
				ui::InputFloat2("Thickness", spth, -1, ImGuiInputTextFlags_ReadOnly);

				ui::Columns(3, "actuator_columns");
				ui::Separator();
				ui::Text("length"); ui::NextColumn();
				ui::Text("restlen"); ui::NextColumn();
				ui::Text("force"); ui::NextColumn();
				ui::Separator();

				if (mPickedActuator < 0)
				{
					for (int i = 0; i < 4 * 3; ++i)
					{
						ui::Text("-"); ui::NextColumn();
					}
				}
				else
				{
					static Eigen::Vector4d len;
					mSelectedSimHandler->simulator->get_bracket_length_at_frame(mPickedActuator, mCurrentFrame, len);

					static Eigen::Vector4d restlen;
					mSelectedSimHandler->simulator->get_bracket_restlen_at_frame(mPickedActuator, mCurrentFrame, restlen);

					static Eigen::Vector4d bracket_forces;
					mSelectedSimHandler->simulator->get_bracket_forces_at_frame(mPickedActuator, mCurrentFrame, bracket_forces);

					for (int i = 0; i < 4; ++i)
					{
						ui::Text("%f", len(i)); ui::NextColumn();
						ui::Text("%f", restlen(i)); ui::NextColumn();
						ui::Text("%f", bracket_forces(i)); ui::NextColumn();
					}
				}
				ui::Columns(1);
				ui::Separator();

				//static float u_mem_energy = 0;
				//if (mPickedActuator >= 0) u_mem_energy = ssym->get_actuator_membrane_energy(mPickedActuator, mCurrentFrame);
				//else u_mem_energy = 0;
				//ui::InputFloat("Membrane energy", &u_mem_energy);
				ui::TreePop();
			}

			if (shandler_erase) ui::OpenPopup("Remove simulation");
			if (ui::BeginPopupModal("Remove simulation", NULL, ImGuiWindowFlags_AlwaysAutoResize))
			{
				ui::Text("%s", ("Simulation '" + shandler_erase->simulator->get_name() + "' will be removed. Are you sure?").c_str());
				if (ui::Button("Remove"))
				{
					removeSimulationHandler(shandler_erase);
					ui::CloseCurrentPopup();
				}
				ui::SameLine();
				if (ui::Button("Cancel"))
				{
					shandler_erase.reset();
					ui::CloseCurrentPopup();
				}
				ui::EndPopup();
			}
		}

	}
	ui::End();
}

void SmpApp::GuiSurfaceImport()
{
	static bool popen = true;
	if (!popen)
	{
		mShowSurfaceImport = false;
		popen = true;
	}
	if (!ui::Begin("Surface import", &popen)) { ui::End(); return; }

#ifdef USE_GEOGRAM
	shared_ptr<GeoRemesher> geo_remesher = mSmpRemeshHandler.geo_remesher;

	if (ui::CollapsingHeader("Visualization", ImGuiTreeNodeFlags_DefaultOpen))
	{
		ui::Checkbox("Show imported surface", &mShowSurfaceImported);
		ui::Checkbox("Show flattened surface", &mShowSurfaceFlattened);
		ui::Checkbox("Show remeshed surface", &mShowSurfaceRemeshed);
	}
	if (ui::CollapsingHeader("Geometry", ImGuiTreeNodeFlags_DefaultOpen))
	{
		if (ui::Button("Import surface from file..."))
		{
			import_surface();
		}
		if (geo_remesher->has_imported_mesh())
		{
			ui::PushItemWidth(80);

			ui::Text("Imported mesh");
			int num_vertices = geo_remesher->get_V_imported().rows();
			ui::InputInt("Vertices", &num_vertices, 0, 0, ImGuiInputTextFlags_ReadOnly);
			int num_faces = geo_remesher->get_F_imported().rows();
			ui::InputInt("Faces", &num_faces, 0, 0, ImGuiInputTextFlags_ReadOnly);
			
			ui::Separator();
			ui::Text("Flattening");
			if (ui::Button("Recompute flattening"))
			{
				geo_remesher->flatten_imported_mesh();
				mSmpRemeshHandler.update_ImportSurfacesBatch();
			}
			ui::InputFloat("Max flat stretch", geo_remesher->get_max_stretch());

			ui::Separator();

			ui::Text("Remeshing");
			ui::Text("N.B.: consider remesh scaling factors");
			static float u_edgelen = 12;
			static int u_vertices = 100;
			if (ui::Button("By target edgelen")) u_vertices = geo_remesher->vertices_by_target_edgelen(u_edgelen);
			ui::SameLine();
			ui::InputFloat("", &u_edgelen);
			ui::InputInt("Target number of vertices", &u_vertices);

			static int u_lloyd_iter = 50;
			static int u_newton_iter = 3000;
			static int u_newton_m = 30;
			if (ui::TreeNode("Remeshing parameters"))
			{
				ui::InputInt("Lloyd iterations", &u_lloyd_iter, 0, 0);
				u_lloyd_iter = std::max(u_lloyd_iter, 1);
				ui::InputInt("Newton iterations", &u_newton_iter, 0, 0);
				u_newton_iter = std::max(u_newton_iter, 1);
				ui::InputInt("Newton M", &u_newton_m, 0, 0);
				u_newton_m = std::max(u_newton_m, 1);
				ui::TreePop();
			}

			if (ui::Button("Remesh"))
			{
				geo_remesher->remesh(u_vertices, u_lloyd_iter, u_newton_iter, u_newton_m);
				mSmpRemeshHandler.update_ImportSurfacesBatch();
			}

			if (geo_remesher->has_remeshed_mesh())
			{
				ui::Separator();
				ui::Text("Post-processing");
				if (ui::Button("Flip faces"))
				{
					geo_remesher->flip_remeshed_faces();
					mSmpRemeshHandler.update_ImportSurfacesBatch();
				}
				if (ui::Button("Remove isolated triangles"))
				{
					geo_remesher->remove_isolated_trangles();
					mSmpRemeshHandler.update_ImportSurfacesBatch();
				}

				ui::Separator();
				ui::Text("Remeshed surface");
				int num_vertices = geo_remesher->get_V_remeshed().rows();
				ui::InputInt("Vertices", &num_vertices, 0, 0, ImGuiInputTextFlags_ReadOnly);
				int num_faces = geo_remesher->get_F_remeshed().rows();
				ui::InputInt("Faces", &num_faces, 0, 0, ImGuiInputTextFlags_ReadOnly);
				float mean_edgelen = geo_remesher->get_remeshed_mean_edgelen();
				ui::InputFloat("Mean edge length", &mean_edgelen, 0, 0, 0, ImGuiInputTextFlags_ReadOnly);

				static bool auto_load = false;
				if (ui::Button("Save"))
				{
					fs::path path = app::getSaveFilePath(fs::path(), vector<string>{"obj"});
					if (path.empty()) return;
					string fname = path.string();
					if (fname.substr(fname.size() - 4, 4).compare(".obj") != 0) fname += ".obj";
					geo_remesher->write_remeshed(fname);
					if (auto_load)
					{
						read_file_info(fname);
						model_open_stencil_request = true;
					}
				}
				ui::SameLine();
				ui::Checkbox("load", &auto_load); tooltip("Automatically load the structure after saving.");
			}
		}
	}
#else
	ui::Text("This version is compiled without Geogram. Remesh not possible.");
#endif

	ui::End();
}

void SmpApp::GuiLandscapeEditor()
{
	if (!mSelectedSimHandler) return;

	static bool popen = true;
	if (!popen)
	{
		mShowTimeLandscapeEditor = false;
		popen = true;
	}
	if (!ui::Begin("Time landscape editor", &popen)) { ui::End(); return; }

	if (ui::CollapsingHeader("Visualization", ImGuiTreeNodeFlags_DefaultOpen))
	{
		ui::Checkbox("Show time landscape", &mShowTimeLandscape);
		ui::Checkbox("Show time ranges", &mShowTimeRanges);

		if (ui::Button("Compute time ranges"))
		{
			const Smp& smp = *mSelectedSimHandler->simulator->get_smp();
			Eigen::MatrixXd& time_ranges = mSelectedSimHandler->time_ranges;
			time_ranges.resize(smp.FF.rows(), 4);
			time_ranges.col(0).setConstant(0);
			time_ranges.col(1).setConstant(120);
			time_ranges.col(2).setConstant(0);
			time_ranges.col(3).setConstant(120);
			Eigen::MatrixXd Vmin; Vmin.setConstant(smp.FF.rows(), 3, 0);
			Eigen::MatrixXd Vmax; Vmax.setConstant(smp.FF.rows(), 3, 120);
			Eigen::MatrixXd V1(smp.FF.rows(), 3);
			for (int i = 0; i < smp.FF.rows(); ++i)
			{
				Vmin.row(i).head(2) = 0.5 * (smp.Fc2d.row(smp.FF(i, 0)) + smp.Fc2d.row(smp.FF(i, 0)));
				Vmax.row(i).head(2) = Vmin.row(i).head(2);
				double actuator_length = mSelectedSimHandler->simulator->get_actuator_length_at_frame(i, 0);
				for (int j = 0; j < 2; ++j)
				{
					Eigen::Vector2i seg_id = mSelectedSimHandler->simulator->get_actuators()[i]->segment_springs.segment(j * 2, 2);
					double target_length = (mSelectedSimHandler->simulator->psim()->bumpers(seg_id(0)) + mSelectedSimHandler->simulator->psim()->bumpers(seg_id(1))) * 0.5;
					bool weak = false;
					if (mSelectedSimHandler->simulator->get_settings().mem_model == StaticSimulatorStruct::MEM_FINE && smp.FFb(i)) weak = true;
					if (mSelectedSimHandler->simulator->get_settings().mem_model == StaticSimulatorStruct::MEM_LINEAR &&
						mSelectedSimHandler->simulator->get_settings().mimic_bounary_weakening && smp.FFb(i)) weak = true;
					Eigen::Vector2d time_range = mBracketConfigurer->compute_time_range(actuator_length, target_length, weak);
					Vmin(i, 2) = max(Vmin(i, 2), time_range(0));
					Vmax(i, 2) = min(Vmax(i, 2), time_range(1));
					time_ranges(i, 0 + j * 2) = time_range(0);
					time_ranges(i, 1 + j * 2) = time_range(1);
				}

				//time_ranges(i, 0) = max(time_ranges(i, 0), time_ranges(i, 2));//
				//time_ranges(i, 1) = max(time_ranges(i, 1), time_ranges(i, 3));//
				//time_ranges(i, 2) = time_ranges(i, 3) = 0;
			}
		}
	}

	if (ui::CollapsingHeader("Editing", ImGuiTreeNodeFlags_DefaultOpen))
	{
		bool time_edited = false;
		if (ui::TreeNodeEx("All actuators", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ui::PushItemWidth(80);
			if (ui::InputInt("Enforce max", &mTimeLandscapeMax, 1, 100, ImGuiInputTextFlags_EnterReturnsTrue)) time_edited = true;
			if (ui::InputInt("Enforce min", &mTimeLandscapeMin, 1, 100, ImGuiInputTextFlags_EnterReturnsTrue)) time_edited = true;
			ui::PopItemWidth();

			ui::Spacing();

			static int u_time_to_set = 30;
			ui::PushItemWidth(80);
			ui::InputInt("##input_time_to_set", &u_time_to_set);
			ui::PopItemWidth();
			ui::SameLine();
			if (ui::Button("Set constant"))
			{
				mSelectedSimHandler->time_landscape.setConstant(u_time_to_set);
				time_edited = true;
			}

			static int u_time_to_add = 3;
			ui::PushItemWidth(80);
			ui::InputInt("##input_time_to_add", &u_time_to_add);
			ui::PopItemWidth();
			ui::SameLine();
			if (ui::Button("Add"))
			{
				mSelectedSimHandler->time_landscape.array() += u_time_to_add;
				time_edited = true;
			}
			ui::SameLine();
			if (ui::Button("Sub"))
			{
				mSelectedSimHandler->time_landscape.array() -= u_time_to_add;
				time_edited = true;
			}

			ui::TreePop();
		}

		if (ui::TreeNodeEx("Time brush", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ui::PushItemWidth(100);
			ui::Combo("Brush mode", (int*)&mTimeBrushMode, "Sculpt\0Smoothen\0");

			ui::PushItemWidth(40);
			ui::InputInt("##time_radius_val", &mTimeLandscapeRadius, 0, 0);
			ui::PopItemWidth();
			ui::SameLine();
			if (ui::Button("-##time_radius")) mTimeLandscapeRadius -= mTimeLandscapeRadius / 3;
			ui::SameLine();
			if (ui::Button("+##time_radius")) mTimeLandscapeRadius += mTimeLandscapeRadius / 2;
			ui::SameLine();
			ui::Text("Radius");
			mTimeLandscapeRadius = min(max(mTimeLandscapeRadius, 2), 1000);

			ui::PushItemWidth(40);
			ui::InputInt("##time_inc_val", &mTimeLandscapeIncrement, 0, 0);
			ui::PopItemWidth();
			ui::SameLine();
			if (ui::Button("-##time_inc")) mTimeLandscapeIncrement -= mTimeLandscapeIncrement / 3;
			ui::SameLine();
			if (ui::Button("+##time_inc")) mTimeLandscapeIncrement += mTimeLandscapeIncrement / 2;
			ui::SameLine();
			if (ui::Button("inv")) mTimeLandscapeIncrement *= -1;
			ui::SameLine();
			ui::Text("Increment");

			ui::TreePop();
		}

		ui::Spacing();

		if (ui::Button("Configure thickness"))
		{
			mSelectedSimHandler->simulator.reset(mSelectedSimHandler->simulator->configure_brackets(mSelectedSimHandler->time_landscape).release());
		}
		ui::SameLine();
		if (ui::Button("Restore"))
		{
			mSelectedSimHandler->time_landscape = mSelectedSimHandler->simulator->get_smp()->time_landscape;
			mSelectedSimHandler->update_TimeLandscapeBatch();
		}

		if (time_edited)
		{
			for (int i = 0; i < mSelectedSimHandler->time_landscape.size(); ++i)
			{
				mSelectedSimHandler->time_landscape(i) = min(mSelectedSimHandler->time_landscape(i), double(mTimeLandscapeMax));
				mSelectedSimHandler->time_landscape(i) = max(mSelectedSimHandler->time_landscape(i), double(mTimeLandscapeMin));
			}
			mSelectedSimHandler->update_TimeLandscapeBatch();
		}
	}
	ui::End();
}

void SmpApp::GuiSolverInfo()
{
	if (!mSelectedSimHandler) return;

	static bool popen = true;
	if (!popen)
	{
		mShowCurrentSolverInfo = false;
		popen = true;
	}
	if (!ui::Begin("Solver info", &popen)) { ui::End(); return; }

	auto simulator = mSelectedSimHandler->simulator;

	static int selected_timestep = -1;
	if (simulator->get_simulated_frames() < 2) selected_timestep = -1;
	if (0 <= selected_timestep)
	{
		if (simulator->get_simulated_frames() < selected_timestep - 1) selected_timestep = -1;
		else if (!simulator->finfo()[selected_timestep]->solver_info) selected_timestep = -1;
	}

	if (ui::CollapsingHeader("Solver time steps", ImGuiTreeNodeFlags_DefaultOpen))
	{
		static bool solver_info_scroll = false;
		ui::Checkbox("Auto scroll##timesteps", &solver_info_scroll);
		ui::BeginChild("scroll_solver_timesteps", ImVec2(0, 200), true, ImGuiWindowFlags_HorizontalScrollbar);

		ui::Columns(10, "solver_info_columns");
		ui::Text("Now"); ui::NextColumn();
		ui::Text("Step#"); ui::NextColumn();
		ui::Text("Status"); ui::NextColumn();
		ui::Text("Iters"); ui::NextColumn();
		ui::Text("Condition"); ui::NextColumn();
		ui::Text("Ini. Energy"); ui::NextColumn();
		ui::Text("Fin. Energy"); ui::NextColumn();
		ui::Text("Ini. Residual"); ui::NextColumn();
		ui::Text("Fin. Residual"); ui::NextColumn();
		ui::Text("Time"); ui::NextColumn();
		ui::Separator();

		static bool initialize_column_widths = true;
		if (initialize_column_widths)
		{
			vector<int> col_init_widths{ 35, 50, 70, 50, 100, -1, -1, -1, -1, -1 };
			for (int i = 0; i < col_init_widths.size(); ++i)
				if (0 <= col_init_widths[i]) ui::SetColumnWidth(i, col_init_widths[i]);
			initialize_column_widths = false;
		};

		for (int i = 1; i < simulator->get_simulated_frames(); ++i)
		{
			auto solver_info = simulator->finfo()[i]->solver_info;
			if (!solver_info) continue;

			string now_check = (i == mCurrentFrame ? "*" : "") + string("##iter_") + to_string(i);
			if (ui::Selectable(now_check.c_str(), selected_timestep == i, ImGuiSelectableFlags_SpanAllColumns))
				selected_timestep = (selected_timestep == i ? -1 : i);
			ui::NextColumn();
			ui::Text("%d", i); ui::NextColumn();
			ui::Text("%s", solver_info->status.c_str()); ui::NextColumn();
			ui::Text("%d", solver_info->num_iter); ui::NextColumn();
			ui::Text("%f", solver_info->condition); ui::NextColumn();
			ui::Text("%f", solver_info->iter.front().energy); ui::NextColumn();
			ui::Text("%f", solver_info->iter.back().residual); ui::NextColumn();
			ui::Text("%f", solver_info->iter.front().residual); ui::NextColumn();
			ui::Text("%f", solver_info->iter.back().residual); ui::NextColumn();
			ui::Text("%.3f", solver_info->time); ui::NextColumn();
		}

		if (solver_info_scroll) ui::SetScrollHere(0.999f);
		ui::EndChild();
	}

	if (ui::CollapsingHeader("Solver iterations", ImGuiTreeNodeFlags_DefaultOpen))
	{
		auto solver_info = simulator->psim()->get_current_solver_info();
		if (selected_timestep < 0) ui::Text("Time step: current");
		else
		{
			ui::Text("Time step: %d", selected_timestep);
			solver_info = simulator->finfo()[selected_timestep]->solver_info;
			ui::SameLine();
			if (ui::SmallButton("switch to current")) selected_timestep = -1;
		}

		static bool solver_iters_scroll = false;
		ui::Checkbox("Auto scroll##iters", &solver_iters_scroll);
		ui::BeginChild("scroll_solver_iters", ImVec2(0, 200), true, ImGuiWindowFlags_HorizontalScrollbar);

		ui::Columns(3, "solver_iters_columns");
		ui::Text("Iter#"); ui::NextColumn();
		ui::Text("Energy"); ui::NextColumn();
		ui::Text("Residual"); ui::NextColumn();
		ui::Separator();
		if (solver_info)
		{
			for (int i = 0; solver_info && i < solver_info->iter.size(); ++i)
			{
				ui::Text("%d", solver_info->iter[i].num); ui::NextColumn();
				ui::Text("%f", solver_info->iter[i].energy); ui::NextColumn();
				ui::Text("%f", solver_info->iter[i].residual); ui::NextColumn();
			}
		}

		if (solver_iters_scroll) ui::SetScrollHere(0.999f);
		ui::EndChild();
	}

	//
	ui::End();
}

void SmpApp::GuiConsole()
{
	static bool popen = true;
	if (!popen)
	{
		mShowConsole = false;
		popen = true;
	}
	if (!ui::Begin("Console", &popen)) { ui::End(); return; }
	ui::PushItemWidth(-1);
	static string u_buf_console;
	u_buf_console = buf_console.str();
	ui::InputTextMultiline("", &u_buf_console, ImVec2(0, ui::GetWindowHeight() - 40), ImGuiInputTextFlags_ReadOnly);
	//ui::TextUnformatted(buf_console.str().c_str());
	static int u_textlen = 0;
	if (u_textlen != buf_console.str().size())
	{
		if (u_textlen < buf_console.str().size())
		{
			ImGui::SetScrollHere(1.0f);
			u_textlen += (buf_console.str().size() - u_textlen + 1) / 2;
		}
		else
			u_textlen = buf_console.str().size();
	}
	ui::End();
}

void SmpApp::GuiBracketExplorer()
{
	static bool popen = true;
	if (!popen)
	{
		mShowBracketExplorer = false;
		popen = true;
	}
	if (!ui::Begin("Bracket explorer", &popen)) { ui::End(); return; }
	static Eigen::VectorXf u_stressstrain(100);
	static float u_time = 0;
	static float u_thickness = 0.4;
	static float u_length = 4.805;
	u_stressstrain.setZero();
	Eigen::MatrixXd poly;
	Eigen::VectorXd th(1); th(0) = u_thickness;
	Eigen::VectorXd le(1); le(0) = u_length;
	mBracket->make_energy_poly_coeff(th, le, u_time, poly);
	for (int i = 0; i < 100; ++i)
	{
		double pos = u_length * i / 100.;
		double pow_pos = 1;
		for (int j = 1; j < poly.cols(); ++j)
		{
			u_stressstrain(i) += pow_pos * poly(0, j) * j;
			pow_pos *= pos;
		}
	}
	u_stressstrain *= 4;
	ui::PushItemWidth(-1);
	ui::Text("Stress-strain for a single actuator");

#ifdef USE_MATLAB
	if (matlab_started)
	{
		ui::SameLine();
		if (ui::SmallButton("Matlab plot"))
		{
			matlab.figure("bracket_stress_strain", "Bracket stress-strain", "stain", "stress", false);
			Eigen::VectorXd strains(u_stressstrain.size());
			for (int i = 0; i < strains.size(); ++i)
				strains(i) = double(i) / strains.size();
			matlab.plot(strains, u_stressstrain.cast<double>());
			matlab.eval("xlim([0 1]);");
			matlab.eval("ylim([0 10]);");
			matlab.eval("grid on");
		}
	}
#endif

	ui::PlotLines("", u_stressstrain.data(), u_stressstrain.size(), 0, "", 0, 10, ImVec2(0, 300));

	ui::PushItemWidth(ui::GetWindowWidth() - 100);
	ui::SliderFloat("time", &u_time, mBracket->get_grid_min()(2) * 60, (mBracket->get_grid_min()(2) + mBracket->get_grid_step()(2) * (mBracket->get_grid_num()(2) - 1)) * 60);
	ui::SliderFloat("thickness", &u_thickness, mBracket->get_grid_min()(0), mBracket->get_grid_min()(0) + mBracket->get_grid_step()(0) * (mBracket->get_grid_num()(0) - 1));
	ui::SliderFloat("legnth", &u_length, mBracket->get_grid_min()(1), mBracket->get_grid_min()(1) + mBracket->get_grid_step()(1) * (mBracket->get_grid_num()(1) - 1));
	//ui::SliderFloat("time", &u_time, -60, 150);
	//ui::SliderFloat("thickness", &u_thickness, 0, 1);
	//ui::SliderFloat("legnth", &u_length, 0, 15);

	ui::Separator();

	static float u_strain = 0;
	ui::SliderFloat("Strain", &u_strain, 0, 1);
	static float u_displacement_force[2];
	u_displacement_force[0] = u_length * u_strain;
	u_displacement_force[1] = 0;
	double pow_pos = 1;
	for (int j = 1; j < poly.cols(); ++j)
	{
		u_displacement_force[1] += pow_pos * poly(0, j) * j;
		pow_pos *= u_displacement_force[0];
	}
	u_displacement_force[1] *= 4;
	ui::InputFloat2("displ/force", u_displacement_force, -1, ImGuiInputTextFlags_ReadOnly);

	ui::End();
}

void SmpApp::GuiPlottingWindow()
{
	static bool popen = true;
	if (!popen)
	{
		mShowPlottingWindow = false;
		popen = true;
	}
	if (!ui::Begin("Plotting", &popen, ImGuiWindowFlags_NoCollapse)) { ui::End(); return; }

#ifdef USE_MATLAB
	if (matlab_started)
	{

		static Eigen::MatrixXd cmap;
		if (cmap.size() == 0)
		{
			ifstream ss(mCoreSettings->get_smpup_path() / "cmap_GnBu.json");
			if (ss.good())
			{
				cereal::JSONInputArchive archive(ss);
				archive(cereal::make_nvp("cmap", cmap));
			}
			else clog << "Error: Could not load color map\n";
		}

		if (ui::CollapsingHeader("General plots", ImGuiTreeNodeFlags_DefaultOpen))
		{

			static float u_init_len = 8.0f;
			static float u_max_strain = 0.7f;
			ui::InputFloat("Actuator length", &u_init_len);
			ui::InputFloat("Clip max strain", &u_max_strain);
			u_max_strain = min(max(u_max_strain, 0.0f), 1.0f);
			if (ui::Button("Plot time ranges"))
			{
				double bracket_length_correction = 2 * mSmpOptions.base_r * cos(37. / 180. * M_PI); //TODO: compute correctly via smp options
				double init_mem_length = u_init_len - bracket_length_correction;
				Eigen::VectorXd mem_poly = mLinearMemModel->make_force_poly_coeff(init_mem_length);

				double th_min = mBracket->get_grid_min()(BracketModel::GRID_THICK);
				double th_stp = mBracket->get_grid_step()(BracketModel::GRID_THICK);
				int th_num = mBracket->get_grid_num()(BracketModel::GRID_THICK);

				int ncol = cmap.rows();
				Eigen::MatrixXd sub_cmap(th_num, 3);

				matlab.eval("set(0,'defaulttextinterpreter','latex')");

				stringstream ss;
				ss << "$l = " << setprecision(2) << u_init_len << "$ mm";
				matlab.figure("time_ranges", ss.str(), "time, s", "displacement, mm", true);
				for (int j = 0; j < th_num; ++j)
				{
					double thick = th_min + j * th_stp;
					Eigen::VectorXd times = natural_seriesd(0, 40) * 3.;
					Eigen::VectorXd disp(times.size());
					for (int i = 0; i < times.size(); ++i)
					{
						disp(i) = mBracket->compute_mem_deformation(thick, init_mem_length, times(i), mem_poly, init_mem_length);
						if (u_max_strain < disp(i) / init_mem_length)
						{
							disp.conservativeResize(i);
							times.conservativeResizeLike(disp);
							break;
						}
					}
					sub_cmap.row(j) = cmap.row(floor(double(j) / th_num  * ncol));
					matlab.plot(times, disp, "", sub_cmap.row(j), 2);
				}

				matlab.eval("ax = gca; ax.FontSize = 8; set(gca,'TickLabelInterpreter', 'latex');");
				matlab.put_variable("cmap", sub_cmap);
				matlab.eval("caxis([0.3 0.69999]);");
				matlab.eval("xlim([0 120]);");
				matlab.eval("set(gca, 'xtick', 0:30:120);");
				matlab.eval("grid on");
				matlab.eval("colormap(cmap);");
				matlab.eval("hcb = colorbar;");
				matlab.eval("colorTitleHandle = get(hcb,'Title');");
				matlab.eval("colorTitleHandle.Interpreter = 'latex';");
				matlab.eval("set(colorTitleHandle ,'String','$h$, mm');");
				matlab.eval("hcb.TickLabelInterpreter = 'latex';");
				stringstream ss_size;
				double width = 6;
				double height = 5;
				ss_size << "set(gcf, 'Units', 'centimeters', 'Position', [0, 0, " << width << ", " << height << "], 'PaperUnits', 'centimeters', 'PaperSize', [" << width << ", " << height << "]);";
				matlab.eval(ss_size.str());

				// Infer time ranges for target displacements
				//Eigen::VectorXd target_deformations = init_mem_length * natural_seriesd(0, 20) / 20 * u_max_strain;
				//Eigen::VectorXd target_lengths = (double)init_mem_length - target_deformations.array();
				//Eigen::MatrixXd time_ranges(target_lengths.size(), th_num);
				//for (int i = 0; i < time_ranges.rows(); ++i)
				//{
				//	double target_length = target_lengths(i);
				//	double target_mem_length = target_length;
				//	Eigen::VectorXd target_mem_length_pow = Eigen::ArrayXd::Constant(2, -target_mem_length).pow(natural_seriesd(2));
				//	double target_mem_force = -mem_poly.dot(target_mem_length_pow); // negative since the force is along negative direction
				//	for (int j = 0; j < time_ranges.cols(); ++j)
				//		time_ranges(i, j) = mBracket->compute_time(init_mem_length, target_length, target_mem_force / 4.0, th_min + th_stp * j);
				//}
				//for (int j = 0; j < time_ranges.cols(); ++j)
				//{
				//	int i = time_ranges.rows() - 1;
				//	while (0 < i && time_ranges(i, j) == time_ranges(i - 1, j)) --i;
				//	//++i;
				//	matlab.plot(time_ranges.col(j).head(i), target_deformations.head(i), "", Eigen::Vector3d{ 0.8, 0.6, 0.4 });
				//}
			}

			ui::Spacing();

			if (ui::Button("Plot linearized membrane"))
			{
				double bracket_length_correction = 2 * mSmpOptions.base_r * cos(37. / 180. * M_PI);

				double len_min = 5;
				double len_stp = 1;
				int len_num = 5;

				int ncol = cmap.rows();
				Eigen::MatrixXd sub_cmap(len_num, 3);

				matlab.eval("set(0,'defaulttextinterpreter','latex')");
				matlab.figure("lin_membrane", "Linearized membrane", "displacement, mm", "force, N", true);

				for (int j = 0; j < len_num; ++j)
				{
					double init_len = len_min + j * len_stp;
					sub_cmap.row(j) = cmap.row(floor(double(j) / len_num  * ncol));

					double init_mem_length = init_len - bracket_length_correction;
					Eigen::VectorXd mem_poly = mLinearMemModel->make_force_poly_coeff(init_mem_length);
					Eigen::VectorXd disp = natural_seriesd(20) / 20. * init_mem_length * 0.5;
					Eigen::MatrixXd mem_length_pow = (disp.array() - init_mem_length).matrix() * Eigen::VectorXd::Ones(2).transpose();
					for (int i = 0; i < mem_length_pow.cols(); ++i)
						mem_length_pow.col(i) = mem_length_pow.array().col(i).pow(i);

					Eigen::VectorXd mem_forces = -mem_length_pow * mem_poly;
					matlab.plot(disp, mem_forces, "--", sub_cmap.row(j), 1);

					SmpSettings settings;
					settings.mem_model = StaticSimulatorStruct::MEM_FINE;
					settings.timestep_mul = 1;
					settings.boundary_weakening = 0.3;

					double scissor_tmp = mSmpOptions.scissor;
					mSmpOptions.scissor = 0.1;
					mSmpOptions.save();
					new_regular(7, 4, init_len, 0.4, settings);
					mSmpOptions.scissor = scissor_tmp;
					mSmpOptions.save();

					mSelectedSimHandler->time_landscape.setConstant(30);
					mSelectedSimHandler->simulator.reset(mSelectedSimHandler->simulator->configure_brackets(mSelectedSimHandler->time_landscape).release());
					mSelectedSimHandler->simulator->simulate();

					while (mSelectedSimHandler->simulator->is_running_task())
					{
						//if (30 < mSelectedSimHandler->simulator->get_time_simulated()) mSelectedSimHandler->simulator->stop_simulation(); // according to time landscape
						sleep(100);
					}

					int central_actuator = 65; // for regular 7x4
					Eigen::VectorXd lengths(mSelectedSimHandler->simulator->get_simulated_frames() - 1);
					Eigen::VectorXd forces(lengths.size());
					for (int i = 0; i < forces.rows(); ++i)
					{
						Eigen::Vector4d bracket_forces;
						mSelectedSimHandler->simulator->get_bracket_forces_at_frame(central_actuator, i + 1, bracket_forces);
						forces(i) = bracket_forces.sum();
						lengths(i) = mSelectedSimHandler->simulator->get_actuator_length_at_frame(central_actuator, i + 1);
					}
					int c = lengths.size() - 1;
					while (init_mem_length * 0.5 < init_len - lengths(c)) --c;
					lengths.conservativeResize(c);
					forces.conservativeResize(c);

					matlab.plot(init_len - lengths.array(), forces, "", sub_cmap.row(j), 2);
				}

				matlab.eval("ax = gca; ax.FontSize = 8; set(gca,'TickLabelInterpreter', 'latex');");
				matlab.put_variable("cmap", sub_cmap);
				matlab.eval("caxis([5 9.9999]);");
				//matlab.eval("xlim([0 120]);");
				//matlab.eval("set(gca, 'xtick', 0:30:120);");
				matlab.eval("grid on");
				matlab.eval("colormap(cmap);");
				matlab.eval("hcb = colorbar;");
				matlab.eval("colorTitleHandle = get(hcb,'Title');");
				matlab.eval("colorTitleHandle.Interpreter = 'latex';");
				matlab.eval("set(colorTitleHandle ,'String','$l$, mm');");
				matlab.eval("hcb.TickLabelInterpreter = 'latex';");
				stringstream ss_size;
				double width = 6;
				double height = 5;
				ss_size << "set(gcf, 'Units', 'centimeters', 'Position', [0, 0, " << width << ", " << height << "], 'PaperUnits', 'centimeters', 'PaperSize', [" << width << ", " << height << "]);";
				matlab.eval(ss_size.str());

			}
		}

		if (ui::CollapsingHeader("Current simulation plots", ImGuiTreeNodeFlags_DefaultOpen))
		{

			if (mSelectedSimHandler && ui::Button("Plot time landscape"))
			{
				matlab.eval("set(0,'defaulttextinterpreter','latex')");
				matlab.figure("time_landscape", "Time landscape", "", "", true);
				matlab.eval("ax = gca; ax.FontSize = 8; set(gca,'TickLabelInterpreter', 'latex');");
				matlab.eval("set(gca, 'xtick', [0 100]);");
				matlab.eval("set(gca, 'ytick', [0 100]);");
				matlab.eval("set(gca, 'ztick', [30 60 80]);");
				matlab.eval("axis equal");
				matlab.eval("grid on");

				Eigen::MatrixXd sub_cmap = cmap.topRows(cmap.rows() / 6 * 5).bottomRows(cmap.rows() / 3 * 2);
				matlab.put_variable("cmap", sub_cmap);
				matlab.eval("colormap(cmap);");
				matlab.eval("colorbar");
				matlab.eval("caxis([30 80]);");

				matlab.put_variable("T", (mSelectedSimHandler->simulator->get_smp()->F.array() + 1).cast<double>());
				matlab.put_variable("X", mSelectedSimHandler->simulator->get_smp()->uv.col(0));
				matlab.put_variable("Y", mSelectedSimHandler->simulator->get_smp()->uv.col(1));
				matlab.put_variable("Z", mSelectedSimHandler->simulator->get_smp()->time_landscape);
				matlab.eval("U = [X Y] * [0 1; -1 0]; X = U(:,1); Y = U(:,2);"); // for Cover
				matlab.eval("X = X - min(X); Y = Y - min(Y);");
				matlab.eval("xlim([0 max([X;Y])]);");
				matlab.eval("ylim([0 max([X;Y])]);");
				matlab.eval("zlim([0 max(Z)]);");
				matlab.eval("trisurf(T,X,Y,Z);");
				matlab.eval("shading interp");
				matlab.eval("daspect([1 1 0.5]);");

				matlab.eval("trisurf(T,X,Y,Z * 0, 'Facecolor',[0.5 0.5 0.5] , 'FaceAlpha',.1 , 'EdgeColor','none');");
				matlab.eval("trisurf(T,X * 0,Y,Z, 'Facecolor',[0.5 0.5 0.5] , 'FaceAlpha',.1 , 'EdgeColor','none');");
				matlab.eval("trisurf(T,X,Y * 0,Z, 'Facecolor',[0.5 0.5 0.5] , 'FaceAlpha',.1 , 'EdgeColor','none');");

				matlab.eval("hcb = colorbar;");
				matlab.eval("colorTitleHandle = get(hcb,'Title');");
				matlab.eval("colorTitleHandle.Interpreter = 'latex';");
				matlab.eval("set(colorTitleHandle ,'String','$t$, s');");
				matlab.eval("hcb.TickLabelInterpreter = 'latex';");
				matlab.eval("set(hcb,'YTick',[30:30:90])");
				stringstream ss_size;
				double width = 6;
				double height = 5;
				ss_size << "set(gcf, 'Units', 'centimeters', 'Position', [0, 0, " << width << ", " << height << "], 'PaperUnits', 'centimeters', 'PaperSize', [" << width << ", " << height << "]);";
				matlab.eval(ss_size.str());
			}

			static int u_side = 0;
			ui::Combo("", &u_side, "front\0back\0\0");
			ui::SameLine();
			if (mSelectedSimHandler && ui::Button("Plot thicknesses"))
			{
				matlab.figure("bracket_thickness", "Bracket thicknesses", "", "", true);

				Eigen::MatrixXd sub_cmap = cmap.topRows(cmap.rows() / 6 * 5);
				matlab.put_variable("cmap", sub_cmap);
				matlab.eval("colormap(cmap);");
				matlab.eval("colorbar");
				matlab.eval("caxis([0.3 0.65]);");

				auto smp = mSelectedSimHandler->simulator->get_smp();
				Eigen::MatrixXd p(smp->FF.rows(), 4);
				for (int i = 0; i < p.rows(); ++i)
				{
					int fi = smp->FF(i, 0);
					int ind = smp->fi(i, 0);
					p.row(i).head(2) = ( smp->uv.row(smp->F(fi, next3(ind))) + smp->uv.row(smp->F(fi, prev3(ind))) ) * 0.5;
					p.row(i).tail(2) = smp->spth.row(i);
				}

				matlab.put_variable("X", p.col(0));
				matlab.put_variable("Y", p.col(1));
				matlab.put_variable("Z1", p.col(2));
				matlab.put_variable("Z2", p.col(3));
				matlab.put_variable("S", Eigen::VectorXd::Constant(p.rows(), 3.0));

				matlab.eval("tri = alphaTriangulation(alphaShape(X,Y,10));");
				if (u_side == 0) matlab.eval("trisurf(tri,X,Y,Z1);");
				else matlab.eval("trisurf(tri,X,Y,Z2);");
				matlab.eval("shading interp");

				//matlab.eval("scatter3(X,Y,Z1,S,C1);");
				//matlab.eval("scatter3(X,Y,Z2,S);");
				matlab.eval("daspect([1 1 0.01]);");
				//matlab.eval("daspect([max(daspect)*[1 1] 1]);");
			}

			ui::Text("For selected actuator:");

			if (mSelectedSimHandler && -1 < mPickedActuator && ui::Button("Plot thickness config"))
			{
				matlab.eval("set(0,'defaulttextinterpreter','latex')");
				matlab.figure("thick_config", "Thickness configuration", "displacement, mm", "force, N", true);
				matlab.eval("grid on");

				const auto& smp = mSelectedSimHandler->simulator->get_smp();
				const auto& time_landscape = mSelectedSimHandler->time_landscape;
				const auto& simulator = mSelectedSimHandler->simulator;

				Eigen::VectorXd actuator_times(smp->FF.rows());
				for (int i = 0; i < smp->FF.rows(); ++i)
					actuator_times(i) = (time_landscape(smp->F(smp->FF(i, 0), next3(smp->fi(i, 0)))) + time_landscape(smp->F(smp->FF(i, 0), prev3(smp->fi(i, 0))))) * 0.5;

				Eigen::Vector2i act_edg = simulator->psim()->edges.row(simulator->psim()->linear(simulator->get_actuators()[mPickedActuator]->segment_springs(0))).transpose(); // since all actuator brackets are equal length in flat config
				double init_mem_length = (simulator->finfo()[0]->nodes.row(act_edg(0)) - simulator->finfo()[0]->nodes.row(act_edg(1))).norm();
				Eigen::VectorXd mem_poly = simulator->psim()->get_linear_mem_model()->make_force_poly_coeff(init_mem_length);

				Eigen::Vector2d target_length;
				for (int j = 0; j < 2; ++j)
				{
					Eigen::Vector2i seg_id = simulator->get_actuators()[mPickedActuator]->segment_springs.segment(j * 2, 2);
					target_length(j) = (simulator->psim()->bumpers(seg_id(0)) + simulator->psim()->bumpers(seg_id(1))) * 0.5;
				}

				double target_mem_length = target_length.mean();
				Eigen::VectorXd target_mem_length_pow = Eigen::ArrayXd::Constant(2, -target_mem_length).pow(natural_seriesd(2));
				double target_mem_force = -mem_poly.dot(target_mem_length_pow); // negative since the force is along negative direction

				Eigen::Vector2d X{ 0, 1 };
				Eigen::Vector2d Y{ 1, 1 };
				matlab.plot(X * init_mem_length, Y * target_mem_force / 2, "k--");

				Eigen::Vector2d th{
					simulator->psim()->get_bracket_model()->configure_thickness(init_mem_length, target_length(0), actuator_times(mPickedActuator), target_mem_force / 4),
					simulator->psim()->get_bracket_model()->configure_thickness(init_mem_length, target_length(1), actuator_times(mPickedActuator), target_mem_force / 4)
				};

				for (int j = 0; j < 2; ++j)
				{
					matlab.plot(Y * (init_mem_length - target_length(j)), X * target_mem_force / 2, "k--");
					
					Eigen::VectorXd br_poly = mBracket->make_force_poly_coeff(th(j), init_mem_length, actuator_times(mPickedActuator));
					Eigen::VectorXd disp = natural_seriesd(20) / 20. * init_mem_length;

					Eigen::MatrixXd disp_pow = disp * Eigen::VectorXd::Ones(br_poly.size()).transpose();
					for (int i = 0; i < disp_pow.cols(); ++i)
						disp_pow.col(i) = disp_pow.array().col(i).pow(i);
					Eigen::VectorXd force = disp_pow * br_poly;

					matlab.plot(disp, force * 2, "", cmap.row(th(j) < th(1 - j) ? cmap.rows() - 1 : 1), 2);
				}

				matlab.eval("ax = gca; ax.FontSize = 8; set(gca,'TickLabelInterpreter', 'latex');");

				stringstream ss_title;
				ss_title << setprecision(2) << "title('$h_1=" << th(0) << "$, " << "$h_2 = " << th(1) << "$');";
				matlab.eval(ss_title.str());

				stringstream ss_size;
				double width = 6;
				double height = 5;
				ss_size << "set(gcf, 'Units', 'centimeters', 'Position', [0, 0, " << width << ", " << height << "], 'PaperUnits', 'centimeters', 'PaperSize', [" << width << ", " << height << "]);";
				matlab.eval(ss_size.str());

			}

		}
	}
	else
	{
		ui::Text("Waiting for Matlab engine");
	}
#else
	ui::Text("No Matlab version");
#endif

	ui::End();
}

void SmpApp::GuiDataExportWindow()
{
	static bool popen = true;
	if (!popen)
	{
		mShowDataExportWindow = false;
		popen = true;
	}
	if (!ui::Begin("Data export", &popen, ImGuiWindowFlags_NoCollapse)) { ui::End(); return; }

	if (mSelectedSimHandler)
	{
		static int u_anchor_base = -1;
		ui::InputInt("Anchor base", &u_anchor_base);
		u_anchor_base = min(max(u_anchor_base, -1), (int)mSelectedSimHandler->simulator->get_smp()->F.rows());

		static bool filter_nomark = false;
		ui::Checkbox("Filter bases without markers", &filter_nomark);

		static bool output_edges = false;
		if (!filter_nomark) ui::Checkbox("Output edges", &output_edges);

		static int u_side = 0;
		ui::Combo("", &u_side, "front\0back\0\0");
		ui::SameLine();
		if (ui::Button("Export current marker locations"))
		{
			auto finfo = mSelectedSimHandler->simulator->finfo()[mCurrentFrame];
			Eigen::MatrixXd markers = finfo->body_centroids;
			for (int i = 0; i < markers.rows(); ++i)
			{
				Eigen::Vector3d axis{ 0, 0, (u_side == 0 ? 1.0 : -1.0) };
				axis = finfo->body_rotation[i] * axis;
				markers.row(i) += axis * (mSelectedSimHandler->simulator->get_smp()->options.base_th + mSelectedSimHandler->simulator->get_smp()->options.base_p);
			}

			if (-1 < u_anchor_base)
			{
				Eigen::Matrix<double, 1, 3> T = -mSelectedSimHandler->simulator->finfo()[mCurrentFrame]->body_centroids.row(u_anchor_base);
				Eigen::Matrix3d R = mSelectedSimHandler->simulator->finfo()[mCurrentFrame]->body_rotation[u_anchor_base];
				//markers += Eigen::VectorXd::Ones(markers.rows()) * T;
				markers.rowwise() += T;
				markers *= R;
			}

			Eigen::MatrixXd markers_filtered(mSelectedSimHandler->simulator->get_smp()->F.rows() - mSelectedSimHandler->simulator->get_smp()->base_align.size(), 3);
			int ind = 0;
			for (int i = 0; i < mSelectedSimHandler->simulator->get_smp()->F.rows(); ++i)
				if (mSelectedSimHandler->simulator->get_smp()->base_align.find(i) == mSelectedSimHandler->simulator->get_smp()->base_align.end())
					markers_filtered.row(ind++) = markers.row(i);
			
			fs::path path = app::getSaveFilePath(mSelectedSimHandler->simulator->get_smp()->working_dir + "/" + mSelectedSimHandler->simulator->get_name() + "_markers.obj",
				vector<string>{"obj"});
			if (!path.empty())
			{
				Eigen::IOFormat fmt_obj_v(-1, 0, " ", "\n", "v ", "", "", "\n");
				Eigen::IOFormat fmt_obj_l(-1, 0, " ", "\n", "l ", "", "", "\n");
				ofstream ofs(path.string());
				ofs << (filter_nomark ? markers_filtered.format(fmt_obj_v) : markers.format(fmt_obj_v));
				if (output_edges) ofs << (mSelectedSimHandler->simulator->get_smp()->FF.array() + 1).format(fmt_obj_l);
			}
		}
	}
	else
	{
		ui::Text("No simulation");
	}

	ui::End();
}

void SmpApp::GuiDebugWindow()
{
	static bool popen = true;
	if (!popen)
	{
		mShowDebugWindow = false;
		popen = true;
	}
	if (!ui::Begin("Debug window", &popen, ImGuiWindowFlags_NoCollapse)) { ui::End(); return; }

	//mSimulationHandlers; // list of all simulations
	//mSelectedSimHandler; // current selected simulation
	//mSelectedSimHandler->simulator->finfo(); // all completed time steps (zero index is fabricated state)
	//mSelectedSimHandler->simulator->psim(); // internal structure interfacing physim
	//mSelectedSimHandler->simulator->psim()->get_current_solver_info(); // current solver iterations info which is not yet in frame info

	if (mSelectedSimHandler)
	{
		ui::PushItemWidth(100);
		static int num_variables;
		num_variables = mSelectedSimHandler->simulator->psim()->get_problem()->GetNumVariables();
		ui::InputInt("Number of variables", &num_variables, 0, 0, ImGuiInputTextFlags_ReadOnly);
		static int num_base_dofs;
		num_base_dofs = mSelectedSimHandler->simulator->get_smp()->F.rows() * 6;
		ui::InputInt("Number of base dofs", &num_base_dofs, 0, 0, ImGuiInputTextFlags_ReadOnly);
		static int num_membrane_dofs;
		num_membrane_dofs = num_variables - num_base_dofs;
		ui::InputInt("Number of membrane dofs", &num_membrane_dofs, 0, 0, ImGuiInputTextFlags_ReadOnly);
		ui::PopItemWidth();
	}

	ui::End();
}
