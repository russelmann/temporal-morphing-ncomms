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

void SmpApp::PopupClearAll()
{
	static string popup_name = "Clear all simulations";
	if (modal_clear_all_request) ui::OpenPopup(popup_name.c_str());
	if (ui::BeginPopupModal(popup_name.c_str(), nullptr, ImGuiWindowFlags_AlwaysAutoResize))
	{
		ui::Text("Are you sure want to remove all simulations?");
		if (ui::Button("Remove all"))
		{
			remove_all();
			ui::CloseCurrentPopup();
		}
		ui::SameLine();
		if (ui::Button("Cancel")) ui::CloseCurrentPopup();
		ui::EndPopup();
		modal_clear_all_request = false;
	}
}

void SmpApp::PopupNewSingle()
{
	static string popup_name = "New single actuator";
	if (modal_new_single_request) ui::OpenPopup(popup_name.c_str());
	if (ui::BeginPopupModal(popup_name.c_str(), nullptr, ImGuiWindowFlags_AlwaysAutoResize))
	{
		static shared_ptr<SmpSettings> u_settings;
		if (modal_new_single_request)
		{
			u_settings = make_shared<SmpSettings>();
			u_settings->mem_model = StaticSimulatorStruct::MEM_LINEAR;
			u_settings->constant_force = 5;
			u_settings->mimic_bounary_weakening = false;
		}

		static float u_length = 8;
		static float u_thickness = 0.4;
		ui::Text("Actuator parameters");
		ui::InputFloat("Actuator length", &u_length);
		ui::InputFloat("Bracket thickness", &u_thickness);

		showGuiSettings(*u_settings);

		if (ui::Button("Ok"))
		{
			new_single(u_length, u_thickness, *u_settings);
			ui::CloseCurrentPopup();
		}
		ui::SameLine();
		if (ui::Button("Cancel")) ui::CloseCurrentPopup();
		ui::EndPopup();
		modal_new_single_request = false;
	}
}

void SmpApp::PopupNewRegular()
{
	static string popup_name = "New regular structure";
	if (modal_new_regular_request) ui::OpenPopup(popup_name.c_str());
	if (ui::BeginPopupModal(popup_name.c_str(), nullptr, ImGuiWindowFlags_AlwaysAutoResize))
	{
		static SmpSettings settings;
		if (modal_new_regular_request)
		{
			settings.mem_model = StaticSimulatorStruct::MEM_FINE;
			settings.triangle = "a5qQ";
		}

		static int u_rows = 7;
		static int u_cols = 4;
		static float u_length = 8;
		static float u_thickness = 0.4;
		ui::Text("Size of the structure");
		ui::InputInt("Rows", &u_rows); u_rows = max(1, u_rows);
		ui::InputInt("Columns", &u_cols); u_cols = max(1, u_cols);
		ui::Spacing();
		ui::InputFloat("Actuator length", &u_length);
		ui::InputFloat("Bracket thickness", &u_thickness);

		showGuiSettings(settings);

		if (ui::Button("Ok"))
		{
			new_regular(u_rows, u_cols, u_length, u_thickness, settings);
			ui::CloseCurrentPopup();
		}
		ui::SameLine();
		if (ui::Button("Cancel")) ui::CloseCurrentPopup();
		ui::EndPopup();
		modal_new_regular_request = false;
	}
}

void SmpApp::PopupOpenStencil()
{
	static string popup_name = "Open structure stencil";
	if (model_open_stencil_request) ui::OpenPopup(popup_name.c_str());
	if (ui::BeginPopupModal(popup_name.c_str(), nullptr, ImGuiWindowFlags_AlwaysAutoResize))
	{
		static SmpSettings settings;
		if (model_open_stencil_request)
		{
			//default modifications to settings
		}

		ui::InputText("", &mCurrentFile);
		ui::SameLine();
		if (ui::Button("Select")) read_file_info();

		showGuiSettings(settings);

		if (ui::Button("Ok"))
		{
			open_current_file(settings);
			ui::CloseCurrentPopup();
		}
		ui::SameLine();
		if (ui::Button("Cancel")) ui::CloseCurrentPopup();
		ui::EndPopup();
		model_open_stencil_request = false;
	}
}

void SmpApp::PopupCloseApp()
{
	static string popup_name = "Close app";
	if (mExitRequest) ui::OpenPopup(popup_name.c_str());
	if (ui::BeginPopupModal(popup_name.c_str(), nullptr, ImGuiWindowFlags_AlwaysAutoResize))
	{
		ui::Text("Are you sure want to exit?");
		if (ui::Button("Exit app")) getWindow()->close();
		ui::SameLine();
		if (ui::Button("Cancel")) ui::CloseCurrentPopup();
		ui::EndPopup();
		mExitRequest = false;
	}
}

void SmpApp::PopupBuildMesh()
{
	static string popup_name = "Build mesh";
	if (modal_build_mesh_request) ui::OpenPopup(popup_name.c_str());
	if (ui::BeginPopupModal(popup_name.c_str(), nullptr, ImGuiWindowFlags_AlwaysAutoResize))
	{
		static int u_build_mode = 0;
		ui::Combo("Build state", &u_build_mode, "Flat\0Actuated\0");

		static SmpPrintOptions u_print_options;
		if (modal_build_mesh_request) u_print_options = mSelectedSimHandler->simulator->get_smp()->options.default_print_options;

		if (u_build_mode == 0)
		{
			ui::PushItemWidth(80);
			ui::Text("NB: heightfield enforced!");
			ui::Checkbox("Markers", &u_print_options.markers);
			ui::InputInt("Handle", &u_print_options.handle);
			u_print_options.handle = max(u_print_options.handle, -1);
			u_print_options.handle = min(u_print_options.handle, (int)mSelectedSimHandler->simulator->get_smp()->F.rows() - 1);
			ui::SameLine();
			if (ui::Button("none")) u_print_options.handle = -1;
			ui::Checkbox("Alignment holes", &u_print_options.alignment);
			ui::InputFloat("Hole diameter", &u_print_options.align_diam);
			u_print_options.align_diam = max(u_print_options.align_diam, 0.0f);
			ui::PopItemWidth();
		}

		static vector<Eigen::VectorXf> u_base_done(2, Eigen::VectorXf::Zero(mSelectedSimHandler->simulator->get_smp()->F.rows()));
		ui::PlotHistogram("Bases built", u_base_done[u_build_mode].data(), u_base_done[u_build_mode].size(), 0, nullptr, 0, 2);

		bool done = u_base_done[u_build_mode].sum() == mSelectedSimHandler->simulator->get_smp()->F.rows() * 2;
		if (mSelectedSimHandler->simulator->get_smp()->is_build_in_progress()) ui::Text("Mesh in progress...");
		else if (done) ui::Text("Mesh is done");
		else ui::Text("No mesh");

		if (mSelectedSimHandler->simulator->get_smp()->is_build_in_progress())
		{
			if (ui::Button("Cancel"))
			{
				mSelectedSimHandler->simulator->get_smp()->stop_build();
				u_base_done[u_build_mode].setZero();
			}
		}
		else
		{
			if (ui::Button("Start"))
			{
				u_base_done[u_build_mode] = Eigen::VectorXf::Zero(mSelectedSimHandler->simulator->get_smp()->F.rows());
				string name = mSelectedSimHandler->simulator->get_smp()->name;
				mSelectedSimHandler->simulator->get_smp()->build_scad(name, (u_build_mode == 0 ? SMP_FLT : SMP_ACT), &u_print_options, u_base_done[u_build_mode].data());
			}
			ui::SameLine();
			if (ui::Button("Open folder")) os_reveal_folder(mSelectedSimHandler->simulator->get_smp()->working_dir);
		}

		ui::Separator();

		if (ui::Button("Build support frame"))
		{
			mSelectedSimHandler->simulator->get_smp()->build_support_frame(0);
			mSelectedSimHandler->simulator->get_smp()->build_support_frame(1);
		}

		if (ui::Button("Generate scad structures"))
		{
			string name = mSelectedSimHandler->simulator->get_smp()->name;
			mSelectedSimHandler->simulator->get_smp()->export_scad(name, &u_print_options);
		}

		if (ui::Button("Save parameterization"))
		{
			string name = mSelectedSimHandler->simulator->get_smp()->name;
			auto& uv = mSelectedSimHandler->simulator->get_smp()->uv;
			Eigen::MatrixXd V(uv.rows(), 3);
			V.leftCols(2) = uv;
			V.col(2).setZero();
			igl::writeOBJ(name + "_map3d.obj", mSelectedSimHandler->simulator->get_smp()->V, mSelectedSimHandler->simulator->get_smp()->F);
			igl::writeOBJ(name + "_map2d.obj", V, mSelectedSimHandler->simulator->get_smp()->F);
			ofstream(name + "_scaling.txt") << mSelectedSimHandler->simulator->get_smp()->eu;
		}

		ui::Separator();

		if (ui::Button("Close"))
		{
			if (!mSelectedSimHandler->simulator->get_smp()->is_build_in_progress())
			{
				u_base_done[0].setZero();
				u_base_done[1].setZero();
				ui::CloseCurrentPopup();
			}
		}
		ui::EndPopup();
		modal_build_mesh_request = false;
	}
}

void SmpApp::PopupBuildVideo()
{
	static string popup_name = "Build video";
	if (modal_build_video_request) ui::OpenPopup(popup_name.c_str());
	if (ui::BeginPopupModal(popup_name.c_str(), nullptr, ImGuiWindowFlags_AlwaysAutoResize))
	{
		static int u_anchor_base = 0;
		ui::InputInt("Anchor base", &u_anchor_base);
		u_anchor_base = min(max(u_anchor_base, 0), (int)mSelectedSimHandler->simulator->get_smp()->F.rows());

		const auto& base_align = mSelectedSimHandler->simulator->get_smp()->base_align;
		static stringstream ss;
		ss.str(std::string());
		ss << *base_align.begin();
		auto b_align = base_align.begin();
		for (++b_align; b_align != base_align.end(); ++b_align)
			ss << ", " << *b_align;
		ui::Text("Align: %s", ss.str().c_str());

		if (ui::Button("Export mesh frames"))
		{
			auto dir = fs::path(mSelectedSimHandler->simulator->get_smp()->working_dir).append("frames");
			string dir_name = dir.string();
			fs::create_directory(dir_name);

			int nframes = mSelectedSimHandler->simulator->get_simulated_frames();
			for (int i = 0; i < nframes; ++i)
			{
				stringstream ssuffix;
				ssuffix << setw(4) << setfill('0') << i;
				string suffix = ssuffix.str();

				string fname_solid = "frame_solid_" + suffix + ".obj";
				string fname_brackets = "frame_brackets_" + suffix + ".obj";
				string fname_membrane = "frame_membrane_" + suffix + ".obj";

				Eigen::Matrix<float,1,3> T = -mSelectedSimHandler->simulator->finfo()[i]->body_centroids.row(u_anchor_base).cast<float>();
				Eigen::Matrix3f R = mSelectedSimHandler->simulator->finfo()[i]->body_rotation[u_anchor_base].cast<float>();

				igl::writeOBJ((dir / fname_solid).string(),
					(mSelectedSimHandler->simulator->get_body_vertices_at_frame(i, 0).rowwise() + T) * R,
					mSelectedSimHandler->simulator->get_body_faces());

				igl::writeOBJ((dir / fname_brackets).string(),
					(mSelectedSimHandler->simulator->get_bracket_vertices_at_frame(i, 0).rowwise() + T) * R,
					mSelectedSimHandler->simulator->get_bracket_faces());

				if (mSelectedSimHandler->simulator->get_membrane_model() == StaticSimulatorStruct::MEM_FINE)
				{
					igl::writeOBJ((dir / fname_membrane).string(),
						(mSelectedSimHandler->simulator->get_positions_at_frame(i, 0).cast<float>().rowwise() + T) * R,
						mSelectedSimHandler->simulator->psim()->memtri);
				}
			}
		}

#if false
		if (ui::Button("Export current frame"))
		{
			auto dir = mCoreSettings->get_smpup_path();
			dir.append("mitsuba");
			string dir_name = dir.string();
			stringstream sprefix;
			sprefix << "frame_" << setw(4) << setfill('0') << mCurrentFrame << "_" << setw(4) << setfill('0') << int(mCurrentFrameAlpha * 1000) << "_";
			string prefix = sprefix.str();

			string fname_solid = prefix + "solid.obj";
			string fname_brackets = prefix + "brackets.obj";
			string fname_membrane = prefix + "membrane.obj";

			igl::writeOBJ(dir_name + "/" + fname_solid,
				mSelectedSimHandler->simulator->get_body_vertices_at_frame(mCurrentFrame, mCurrentFrameAlpha),
				mSelectedSimHandler->simulator->get_body_faces());

			igl::writeOBJ(dir_name + "/" + fname_brackets,
				mSelectedSimHandler->simulator->get_bracket_vertices_at_frame(mCurrentFrame, mCurrentFrameAlpha),
				mSelectedSimHandler->simulator->get_bracket_faces());

			if (mSelectedSimHandler->simulator->get_membrane_model() == StaticSimulatorStruct::MEM_FINE)
			{
				igl::writeOBJ(dir_name + "/" + fname_membrane,
					mSelectedSimHandler->simulator->get_positions_at_frame(mCurrentFrame, mCurrentFrameAlpha),
					mSelectedSimHandler->simulator->psim()->memtri);
			}

			XmlTree scene = XmlTree(loadFile("C:/Research/smpup/source/smpsym/mitsuba/template.xml")).getChild("scene");

			XmlTree transform("transform", "");
			transform.setAttribute("name", "toWorld");
			transform.push_back(XmlTree("rotate", "").setAttribute("x", 1).setAttribute("angle", 180));
			transform.push_back(XmlTree("rotate", "").setAttribute("y", 1).setAttribute("angle", 200));
			transform.push_back(XmlTree("translate", "").setAttribute("x", 3000).setAttribute("y", 1000).setAttribute("z", 900));

			XmlTree bsdf("bsdf", "");
			bsdf.setAttribute("type", "plastic");
			bsdf.push_back(XmlTree("srgb", "").setAttribute("name", "diffuseReflectance").setAttribute("value", "#FFFFFF"));
			bsdf.push_back(XmlTree("float", "").setAttribute("name", "intIOR").setAttribute("value", 1.1));

			// Solid
			XmlTree shape_solid("shape", "");
			shape_solid.setAttribute("type", "obj");
			shape_solid.push_back(XmlTree("string", "").setAttribute("name", "filename").setAttribute("value", fname_solid));
			shape_solid.push_back(XmlTree("float", "").setAttribute("name", "maxSmoothAngle").setAttribute("value", 0));
			shape_solid.push_back(transform);
			shape_solid.push_back(bsdf);
			scene.push_back(shape_solid);

			// Brackets
			XmlTree shape_brackets("shape", "");
			shape_brackets.setAttribute("type", "obj");
			shape_brackets.push_back(XmlTree("string", "").setAttribute("name", "filename").setAttribute("value", fname_brackets));
			shape_brackets.push_back(XmlTree("float", "").setAttribute("name", "maxSmoothAngle").setAttribute("value", 0));
			shape_brackets.push_back(transform);
			shape_brackets.push_back(bsdf);
			scene.push_back(shape_brackets);

			// Membrane
			if (mSelectedSimHandler->simulator->get_membrane_model() == StaticSimulatorStruct::MEM_FINE)
			{
				XmlTree shape_membrane("shape", "");
				shape_membrane.setAttribute("type", "obj");
				shape_membrane.push_back(XmlTree("string", "").setAttribute("name", "filename").setAttribute("value", fname_membrane));
				shape_membrane.push_back(transform);
				XmlTree bsdf_ts("bsdf", "");
				bsdf_ts.setAttribute("type", "twosided");
				bsdf.getChild("srgb").setAttribute("value", "1.0, 1.0, 0.8");
				bsdf_ts.push_back(bsdf);
				shape_membrane.push_back(bsdf_ts);
				scene.push_back(shape_membrane);
			}

			scene.write(writeFile(dir_name + "/" + prefix + ".xml"));
		}

		if (ui::Button("Export frames"))
		{
			auto dir = mCoreSettings->get_smpup_path();
			dir.append("mitsuba");
			string dir_name = dir.string();

			ofstream bat_file(dir_name + "/all_render.bat");
			stringstream tonemap;
			tonemap << "mtsutil tonemap ";
			int nframes = mSelectedSimHandler->simulator->get_simulated_frames();
			for (int i = 0; i < nframes; ++i)
			{
				stringstream sprefix;
				sprefix << "frame_" << setw(4) << setfill('0') << i << "_";
				string prefix = sprefix.str();

				string fname_solid = prefix + "solid.obj";
				string fname_brackets = prefix + "brackets.obj";
				string fname_membrane = prefix + "membrane.obj";

				igl::writeOBJ(dir_name + "/" + fname_solid,
					mSelectedSimHandler->simulator->get_body_vertices_at_frame(i, 0),
					mSelectedSimHandler->simulator->get_body_faces());

				igl::writeOBJ(dir_name + "/" + fname_brackets,
					mSelectedSimHandler->simulator->get_bracket_vertices_at_frame(i, 0),
					mSelectedSimHandler->simulator->get_bracket_faces());

				if (mSelectedSimHandler->simulator->get_membrane_model() == StaticSimulatorStruct::MEM_FINE)
				{
					igl::writeOBJ(dir_name + "/" + fname_membrane,
						mSelectedSimHandler->simulator->get_positions_at_frame(i, 0),
						mSelectedSimHandler->simulator->psim()->memtri);
				}

				XmlTree scene = XmlTree(loadFile("C:/Research/smpup/source/smpsym/mitsuba/template.xml")).getChild("scene");

				XmlTree transform("transform", "");
				transform.setAttribute("name", "toWorld");
				transform.push_back(XmlTree("rotate", "").setAttribute("x", 1).setAttribute("angle", 180));
				transform.push_back(XmlTree("rotate", "").setAttribute("y", 1).setAttribute("angle", 200));
				transform.push_back(XmlTree("translate", "").setAttribute("x", 3000).setAttribute("y", 1000).setAttribute("z", 900));

				XmlTree bsdf("bsdf", "");
				bsdf.setAttribute("type", "plastic");
				bsdf.push_back(XmlTree("srgb", "").setAttribute("name", "diffuseReflectance").setAttribute("value", "#FFFFFF"));
				bsdf.push_back(XmlTree("float", "").setAttribute("name", "intIOR").setAttribute("value", 1.1));

				// Solid
				XmlTree shape_solid("shape", "");
				shape_solid.setAttribute("type", "obj");
				shape_solid.push_back(XmlTree("string", "").setAttribute("name", "filename").setAttribute("value", fname_solid));
				shape_solid.push_back(XmlTree("float", "").setAttribute("name", "maxSmoothAngle").setAttribute("value", 0));
				shape_solid.push_back(transform);
				shape_solid.push_back(bsdf);
				scene.push_back(shape_solid);

				// Brackets
				XmlTree shape_brackets("shape", "");
				shape_brackets.setAttribute("type", "obj");
				shape_brackets.push_back(XmlTree("string", "").setAttribute("name", "filename").setAttribute("value", fname_brackets));
				shape_brackets.push_back(XmlTree("float", "").setAttribute("name", "maxSmoothAngle").setAttribute("value", 0));
				shape_brackets.push_back(transform);
				shape_brackets.push_back(bsdf);
				scene.push_back(shape_brackets);

				// Membrane
				if (mSelectedSimHandler->simulator->get_membrane_model() == StaticSimulatorStruct::MEM_FINE)
				{
					XmlTree shape_membrane("shape", "");
					shape_membrane.setAttribute("type", "obj");
					shape_membrane.push_back(XmlTree("string", "").setAttribute("name", "filename").setAttribute("value", fname_membrane));
					shape_membrane.push_back(transform);
					XmlTree bsdf_ts("bsdf", "");
					bsdf_ts.setAttribute("type", "twosided");
					bsdf.getChild("srgb").setAttribute("value", "1.0, 1.0, 0.8");
					bsdf_ts.push_back(bsdf);
					shape_membrane.push_back(bsdf_ts);
					scene.push_back(shape_membrane);
				}

				scene.write(writeFile(dir_name + "/" + prefix + ".xml"));
				bat_file << "mitsuba " << prefix << ".xml\n";
				tonemap << prefix << ".exr ";
			}
			bat_file << tonemap.str() << endl;
			//bat_file << "ffmpeg -f image2 -i \"frame_%04d_.png\" animation.mp4" << endl;
		}
#endif 

		if (ui::Button("Close")) ui::CloseCurrentPopup();

		ui::EndPopup();
		modal_build_video_request = false;
	}
}

void SmpApp::PopupCoreSettings()
{
	static string popup_name = "Core settings";
	if (modal_core_settings_request) ui::OpenPopup(popup_name.c_str());
	if (ui::BeginPopupModal(popup_name.c_str(), nullptr, ImGuiWindowFlags_AlwaysAutoResize))
	{
		static CoreSettings u_core_settings;
		if (modal_core_settings_request) u_core_settings = *mCoreSettings;
		ui::InputText("OpenSCAD path", &u_core_settings.openscad_path);
		ui::InputText("SmpUp core path", &u_core_settings.smpup_path);
		ui::InputFloat("Time step, sec", &u_core_settings.timestep); tooltip("Restart SmpApp to apply change");
		if (ui::Button("Restore defaults")) u_core_settings = CoreSettings();
		if (ui::Button("Ok"))
		{
			*mCoreSettings = u_core_settings;
			mCoreSettings->save();
			ui::CloseCurrentPopup();
		}
		ui::SameLine();
		if (ui::Button("Cancel")) ui::CloseCurrentPopup();
		ui::EndPopup();
		modal_core_settings_request = false;
	}
}

void SmpApp::PopupMembraneLinearity()
{
	static string popup_name = "Membrane linearity sampling";
	if (modal_membrane_linearity_request) ui::OpenPopup(popup_name.c_str());
	if (ui::BeginPopupModal(popup_name.c_str(), nullptr, ImGuiWindowFlags_AlwaysAutoResize))
	{
		ui::Text("Not tested");
		if (ui::Button("Start"))
		{
			SmpSettings settings;
			settings.mem_model = StaticSimulatorStruct::MEM_FINE;
			settings.timestep_mul = 1;

			Eigen::VectorXd act_lengths = natural_seriesd(5, 9);
			SmpOptions options = mSmpOptions;
			int timesteps = 120; // 120;
			Eigen::VectorXd leng(act_lengths.size());
			Eigen::MatrixXd disp = Eigen::MatrixXd::Zero(timesteps - 1, act_lengths.size());
			Eigen::MatrixXd forc = Eigen::MatrixXd::Zero(timesteps - 1, act_lengths.size());

			for (int i = 0; i < act_lengths.size(); ++i)
			{
				options.min_gap = act_lengths(i) / 3; // assuming 1.5 contraction
				options.min_bumper = act_lengths(i) / 3 - mSmpOptions.base_r;

				Smp smp;
				smp.options = options;
				smp.init_regular(7, 4);
				int actuator_id = 65; // NOTE: hard-coded 65 for 7x4 regular
				smp.spth.setConstant(0.3); // NOTE: hard-coded 0.3

				SmpSimulator* simulator = new SmpSimulator(0, timesteps, 1, settings, mLinearMemModel, mBracket);
				simulator->load_smp(smp);
				simulator->simulate();
				while (simulator->get_simulated_frames() <= timesteps && simulator->is_running_task())
				{
					//if (!simulator->finfo().empty() && simulator->finfo().back()->collided(simulator->get_actuators()[actuator_id]->segment_springs(0)))
					//{
					//	simulator->stop_simulation();
					//	clog << "Stop!!\n";
					//}
					//clog << "Frames: " << simulator->get_simulated_frames() << "\n";
					process_sleep(10);
				}
				clog << "Hop!!\n";

				for (int frame_id = 0; frame_id < simulator->get_simulated_frames(); ++frame_id)
				{
					Eigen::Vector4d bracket_length;
					Eigen::Vector4d bracket_forces;
					simulator->get_bracket_length_at_frame(actuator_id, frame_id, bracket_length);
					simulator->get_bracket_forces_at_frame(actuator_id, frame_id, bracket_forces);
					if (simulator->finfo()[frame_id]->collided(simulator->get_actuators()[actuator_id]->segment_springs(0))) continue; // checking only one spring out of four
					if (frame_id == 0) leng(i) = bracket_length.mean();
					else
					{
						disp(frame_id - 1, i) = bracket_length.mean();
						forc(frame_id - 1, i) = bracket_forces.sum();
					}
				}
			}
			Eigen::MatrixXd M_leng = Eigen::VectorXd::Ones(disp.rows()) * leng.transpose();
			auto var_leng = Eigen::Map<Eigen::ArrayXd>(M_leng.data(), M_leng.size());
			auto var_disp = Eigen::Map<Eigen::ArrayXd>(disp.data(), disp.size());
			auto var_forc = Eigen::Map<Eigen::ArrayXd>(forc.data(), forc.size());

			Eigen::MatrixXd X = Eigen::MatrixXd::Ones(var_forc.size(), 9);
			X.col(1) = var_leng;
			X.col(2) = var_leng.array().square();
			X.block(0, 3, X.rows(), 3) = X.leftCols(3);
			X.block(0, 6, X.rows(), 3) = X.leftCols(3);
			X.col(3).array() *= var_disp;
			X.col(4).array() *= var_disp;
			X.col(5).array() *= var_disp;
			X.col(6).array() *= var_disp.array().square();
			X.col(7).array() *= var_disp.array().square();
			X.col(8).array() *= var_disp.array().square();

			int n_keep = (var_forc != 0).count();
			Eigen::VectorXd y = var_forc;
			Eigen::MatrixXd X_keep(n_keep, X.cols());
			Eigen::VectorXd y_keep(n_keep);
			int ind = 0;
			for (int i = 0; i < X.rows(); ++i)
			{
				if (var_forc(i) != 0)
				{
					X_keep.row(ind) = X.row(i);
					y_keep(ind) = y(i);
					++ind;
				}
			}
			X = X_keep;
			y = y_keep;
			Eigen::VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
			ofstream("beta.txt") << beta;
		}
		ui::SameLine();
		if (ui::Button("Close")) ui::CloseCurrentPopup();
		ui::EndPopup();
		modal_membrane_linearity_request = false;
	}
}

void SmpApp::PopupBuildSpecimenBlock()
{
	static string popup_name = "Build specimen block";
	if (modal_build_specimen_block_request) ui::OpenPopup(popup_name.c_str());
	if (ui::BeginPopupModal(popup_name.c_str(), nullptr, ImGuiWindowFlags_AlwaysAutoResize))
	{
		if (ui::Button("Start")) specimen_block(mCoreSettings->openscad_path, mCoreSettings->smpup_path);
		ui::SameLine();
		if (ui::Button("Cancel")) ui::CloseCurrentPopup();
		ui::EndPopup();
		modal_build_specimen_block_request = false;
	}
}

void SmpApp::PopupBracketSampling()
{
	static string popup_name = "Bracket sampling";
	if (modal_bracket_sampling_request) ui::OpenPopup(popup_name.c_str());
	if (ui::BeginPopupModal(popup_name.c_str(), nullptr))//, ImGuiWindowFlags_AlwaysAutoResize))
	{
		ui::Text("Current bracket configurer");
		ui::BeginChild("scroll_brackets", ImVec2(0, 200), true, ImGuiWindowFlags_HorizontalScrollbar);
		for (int length_id = 0; length_id < mBracketConfigurer->lengths.size(); ++length_id)
		{
			if (ui::TreeNode((string("legnth ") + to_string(mBracketConfigurer->lengths(length_id))).c_str()))
			{
				stringstream ss("thicknesses: ");
				for (int thickness_id = 0; thickness_id < mBracketConfigurer->thicknesses[length_id].size(); ++thickness_id)
				{
					if (0 < thickness_id) ss << "; ";
					ss << mBracketConfigurer->thicknesses[length_id](thickness_id);
				}
				ui::Text("%s", ss.str().c_str());
				ui::TreePop();
			}
		}
		ui::EndChild();

#ifdef USE_MATLAB
		if (matlab_started && ui::Button("Plot everything"))
		{
			matlab.figure("single_sampling", "Single actuator lengths", "time, sec", "length, mm");

			Eigen::VectorXd times(121);
			for (int i = 0; i < 121; ++i)
				times(i) = i;

			for (int length_id = 0; length_id < mBracketConfigurer->lengths.size(); ++length_id)
				for (int thickness_id = 0; thickness_id < mBracketConfigurer->thicknesses[length_id].size(); ++thickness_id)
					matlab.plot(times, mBracketConfigurer->displacements[length_id].col(thickness_id), "blue");
		}
#endif

		static shared_ptr<BracketConfigurer> u_bracket_configurer;
		if (modal_bracket_sampling_request)
		{
			u_bracket_configurer = make_shared<BracketConfigurer>(mLinearMemModel, mBracket);
			u_bracket_configurer->set_default_directory(mBracketConfigurer->get_default_directory());
		}

		ui::Text("Single bracket sampling");
		ui::Text("Lengths");
		ui::ProgressBar(u_bracket_configurer->get_progress_length());
		ui::Text("Thicknesses");
		ui::ProgressBar(u_bracket_configurer->get_progress_thickness());

		if (u_bracket_configurer->valid_computation_done())
		{
			if (ui::Button("Save and close"))
			{
				u_bracket_configurer->save();
				mBracketConfigurer = u_bracket_configurer;
				u_bracket_configurer.reset();
				ui::CloseCurrentPopup();
			}
		}
		else
		{
			if (ui::Button("Start"))
			{
				SmpSettings settings;
				settings.mem_model = StaticSimulatorStruct::MEM_LINEAR;

				Eigen::VectorXd lengths = 5 + 4 * natural_seriesd(0, 20) / 20;
				Eigen::VectorXd thickness = 0.3 + 0.35 * natural_seriesd(0, 20) / 20;
				vector<Eigen::VectorXd> thicknesses(lengths.size(), thickness);

				SmpOptions options = mSmpOptions;

				u_bracket_configurer->timestep = 1;
				u_bracket_configurer->lengths = lengths;
				u_bracket_configurer->thicknesses = thicknesses;

				u_bracket_configurer->single_actuator_sampling(options);
			}
			ui::SameLine();
			if (ui::Button("Cancel"))
			{
				u_bracket_configurer.reset();
				ui::CloseCurrentPopup();
			}
		}
		ui::EndPopup();
		modal_bracket_sampling_request = false;
	}
}